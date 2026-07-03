cf_density_ratios <- function(task, learners, mtp, control, pb) {
  ans <- vector("list", length = length(task$folds))

  if (length(learners) == 1 && learners == "mean") {
    warning("Using 'mean' as the only learner of the density ratios will always result in a misspecified model! If your exposure is randomized, consider using `c('glm', 'cv_glmnet')`.",
            call. = FALSE)
  }

  for (fold in seq_along(task$folds)) {
    ans[[fold]] <- future::future({
      estimate_density_ratios(task, fold, learners, mtp, control, pb)
    },
    seed = TRUE)
  }

  ans <- future::value(ans)

  ans <- list(density_ratios = recombine(rbind_depth(ans, "ratios"), task$folds),
              fits = lapply(ans, function(x) x[["fits"]]))

  ans$density_ratios <- trim(ans$density_ratios, control$.trim)
  ans
}

estimate_density_ratios <- function(task, fold, learners, mtp, control, pb) {
  natural <- get_folded_data(task$natural, task$folds, fold)
  shifted <- get_folded_data(task$shifted, task$folds, fold)

  density_ratios <- matrix(nrow = nrow(natural$valid), ncol = task$time_horizon)
  fits <- vector("list", length = task$time_horizon)

  for (time in seq_len(task$time_horizon)) {
    i_train <- task$observed(natural$train, time - 1) %and% task$is_at_risk(natural$train, time)
    i_valid <- task$observed(natural$valid, time - 1) %and% task$is_at_risk(natural$valid, time)

    A_t <- current_trt(task$vars$A, time)
    V_through_t <- task$vars$cluster_cols_through(time)        # current + past cluster exposure
    V_t <- if (is.null(task$vars$V)) character(0) else current_trt(task$vars$V, time)

    # ---- Unit-level treatment density ratio: g*(A_t | V_t, H) / g(A_t | V_t, H) ----
    # The cluster exposure (current + past) enters as a *covariate* here, so the
    # unit model estimates the conditional ratio given V. Only the A columns are
    # shifted in the stacked data, so V is held at its natural value in both halves.
    trt_vars <- c("..i..lmtp_id", task$vars$history("A", time), V_through_t, A_t, "..i..lmtp_stack_indicator")
    trt_vars <- unique(trt_vars)
    stacked <- stack_data(natural$train, shifted$train, task$vars$A, task$vars$C, time)
    fit_trt <- run_ensemble(stacked[rep(i_train, 2), trt_vars], "..i..lmtp_stack_indicator",
                            learners, "binomial", "..i..lmtp_id",
                            control$.learners_trt_folds,
                            control$.discrete,
                            control$.info)

    # Separate censoring model
    if (!is.null(task$vars$C)) {
      cens_vars <- unique(c("..i..lmtp_id", task$vars$history_cens(time), V_through_t, task$vars$C[time]))
      fit_cens <- run_ensemble(natural$train[i_train, cens_vars], task$vars$C[time],
                               learners, "binomial", "..i..lmtp_id",
                               control$.learners_trt_folds,
                               control$.discrete,
                               control$.info)
      pred_cens <- rep(-999L, nrow(natural$valid))
      pred_cens[i_valid] <- predict(fit_cens, natural$valid[i_valid, ])
    } else {
      fit_cens <- NULL
      pred_cens <- 1
    }

    # ---- Cluster-level treatment density ratio: g*(V_t | H_c) / g(V_t | H_c) ----
    cluster_fit <- NULL
    cluster_ratio <- rep(1, nrow(natural$valid))
    if (!is.null(task$vars$V)) {
      # H_c for the V model: explicit cluster set if supplied, else auto-detected.
      if (task$vars$has_cluster_covs()) {
        cluster_cov_pool <- unique(c(task$vars$history_cluster(time),
                                     task$vars$cluster_cols_through(time - 1)))
      } else {
        cluster_cov_pool <- unique(c(task$vars$history("A", time),
                                     task$vars$cluster_cols_through(time - 1)))
      }
      
      cluster_out <- estimate_cluster_ratio(
        natural_train = natural$train, shifted_train = shifted$train,
        natural_valid = natural$valid, shifted_valid = shifted$valid,
        V_t = V_t,
        cov_pool = cluster_cov_pool,
        i_train = i_train, i_valid = i_valid,
        learners = learners, control = control, mtp = mtp
      )
      cluster_ratio <- cluster_out$ratio
      cluster_fit <- cluster_out$fit
    }
    
    if (control$.return_full_fits) {
      fits[[time]] <- list(treatment = fit_trt, censoring = fit_cens, cluster = cluster_fit)
    } else {
      fits[[time]] <- list(treatment = extract_sl_weights(fit_trt),
                           censoring = if (is.null(fit_cens)) NULL else extract_sl_weights(fit_cens),
                           cluster = if (is.null(cluster_fit)) NULL else extract_sl_weights(cluster_fit))
    }

    pred <- rep(-999L, nrow(natural$valid))
    pred[i_valid] <- predict(fit_trt, natural$valid[i_valid, ])

    obs <- task$observed(natural$valid, time)
    at_risk <- task$is_at_risk(natural$valid, time)
    followed <- followed_rule(natural$valid, shifted$valid, A_t, mtp)

    pred <- ifelse(followed & !mtp, pmax(pred, 0.5), pred)
    unit_ratio <- ((pred * obs * at_risk * followed) / (1 - pmin(pred, 0.999))) * (1 / pred_cens)
    
    # Joint ratio for the combined exposure = unit-level (A) x cluster-level (V)
    density_ratios[, time] <- unit_ratio * cluster_ratio

    pb()
  }

  list(ratios = density_ratios, fits = fits)
}

# Cluster-level density ratio for the cluster exposure V:
#   g*(V_t | H_c) / g(V_t | H_c)
# estimated with the classification ("density-ratio") trick on a data set with
# exactly one row per cluster. Because V is constant within cluster, fitting on
# unit rows would raise each cluster's likelihood contribution to the power of
# its size; de-duplicating to the cluster level fixes this. The estimated ratio
# is then mapped back to every unit in the cluster.
estimate_cluster_ratio <- function(natural_train, shifted_train,
                                   natural_valid, shifted_valid,
                                   V_t, cov_pool, i_train, i_valid,
                                   learners, control, mtp) {
  id_col <- "..i..lmtp_id"
  
  # Restrict to observations contributing at this time point
  nt <- natural_train[i_train, , drop = FALSE]
  st <- shifted_train[i_train, , drop = FALSE]
  
  # Cluster-level covariates = covariates in the pool that are constant within
  # cluster (genuinely cluster level). Detected on the training data.
  is_constant <- function(d, cols) {
    vapply(cols, function(v) {
      all(tapply(d[[v]], d[[id_col]], function(x) length(unique(x[!is.na(x)])) <= 1))
    }, logical(1))
  }
  cov_pool <- setdiff(unique(cov_pool), V_t)
  Ccov <- if (length(cov_pool)) cov_pool[is_constant(nt, cov_pool)] else character(0)
  
  if (length(Ccov) == 0) {
    warning("No cluster-constant covariates available for the cluster-exposure model; ",
            "estimating an unconditional (marginal) cluster density ratio. Include ",
            "cluster-level confounders (e.g., cluster summaries) in `baseline`/`time_vary` ",
            "to condition on them.", call. = FALSE)
  }
  
  # Collapse to one row per cluster (V_t and Ccov are constant within cluster).
  first_by_id <- function(d, cols) {
    keep <- !duplicated(d[[id_col]])
    out <- d[keep, c(id_col, cols), drop = FALSE]
    rownames(out) <- NULL
    out
  }
  
  nat_c <- first_by_id(nt, c(Ccov, V_t))
  shi_c <- nat_c
  shi_c[, V_t] <- first_by_id(st, V_t)[, V_t]
  
  stacked <- rbind(nat_c, shi_c)
  stacked[["..i..lmtp_stack_indicator"]] <- rep(c(0, 1), each = nrow(nat_c))
  
  fit <- run_ensemble(stacked, "..i..lmtp_stack_indicator",
                      learners, "binomial", id_col,
                      control$.learners_trt_folds,
                      control$.discrete,
                      control$.info)
  
  # Predict on the validation clusters (one row per cluster), then broadcast.
  nv <- natural_valid[i_valid, , drop = FALSE]
  sv <- shifted_valid[i_valid, , drop = FALSE]
  
  ratio <- rep(1, nrow(natural_valid))
  if (nrow(nv) > 0) {
    valid_c <- first_by_id(nv, c(Ccov, V_t))
    pred_c <- predict(fit, valid_c)
    
    # Did each cluster follow the (deterministic) V rule? Only relevant if !mtp.
    if (!mtp && length(V_t)) {
      v_nat <- first_by_id(nv, V_t)[, V_t, drop = FALSE]
      v_shi <- first_by_id(sv, V_t)[, V_t, drop = FALSE]
      followed_c <- vapply(seq_len(nrow(v_nat)), function(cc) {
        as.integer(all(mapply(function(x, y) isTRUE(all.equal(x, y)),
                              v_nat[cc, ], v_shi[cc, ])))
      }, integer(1))
      pred_c <- ifelse(followed_c == 1, pmax(pred_c, 0.5), pred_c)
    } else {
      followed_c <- 1
    }
    
    ratio_c <- (pred_c * followed_c) / (1 - pmin(pred_c, 0.999))
    
    # Map cluster ratio back onto every validation row by cluster id
    names(ratio_c) <- valid_c[[id_col]]
    ratio[i_valid] <- ratio_c[as.character(natural_valid[[id_col]][i_valid])]
  }
  
  list(ratio = ratio, fit = fit)
}

stack_data <- function(natural, shifted, trt, cens, time) {
  shifted_half <- natural

  if (length(trt) > 1 || time == 1) {
    shifted_half[, trt[[time]]] <- shifted[, trt[[time]]]
  }

  if (!is.null(cens)) {
    shifted_half[[cens[time]]] <- shifted[[cens[time]]]
  }

  out <- rbind(natural, shifted_half)
  out[["..i..lmtp_stack_indicator"]] <- rep(c(0, 1), each = nrow(natural))
  out
}
