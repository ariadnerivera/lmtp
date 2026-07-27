get_folded_data <- function(data, folds, index) {
  out <- list()
  out[["train"]] <- data[folds[[index]]$training_set, , drop = FALSE]
  out[["valid"]] <- data[folds[[index]]$validation_set, , drop = FALSE]
  out
}

fix_censoring_ind <- function(data, cens) {
  if (is.null(cens)) {
    return(data)
  }

  data <- data.table::copy(data)
  for (cen in cens) {
    data.table::set(data, j = cen, value = ifelse(is.na(data[[cen]]), 0, data[[cen]]))
  }
  data
}

bound <- function(x, p = 1e-05) {
  pmax(pmin(x, 1 - p), p)
}

followed_rule <- function(natural, shifted, A, mtp) {
  if (mtp) {
    return(rep(TRUE, nrow(natural)))
  }

  followed <- matrix(nrow = nrow(natural), ncol = length(A))

  for (i in seq_along(A)) {
    a <- A[i]
    followed[, i] <- mapply(function(x, y) isTRUE(all.equal(x, y)),
                            as.list(natural[, a]),
                            as.list(shifted[, a]))
  }

  apply(followed, 1, prod)
}

trim <- function(x, trim) {
  pmin(x, quantile(x, trim, na.rm = TRUE))
}

is.lmtp <- function(x) {
  class(x) == "lmtp"
}

sw <- function(x) {
  suppressWarnings(x)
}

last <- function(x) {
  x[length(x)]
}

extract_sl_weights <- function(fit) {
  if (inherits(fit, "mlr3superlearner")) {
    return(cbind(Risk = fit$risk))
  }
  fit$coef
}

convert_to_surv <- function(x) {
  data.table::fcase(
    x == 0, 1,
    x == 1, 0
  )
}

is_normalized <- function(x, tolerance = .Machine$double.eps^0.5) {
  # Check if the mean is approximately 1 within the given tolerance
  abs(mean(x) - 1) < tolerance
}

fix_surv_time1 <- function(x) {
  to_fix <- x[[1]]$estimate
  x[[1]]$estimate <- ife::ife(1 - to_fix@x, 1 - to_fix@eif, to_fix@weights, to_fix@id)
  x
}

is_decimal <- function(x) {
  test <- floor(x)
  !(x == test)
}

`%and%` <- function(o, r) {
  i <- vector("logical", length(o))
  for (j in 1:length(o)) {
    if (is.na(r[j]) & !is.na(o[j])) {
      i[j] <- o[j]
    } else if (!is.na(r[j]) & is.na(o[j])) {
      i[j] <- r[j]
    } else {
      i[j] <- o[j] & r[j]
    }
  }
  i
}

current_trt <- function(trt, time) {
  if (length(trt) > 1) {
    return(trt[[time]])
  }
  trt[[1]]
}

  
# Outcome-regression conditioning set including the cluster exposure (V) active
# through the time point whose parents are being modelled.
outcome_history <- function(task, time) {
  unique(c(task$vars$history("L", time + 1), task$vars$cluster_cols_through(time)))
}
  
# Columns set to their shifted values when predicting the outcome under the
# intervention: the unit exposure A_t and (if present) the cluster exposure V_t.
shift_targets <- function(task, time) {
  A_t <- current_trt(task$vars$A, time)
  if (is.null(task$vars$V)) return(A_t)
  c(A_t, current_trt(task$vars$V, time))
}


# Coerce a cross-sectional cluster exposure to lmtp's one-time-point node format.
# `cluster = c("V1","V2")` and `cluster = list(c("V1","V2"))` both become
# list(c("V1","V2")). Applied to `cluster` ONLY - never to `trt`, because a bare
# character vector in `trt` legitimately means one column per time point.
normalize_cluster_exposure <- function(x) {
  if (is.null(x)) return(NULL)
  if (is.character(x)) return(list(x))
  if (is.list(x) && length(x) == 1L && is.character(x[[1]])) return(x)
  stop("`cluster` must be a character vector of column names, or a one-element ",
       "list containing one.", call. = FALSE)
}

# Referenced by Task.R's make_folds(); was missing from the package.
final_outcome <- function(Y) Y[length(Y)]
