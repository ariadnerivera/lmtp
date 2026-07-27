theta_dr <- function(task, sequential_regressions, density_ratios, fits_m, fits_r, shift, is_sdr) {
  influence_function <- eif(density_ratios, sequential_regressions$shifted, sequential_regressions$natural)

  if (is_sdr) {
    theta <- weighted.mean(influence_function, task$weights)
  } else {
    theta <- weighted.mean(sequential_regressions$shifted[, 1], task$weights)
  }

  influence_function <- task$rescale(influence_function)
  theta <- task$rescale(theta)

  # With few independent clusters a normal critical value is optimistic.
  n_clusters <- length(unique(task$id))
  crit <- if (n_clusters > 1 && n_clusters < length(task$id)) {
    stats::qt(0.975, df = n_clusters - 1L)
  } else {
    stats::qnorm(0.975)
  }
  
  out <- list(
    estimator = ifelse(is_sdr, "SDR", "TMLE"),
    estimate = ife::ife(theta, influence_function, task$weights, as.character(task$id),
                        critical_value = crit),
    shift = shift,
    outcome_reg = task$rescale(sequential_regressions$shifted),
    density_ratios = density_ratios,
    fits_m = fits_m,
    fits_r = fits_r,
    outcome_type = task$outcome_type
  )

  class(out) <- "lmtp"
  out
}
