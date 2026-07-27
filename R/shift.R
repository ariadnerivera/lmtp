make_shifted <- function(data, trt, cens, shift, shifted, cluster = NULL,
                         shift_cluster = NULL,
                         patient_shift_uses_shifted_cluster = FALSE) {
  assert_function(shift, nargs = 2, null.ok = TRUE)
  assert_function(shift_cluster, nargs = 2, null.ok = TRUE)
  
  if (!is.null(shifted)) {
    assert_correctly_shifted(data, shifted, c(unlist(trt), unlist(cluster)), cens)
    return(shifted)
  }
  
  if (is.null(cluster)) {
    return(shift_data(data, trt, shift))
  }
  
  if (is.null(shift_cluster)) {
    stop("`shift_cluster` must be supplied when `cluster` is supplied, unless a ",
         "precomputed `shifted` data frame is used.", call. = FALSE)
  }
  
  if (isTRUE(patient_shift_uses_shifted_cluster)) {
    # V is intervened on first; the A policy sees V^d.
    out <- shift_data(data, cluster, shift_cluster)
    out <- shift_data(out, trt, shift)
  } else {
    # A policy uses the natural V (lmtp's usual convention).
    out <- shift_data(data, trt, shift)
    shifted_V <- shift_data(data, cluster, shift_cluster)
    out[, unlist(cluster)] <- shifted_V[, unlist(cluster), drop = FALSE]
  }
  out
}

shift_data <- function(data, trt, shift) {
  if (is.null(shift)) {
    return(data)
  }

  is_multivariate <- is.list(trt)
  if (isTRUE(is_multivariate)) {
    return(shift_trt_list(data, trt, shift))
  }

  shift_trt_character(data, trt, shift)
}


shift_trt_character <- function(data, trt, .f) {
  out <- as.list(data)
  for (a in trt) {
    out[[a]] <- .f(data, a)
  }
  as.data.frame(out, check.names = FALSE)
}

shift_trt_list <- function(data, trt, .f) {
  out <- as.list(data)
  for (a in trt) {
    new <- .f(data, a)
    if (!is.list(new) && !is.data.frame(new)) {
      if (length(a) != 1L) {
        stop("A shift function for a multivariate exposure node must return a named ",
             "list or data.frame with one element per column.", call. = FALSE)
      }
      new <- stats::setNames(list(new), a)
    }
    for (col in a) {
      out[[col]] <- new[[col]]
    }
  }
  as.data.frame(out, check.names = FALSE)
}

#' Turn All Treatment Nodes On
#'
#' A pre-packaged shift function for use with provided estimators when the exposure is binary.
#' Used to estimate the population intervention effect when all treatment variables are set to 1.
#'
#' @param data A dataframe containing the treatment variables.
#' @param trt The name of the current treatment variable.
#'
#' @seealso [lmtp_tmle()], [lmtp_sdr()]
#' @return A dataframe with all treatment nodes set to 1.
#' @export
#'
#' @examples
#' \donttest{
#' data("iptwExWide", package = "twang")
#' a <- paste0("tx", 1:3)
#' baseline <- c("gender", "age")
#' tv <- list(c("use0"), c("use1"), c("use2"))
#' lmtp_sdr(iptwExWide, a, "outcome", baseline = baseline, time_vary = tv,
#'          shift = static_binary_on, outcome_type = "continuous", folds = 2)
#' }
static_binary_on <- function(data, trt) {
  rep(1, length(data[[trt]]))
}

#' Turn All Treatment Nodes Off
#'
#' A pre-packaged shift function for use with provided estimators when the exposure is binary.
#' Used to estimate the population intervention effect when all treatment variables are set to 0.
#'
#' @param data A dataframe containing the treatment variables.
#' @param trt The name of the current treatment variable.

#' @seealso [lmtp_tmle()], [lmtp_sdr()]
#' @return A dataframe with all treatment nodes set to 0.
#' @export
#'
#' @examples
#' \donttest{
#' data("iptwExWide", package = "twang")
#' a <- paste0("tx", 1:3)
#' baseline <- c("gender", "age")
#' tv <- list(c("use0"), c("use1"), c("use2"))
#' lmtp_sdr(iptwExWide, a, "outcome", baseline = baseline, time_vary = tv,
#'          shift = static_binary_off, outcome_type = "continuous", folds = 2)
#' }
static_binary_off <- function(data, trt) {
  rep(0, length(data[[trt]]))
}

#' IPSI Function Factory
#'
#' A function factory that returns a shift function for increasing or decreasing
#' the probability of exposure when exposure is binary.
#'
#' @param delta \[\code{numeric(1)}\]\cr
#'  A risk ratio between 0 and Inf.
#'
#' @seealso [lmtp_tmle()], [lmtp_sdr()]
#' @return A shift function.
#' @export
#'
#' @examples
#' \donttest{
#' data("iptwExWide", package = "twang")
#' a <- paste0("tx", 1:3)
#' baseline <- c("gender", "age")
#' tv <- list(c("use0"), c("use1"), c("use2"))
#' lmtp_sdr(iptwExWide, a, "outcome", baseline = baseline, time_vary = tv,
#'          shift = ipsi(0.5), outcome_type = "continuous", folds = 2)
#' }
ipsi <- function(delta) {
  if (delta > 1) {
    return(ipsi_up(1 / delta))
  }
  ipsi_down(delta)
}

ipsi_up <- function(delta) {
  function(data, trt) {
    eps <- runif(nrow(data), 0, 1)
    ifelse(eps < delta, data[[trt]], 1)
  }
}

ipsi_down <- function(delta) {
  function(data, trt) {
    eps <- runif(nrow(data), 0, 1)
    ifelse(eps < delta, data[[trt]], 0)
  }
}
