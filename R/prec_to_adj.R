#' Adjacency Matrix from Precision Matrix
#'
#' @description
#' Convert a precision matrix to a partial-correlation-based adjacency matrix.
#'
#' @param prec.mat A numeric symmetric \eqn{p \times p} matrix, with \eqn{p} > 1,
#' specifying the precision matrix.
#'
#' @param diag.zero A logical value (default = TRUE) specifying whether to
#' set the diagonal entries of the adjacency matrix to 0.
#' If \code{diag.zero = FALSE}, the diagonal entries are set to 1 for a weighted
#' network. For an unweighted network (\code{weighted = FALSE}), the diagonal is
#' always forced to 0 to avoid self-loops.
#'
#' @param absolute A logical value (default = FALSE) specifying whether to
#' take the absolute values of the partial correlations.
#'
#' @param threshold A non-negative numeric value (default = \code{NULL})
#' specifying the threshold for edge filtering.
#' Entries with absolute values smaller than the threshold are set to 0.
#'
#' @param weighted A logical value (default = TRUE) specifying whether to
#' return a weighted adjacency matrix.
#' If \code{weighted = FALSE}, the matrix is a binary adjacency matrix with
#' entries equal to 0 or 1.
#'
#' @return
#' A numeric adjacency matrix with S3 class \code{"adjmat"}.
#'
#' @details
#' For a precision matrix \eqn{\Omega}, the partial correlation between nodes
#' \eqn{i} and \eqn{j} is computed as
#' \deqn{\rho_{ij} = - \frac{\Omega_{ij}}{\sqrt{\Omega_{ii}\Omega_{jj}}}.}
#'
#' The input matrix is not required to be positive definite.
#' If \code{prec.mat} is not positive definite, the function issues a warning
#' but continues the conversion.
#' In that case, the result may not represent valid partial correlations.
#'
#' @example
#' inst/example/ex-grasps.R
#'
#' @export

prec_to_adj <- function(prec.mat, diag.zero = TRUE, absolute = FALSE,
                        threshold = NULL, weighted = TRUE) {

  ## validate `prec.mat`
  validate_symmat(prec.mat, "prec.mat", min.dim = 2)
  d <- diag(prec.mat)
  if (any(d <= 0)) {
    stop("The diagonal entries of `prec.mat` must be positive.", call. = FALSE)
  }

  ## validate `diag.zero`, `absolute`, `weighted`
  controls <- list(diag.zero = diag.zero,
                   absolute = absolute,
                   weighted = weighted)
  for (name in names(controls)) {
    value <- controls[[name]]
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(sprintf("`%s` must be TRUE or FALSE.", name), call. = FALSE)
    }
  }

  ## validate `threshold`
  if (!is.null(threshold)) {
    validate_scalar(threshold, "threshold", lower = 0, lower.inclusive = TRUE)
  }

  ## check positive definiteness
  chol.prec <- tryCatch(chol(prec.mat), error = function(e) NULL)
  if (is.null(chol.prec)) {
    warning("`prec.mat` is not positive definite. ",
            "The result may not represent valid partial correlations.",
            call. = FALSE)
  }

  ## partial-correlation-based adjacency matrix
  norm_mat <- outer(sqrt(d), sqrt(d), "*")
  adj_mat <- -prec.mat / norm_mat

  ## optional: absolute value
  if (absolute) {
    adj_mat <- abs(adj_mat)
  }

  ## optional: thresholding
  if (!is.null(threshold)) {
    adj_mat[abs(adj_mat) < threshold] <- 0
  }

  ## optional: binary adjacency matrix
  ## optional: diagonal entries
  if (weighted) {
    diag(adj_mat) <- if (diag.zero) 0 else 1
  } else {
    adj_mat <- 1 * (adj_mat != 0)
    diag(adj_mat) <- 0
    if (!diag.zero) {
      message("For an unweighted network, the diagonal is forced to 0 ",
              "to avoid self-loops.")
    }
  }

  class(adj_mat) <- c("adjmat")
  return(adj_mat)
}
