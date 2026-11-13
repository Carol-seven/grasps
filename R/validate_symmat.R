#' Validate a symmetric matrix
#'
#' @description
#' Validate that an object is a finite numeric square symmetric matrix.
#'
#' @param x The object to validate.
#'
#' @param name A character string specifying the argument name to display
#' in the error message.
#'
#' @param min.dim An integer (default = 1) specifying the minimum matrix
#' dimension.
#'
#' @return
#' Invisibly returns \code{NULL} if \code{x} is valid; otherwise, an error is
#' raised.
#'
#' @noRd

validate_symmat <- function(x, name, min.dim = 1) {

  if (!is.matrix(x) || !is.numeric(x) || is.complex(x) ||
      nrow(x) < min.dim || nrow(x) != ncol(x) ||
      anyNA(x) || any(!is.finite(x))) {
    stop(sprintf(paste("`%s` must be a finite numeric square matrix",
                       "with dimension at least %d."), name, min.dim),
         call. = FALSE)
  }

  if (!isTRUE(all.equal(x, t(x)))) {
    stop(sprintf("`%s` must be symmetric.", name), call. = FALSE)
  }

  invisible(NULL)
}
