#' Validate a scalar control argument
#'
#' @description
#' Validates a scalar numeric argument against finiteness, range, and optional
#' integer-value requirements.
#' An error is raised if validation fails.
#'
#' @param value The value to validate.
#'
#' @param name A character string specifying the argument name to display
#' in the error message.
#'
#' @param lower A numeric scalar (default = 0) specifying the lower bound.
#'
#' @param upper A numeric scalar (default = Inf) specifying the inclusive upper
#' bound.
#'
#' @param integer A logical value (default = FALSE) specifying whether
#' \code{value} must be integer-valued and no greater than
#' \code{.Machine$integer.max}.
#'
#' @param lower.inclusive A logical value (default = FALSE) specifying whether
#' \code{lower} is included in the valid range.
#'
#' @return
#' Invisibly returns \code{NULL} if \code{value} is valid; otherwise, an error
#' is raised.
#'
#' @details
#' The input \code{value} must be a finite, non-complex numeric scalar.
#' Its valid range is \code{(lower, upper]} when \code{lower.inclusive = FALSE},
#' and \code{[lower, upper]} when \code{lower.inclusive = TRUE}.
#' When \code{integer = TRUE}, \code{value} must additionally be integer-valued
#' and representable as an R integer.
#'
#' @noRd

validate_scalar <- function(value, name, lower = 0, upper = Inf,
                            integer = FALSE, lower.inclusive = FALSE) {
  valid <- is.numeric(value) && !is.complex(value) && length(value) == 1L &&
    is.finite(value)
  if (valid) {
    valid <- if (lower.inclusive) value >= lower else value > lower
    valid <- valid && value <= upper
    if (integer) {
      valid <- valid && value == floor(value) && value <= .Machine$integer.max
    }
  }
  if (!valid) {
    type <- if (integer) "integer" else "numeric scalar"
    bound <- if (lower.inclusive) "at least" else "greater than"
    stop(sprintf("`%s` must be a finite %s %s %s%s.", name, type, bound, lower,
                 if (is.finite(upper)) paste0(" and at most ", upper) else ""),
         call. = FALSE)
  }
  invisible(NULL)
}
