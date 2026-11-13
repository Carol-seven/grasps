#' Validate a numeric vector argument
#'
#' @description
#' Validate a numeric vector argument against finiteness, range, and optional
#' missing-value requirements.
#' An error is raised if validation fails.
#'
#' @param value The value to validate.
#'
#' @param name A character string specifying the argument name to display
#' in the error message.
#'
#' @param lower A numeric scalar (default = -Inf) specifying the lower bound.
#'
#' @param upper A numeric scalar (default = Inf) specifying the inclusive upper
#' bound.
#'
#' @param lower.inclusive A logical value (default = TRUE) specifying whether
#' \code{lower} is included in the valid range.
#'
#' @param allow.na A logical value (default = FALSE) specifying whether missing
#' values are allowed.
#'
#' @return
#' Invisibly returns \code{NULL} if \code{value} is valid; otherwise, an error
#' is raised.
#'
#' @details
#' The input \code{value} must be a non-empty, non-complex numeric vector.
#' Non-missing values must be finite and lie in \code{(lower, upper]} when
#' \code{lower.inclusive = FALSE}, or in \code{[lower, upper]} when
#' \code{lower.inclusive = TRUE}.
#' Missing values are permitted only when \code{allow.na = TRUE}.
#' NaN values are not permitted.
#'
#' @noRd

validate_vector <- function(value, name, lower = -Inf, upper = Inf,
                            lower.inclusive = TRUE, allow.na = FALSE) {

  numeric <- is.numeric(value) && !is.complex(value)
  missing.only <- allow.na && is.logical(value) && all(is.na(value))
  valid <- (numeric || missing.only) &&
    is.null(dim(value)) && length(value) > 0L &&
    !any(is.nan(value)) && (allow.na || !anyNA(value))
  if (valid) {
    observed <- value[!is.na(value)]
    valid <- all(is.finite(observed)) && all(observed <= upper)
    valid <- valid && if (lower.inclusive) {
      all(observed >= lower)
    } else {
      all(observed > lower)
    }
  }
  if (!valid) {
    lower.bound <- if (is.finite(lower)) {
      paste(if (lower.inclusive) "at least" else "greater than", lower)
    } else {
      ""
    }
    upper.bound <- if (is.finite(upper)) paste("at most", upper) else ""
    range <- paste(c(lower.bound, upper.bound)[c(lower.bound, upper.bound) != ""],
                   collapse = " and ")
    range <- if (nzchar(range)) paste(" that are", range) else ""
    missing.text <- if (allow.na) " or `NA`" else ""
    stop(sprintf(paste("`%s` must be a non-empty numeric vector containing",
                       "finite values%s%s."), name, missing.text, range),
         call. = FALSE)
  }
  invisible(NULL)
}
