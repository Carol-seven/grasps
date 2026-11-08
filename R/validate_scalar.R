# Internal validation shared by the scalar controls in grasps().
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
    stop(sprintf("`%s` must be a finite %s %s %s%s.", name, type, bound,
                 lower, if (is.finite(upper)) paste0(" and at most ", upper) else ""),
         call. = FALSE)
  }
  invisible(NULL)
}
