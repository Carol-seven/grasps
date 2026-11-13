#' Draw Random Samples from a Specified Sampler (Random Variate Generator Wrapper)
#'
#' @description
#' Generate random samples from either a base R distribution or a user-supplied
#' sampling function.
#'
#' @param sampler A specification of the sampling distribution, either:
#' \enumerate{
#' \item A user-supplied sampling function. The function must accept an argument
#' \code{n} specifying the number of samples.
#' \item A character string specifying the distribution family.
#' Supported distribution families (corresponding base R generators
#' in parentheses) include:
#' \itemize{
#' \item "beta": Beta distribution (\code{\link[stats]{rbeta}}).
#' \item "cauchy": Cauchy distribution (\code{\link[stats]{rcauchy}}).
#' \item "chisq": Chi-squared distribution (\code{\link[stats]{rchisq}}).
#' \item "exp": Exponential distribution (\code{\link[stats]{rexp}}).
#' \item "f": F distribution (\code{\link[stats]{rf}}).
#' \item "gamma": Gamma distribution (\code{\link[stats]{rgamma}}).
#' \item "lnorm": Log normal distribution (\code{\link[stats]{rlnorm}}).
#' \item "norm": Normal distribution (\code{\link[stats]{rnorm}}).
#' \item "t": Student's t distribution (\code{\link[stats]{rt}}).
#' \item "unif": Uniform distribution (\code{\link[stats]{runif}}).
#' \item "weibull": Weibull distribution (\code{\link[stats]{rweibull}}).
#' }
#' }
#'
#' @param para A list or vector (default = \code{list()}) specifying additional
#' parameters passed to \code{sampler}.
#'
#' @param n An integer specifying the number of samples to generate.
#'
#' @return
#' A finite numeric vector of length \eqn{n}.
#'
#' @noRd

draw_sample <- function(sampler, para = list(), n) {

  ## construct argument list
  args <- c(list(n = n), para)

  if (is.function(sampler)) { ## user-defined sampler
    result <- tryCatch(
      do.call(sampler, args),
      error = function(e) {
        stop("Error in custom sampler: ", conditionMessage(e),
             call. = FALSE)
      }
    )
  } else if (is.character(sampler) && length(sampler) == 1L &&
             !is.na(sampler) && nzchar(sampler)) { ## base R distribution
    ## match distribution keyword to base R generator (e.g., "norm" -> "rnorm")
    rfun_name <- paste0("r", sampler)
    rfun <- get0(rfun_name, envir = asNamespace("stats"),
                 mode = "function", inherits = FALSE)
    if (is.null(rfun)) {
      stop(sprintf("Unknown sampler '%s' (expected base R generator '%s').",
                   sampler, rfun_name), call. = FALSE)
    }
    result <- tryCatch(
      do.call(rfun, args),
      error = function(e) {
        stop("Error calling ", rfun_name, ": ", conditionMessage(e),
             call. = FALSE)
      }
    )
  } else {
    stop("`sampler` must be a nonempty character scalar or a function.",
         call. = FALSE)
  }

  if (!is.numeric(result) || is.complex(result) || !is.null(dim(result)) ||
      length(result) != n || anyNA(result) || any(!is.finite(result))) {
    stop(sprintf("The sampler must return %d finite numeric value%s.", n,
                 if (n == 1L) "" else "s"), call. = FALSE)
  }

  return(result)
}
