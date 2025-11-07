#' Draw Random Samples from a Specified Distribution (Random Variate Generator Wrapper)
#'
#' @description
#' Generate random samples from a base R distribution.
#'
#' @param dist A string specifying the distribution family. Supported
#' distributions include (corresponding base R functions in parentheses):
#' \enumerate{
#' \item "beta": Beta distribution (\code{\link[stats]{rbeta}})
#' \item "cauchy": Cauchy distribution (\code{\link[stats]{rcauchy}}).
#' \item "chisq": Chi-squared distribution (\code{\link[stats]{rchisq}}).
#' \item "exp": Exponential distribution (\code{\link[stats]{rexp}}).
#' \item "f": F distribution (\code{\link[stats]{rf}}).
#' \item "gamma": Gamma distribution (\code{\link[stats]{rgamma}}).
#' \item "lnorm": Log normal distribution (\code{\link[stats]{rlnorm}}).
#' \item "norm: Normal distribution (\code{\link[stats]{rnorm}}).
#' \item "t": Student's t distribution (\code{\link[stats]{rt}}).
#' \item "unif: Uniform distribution (\code{\link[stats]{runif}}).
#' \item "weibull": Weibull distribution (\code{\link[stats]{rweibull}}).
#' }
#'
#' @param para A list specifying the required parameters for the chosen
#' distribution.
#'
#' @param n An integer specifying the number of samples to generate.
#'
#' @return
#' A numeric vector of length \code{n} containing simulated data.
#'
#' @noRd

draw_dist <- function(dist, para, n) {
  ## match distribution keyword to base R generator (e.g., "norm" -> "rnorm")
  rfun <- match.fun(paste0("r", dist))
  ## construct argument list
  args <- c(list(n = n), para)
  ## call the random generator with supplied arguments
  do.call(rfun, args)
}

