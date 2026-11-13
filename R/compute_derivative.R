#' Penalty Derivative Computation
#'
#' @description
#' Compute one or more derivative values for a given \code{omega}, allowing
#' vectorized specifications of \code{penalty}, \code{lambda}, and \code{gamma}.
#'
#' @param omega A numeric vector specifying the values at which the penalty
#' derivative is evaluated.
#'
#' @param penalty A character vector specifying one or more penalty types.
#' Available options include:
#' \enumerate{
#' \item "lasso": Least absolute shrinkage and selection operator
#' \insertCite{tibshirani1996regression,friedman2008sparse}{grasps}.
#' \item "atan": Arctangent type penalty
#' \insertCite{wang2016variable}{grasps}.
#' \item "exp": Exponential type penalty
#' \insertCite{wang2018variable}{grasps}.
#' \item "lq": Lq penalty
#' \insertCite{frank1993statistical,fu1998penalized,fan2001variable}{grasps}.
#' \item "lsp": Log-sum penalty
#' \insertCite{candes2008enhancing}{grasps}.
#' \item "mcp": Minimax concave penalty
#' \insertCite{zhang2010nearly}{grasps}.
#' \item "scad": Smoothly clipped absolute deviation
#' \insertCite{fan2001variable,fan2009network}{grasps}.
#' }
#' If \code{penalty} has length 1, it is recycled to the common length
#' determined by \code{penalty}, \code{lambda}, and \code{gamma}.
#'
#' @param lambda A non-negative numeric vector specifying the regularization
#' parameter.
#' If \code{lambda} has length 1, it is recycled to the common length
#' determined by \code{penalty}, \code{lambda}, and \code{gamma}.
#'
#' @param gamma A numeric vector specifying the additional parameter for
#' the penalty function.
#' If \code{gamma} has length 1, it is recycled to the common length
#' determined by \code{penalty}, \code{lambda}, and \code{gamma}.
#' The default value and valid range depend on the penalty:
#' \enumerate{
#' \item "atan": default 0.005; must be positive.
#' \item "exp": default 0.01; must be positive.
#' \item "lq": default 0.5; must be in \eqn{(0, 1)}.
#' \item "lsp": default 0.1; must be positive.
#' \item "mcp": default 3; must be greater than 1.
#' \item "scad": default 3.7; must be greater than 2.
#' }
#' For \code{"lasso"}, \code{gamma} is ignored.
#'
#' @return
#' A data frame with S3 class \code{"derivative"} containing:
#' \describe{
#' \item{omega}{The input \code{omega} values.}
#' \item{penalty}{The penalty type for each row.}
#' \item{lambda}{The regularization parameter used.}
#' \item{gamma}{The additional penalty parameter used.}
#' \item{value}{The computed derivative value.}
#' }
#' The number of rows equals \code{length(omega)} multiplied by
#' \code{max(length(penalty), length(lambda), length(gamma))}.
#' Each of \code{penalty}, \code{lambda}, and \code{gamma} must have length 1
#' or the same common length.
#' Arguments of length 1 are recycled to this common length.
#'
#' @references
#' \insertAllCited{}
#'
#' @example
#' inst/example/ex-compute_derivative.R
#'
#' @export

compute_derivative <- function(omega, penalty, lambda, gamma = NA) {

  validate_vector(omega, "omega")
  penalties <- c("lasso", "atan", "exp", "lq", "lsp", "mcp", "scad")
  if (!is.character(penalty) || length(penalty) == 0L ||
      anyNA(penalty) || any(!penalty %in% penalties)) {
    stop(sprintf("Invalid `penalty`.\nAvailable options: %s.",
                 paste(sprintf('"%s"', penalties), collapse = ", ")),
         call. = FALSE)
  }
  validate_vector(lambda, "lambda", lower = 0)
  validate_vector(gamma, "gamma", allow.na = TRUE)

  para_len <- lengths(list(penalty, lambda, gamma))
  n <- max(para_len)
  if (any(!para_len %in% c(1L, n))) {
    stop("`penalty`, `lambda`, and `gamma` must each have length 1 ",
         "or the same common length.", call. = FALSE)
  }
  penalty <- rep(penalty, length.out = n)
  lambda <- rep(lambda, length.out = n)
  gamma <- rep(gamma, length.out = n)

  res <- do.call(rbind, lapply(seq_len(n), function(k) {
    deriv_internal(omega = omega,
                   penalty = penalty[k],
                   lambda = lambda[k],
                   gamma = gamma[k])
  }))
  class(res) <- c("derivative", "penderiv", class(res))

  return(res)
}


#' @noRd

deriv_internal <- function(omega, penalty, lambda, gamma) {

  ## default gamma by penalty
  if (penalty != "lasso") {
    if (missing(gamma) || is.na(gamma)) {
      gamma <- switch(penalty,
                      "atan" = 0.005, "exp" = 0.01, "lq" = 0.5,
                      "lsp" = 0.1, "mcp" = 3, "scad" = 3.7)
    }
    lower <- switch(penalty, "mcp" = 1, "scad" = 2, 0)
    validate_scalar(gamma, "gamma", lower = lower)
    if (penalty == "lq" && gamma >= 1) {
      stop('For "lq", `gamma` must be in (0,1).', call. = FALSE)
    }
  } else {
    gamma <- NA_real_
  }

  a <- abs(omega)

  if (penalty == "atan") {
    res <- lambda * gamma * (gamma + 2 / pi) / (gamma^2 + a^2)

  } else if (penalty == "exp") {
    res <- (lambda / gamma) * exp(-a / gamma)

  } else if (penalty == "lasso") {
    res <- lambda

  } else if (penalty == "lq") {
    epsilon <- 1e-10
    res <- lambda * gamma * (a + epsilon)^(gamma - 1)

  } else if (penalty == "lsp") {
    res <- lambda / (gamma + a)

  } else if (penalty == "mcp") {
    res <- (lambda - a / gamma) * (a <= gamma * lambda)

  } else if (penalty == "scad") {
    res <- lambda * (a <= lambda) +
      pmax(gamma * lambda - a, 0) / (gamma - 1) * (a > lambda)
  }

  return(data.frame(omega = omega,
                    penalty = penalty,
                    lambda = lambda,
                    gamma = gamma,
                    value = res))
}
