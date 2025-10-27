#' Sparse-Group Graphical Model
#'
#' @description
#' Provide a collection of statistical methods that incorporate both
#' element-wise and group-wise penalties to estimate a precision matrix.
#'
#' @param X \enumerate{
#' \item An n-by-p data matrix with sample size n and dimension p.
#' \item A p-by-p sample covariance/correlation matrix with dimension p.
#' }
#'
#' @param n An integer (default = the number of columns of \code{X}) specifying
#' the sample size. This is only required when the input matrix \code{X} is a
#' p-by-p sample covariance/correlation matrix with dimension p.
#'
#' @param groups An integer vector specifying the group membership. The group
#' sizes must be consistent with the dimension p.
#'
#' @param penalty A character string specifying the penalty for estimating
#' precision matrix. Available options include: \enumerate{
#' \item "adapt": adaptive lasso \insertCite{zou2006adaptive,fan2009network}{grasps}.
#' \item "lasso": lasso \insertCite{tibshirani1996regression,friedman2008sparse}{grasps}.
#' \item "mcp": minimax concave penalty \insertCite{zhang2010nearly}{grasps}.
#' \item "scad": smoothly clipped absolute deviation \insertCite{fan2001variable,fan2009network}{grasps}.
#' }
#'
#' @param diag.ind A boolean (default = TRUE) specifying whether to penalize
#' the diagonal elements.
#'
#' @param diag.grp A boolean (default = TRUE) specifying whether to penalize
#' the within-group blocks.
#'
#' @param diag.include A boolean (default = FALSE) specifying whether to include
#' the diagonal entries in the penalty for within-group blocks when
#' \code{diag.grp = TRUE}.
#'
#' @param lambda A grid of non-negative scalars for the regularization parameter.
#' The default is \code{NULL}, which generates its own \code{lambda} sequence
#' based on \code{nlambda} and \code{lambda.min.ratio}.
#'
#' @param alpha A grid of scalars in [0,1] specifying the parameter leveraging
#' the element-wise individual L1 penalty and the block-wise group L2 penalty.
#' An alpha of 1 corresponds to the element penalty only; an alpha of 0
#' corresponds to the group penalty only. The default values is a sequence from
#' 0.05 to 0.95 with increments of 0.05.
#'
#' @param gamma A scalar specifying the hyperparameter for the chosen
#' \code{penalty}. Default values: \enumerate{
#' \item "adapt": 0.5
#' \item "mcp": 3
#' \item "scad": 3.7
#' }
#'
#' @param nlambda An integer (default = 10) specifying the number of
#' \code{lambda} values to be generated when \code{lambda = NULL}.
#'
#' @param lambda.min.ratio A scalar (default = 0.01) specifying the fraction of
#' the maximum \code{lambda} value \eqn{\lambda_{max}} to generate the minimum
#' \code{lambda} \eqn{\lambda_{min}}. If \code{lambda = NULL}, the program
#' automatically generates a \code{lambda} grid as a sequence of length
#' \code{nlambda} in log scale, starting from \eqn{\lambda_{min}} to
#' \eqn{\lambda_{max}}.
#'
#' @param growiter.lambda An integer (default = 30) specifying the maximum
#' number of exponential growth steps during the initial search for an
#' admissible upper bound \eqn{\lambda_{\max}}.
#'
#' @param tol.lambda A scalar (default = 1e-03) specifying the relative
#' tolerance for the bisection stopping rule on the interval width.
#'
#' @param maxiter.lambda An integer (default = 50) specifying the maximum number
#' of bisection iterations in the line search for  \eqn{\lambda_{\max}}.
#'
#' @param rho A scalar > 0 (default = 2) specifying the ADMM
#' augmented-Lagrangian penalty parameter (often called the ADMM step size).
#' Larger values typically put more weight on enforcing the consensus
#' constraints each iteration; smaller values yield more conservative updates.
#'
#' @param tau.incr A scalar > 1 (default = 2) specifying the multiplicative
#' factor used to increase \code{rho} when the primal residual dominates the
#' dual residual in ADMM.
#'
#' @param tau.decr A scalar > 1 (default = 2) specifying the multiplicative
#' factor used to decrease \code{rho} when the dual residual dominates the
#' primal residual in ADMM.
#'
#' @param nu A scalar > 1 (default = 10) controlling how aggressively \code{rho}
#' is rescaled in the adaptive-\code{rho} scheme (residual balancing).
#'
#' @param tol.abs A scalar > 0 (default = 1e-04) specifying the absolute
#' tolerance for ADMM stopping (applied to primal/dual residual norms).
#'
#' @param tol.rel A scalar > 0 (default = 1e-04) specifying the relative
#' tolerance for ADMM stopping (applied to primal/dual residual norms).
#'
#' @param maxiter An integer (default = 1e+04) specifying the maximum number of
#' ADMM iterations.
#'
#' @param crit A string (default = "BIC") specifying the parameter selection
#' method to use. Available options include: \enumerate{
#' \item "AIC": Akaike information criterion \insertCite{akaike1973information}{grasps}.
#' \item "BIC": Bayesian information criterion \insertCite{schwarz1978estimating}{grasps}.
#' \item "EBIC": extended Bayesian information criterion \insertCite{foygel2010extended}{grasps}.
#' \item "HBIC": high dimensional Bayesian information criterion \insertCite{wang2013calibrating,fan2017high}{grasps}.
#' \item "CV": k-fold cross validation with negative log-likelihood loss.
#' }
#'
#' @param kfold An integer (default = 5) specifying the number of folds used for
#' \code{crit = "CV"}.
#'
#' @param ebic.tuning A scalar (default = 0.5) specifying the tuning parameter to
#' calculate for \code{crit = "EBIC"}.
#'
#' @importFrom stats cov
#' @importFrom Rdpack reprompt
#'
#' @return
#' A list containing the following components:
#' \describe{
#' \item{hatOmega}{The estimated precision matrix.}
#' \item{lambda}{The optimal regularization parameter.}
#' \item{alpha}{The optimal penalty balancing parameter.}
#' \item{initial}{The initial estimate of \code{hatOmega} when \code{penalty} is
#' set to \code{"adapt"}, \code{"mcp"}, or \code{"scad"}.}
#' \item{gamma}{The optimal hyperparameter when \code{penalty} is set to
#' \code{"adapt"}, \code{"mcp"}, or \code{"scad"}.}
#' \item{iterations}{The number of ADMM iterations.}
#' \item{lambda.grid}{The actual lambda grid used in the program.}
#' \item{alpha.grid}{The actual alpha grid used in the program.}
#' \item{lambda.safe}{The bisection-refined upper bound \eqn{\lambda_{\max}},
#' corresponding to \code{alpha.grid}, when \code{lambda = NULL}.}
#' \item{loss}{The optimal k-fold loss when \code{crit = "CV"}.}
#' \item{CV.loss}{Matrix of CV losses, with rows for parameter combinations and
#' columns for CV folds, when \code{crit = "CV"}.}
#' \item{score}{The optimal information criterion score when \code{crit} is set
#' to \code{"AIC"}, \code{"BIC"}, \code{"EBIC"}, or \code{"HBIC"}.}
#' \item{IC.score}{The information criterion score for each parameter
#' combination when \code{crit} is set to \code{"AIC"}, \code{"BIC"},
#' \code{"EBIC"}, or \code{"HBIC"}.}
#' }
#'
#' @references
#' \insertAllCited{}
#'
#' @export

sggm <- function(X, n = nrow(X), groups, penalty,
                 diag.ind = TRUE, diag.grp = TRUE, diag.include = FALSE,
                 lambda = NULL, alpha = NULL, gamma = NULL,
                 nlambda = 10, lambda.min.ratio = 0.01,
                 growiter.lambda = 30, tol.lambda = 1e-03, maxiter.lambda = 50,
                 rho = 2, tau.incr = 2, tau.decr = 2, nu = 10,
                 tol.abs = 1e-04, tol.rel = 1e-04, maxiter = 1e+04,
                 crit = "BIC", kfold = 5, ebic.tuning = 0.5) {

  p <- ncol(X)

  if (length(groups) != p) {
    stop("The length of 'groups' must be the column dimension of X!\n");
  }
  if (!penalty %in% c("lasso", "adapt", "mcp", "scad")) {
    stop("Error in penalty! Available options: 'lasso', 'adapt', 'mcp', 'scad'.\n")
  }
  if (!all(lambda > 0)) {
    stop("The parameter 'lambda' must be positive!\n")
  }
  if (!all(alpha >= 0 & alpha <= 1)) {
    stop("The parameter 'alpha' must be in [0,1]!\n")
  }
  if (rho <= 0) {
    stop("The parameter 'rho' must be positive!\n")
  }
  if (tau.incr <= 1) {
    stop("The parameter 'tau.incr' must be greater than 1!\n");
  }
  if (tau.decr <= 1) {
    stop("The parameter 'tau.decr' must be greater than 1!\n");
  }
  if (nu <= 1) {
    stop("The parameter 'nu' must be greater than 1!\n");
  }
  if (tol.abs <= 0) {
    stop("The parameter 'tol.abs' must be positive!\n");
  }
  if (tol.rel <= 0) {
    stop("The parameter 'tol.rel' must be positive!\n");
  }

  if (isSymmetric(X)) {
    S <- X
    X <- NULL
  } else {
    S <- (n-1)/n*cov(X)
  }

  grp <- sort(unique(groups))
  group.idx <- lapply(grp, function(g) which(groups == g)-1)
  names(group.idx) <- grp

  alpha_null <- is.null(alpha)
  lambda_null <- is.null(lambda)

  if (alpha_null) {
    alpha <- seq(0.05, 0.95, 0.05)
  }

  if (lambda_null) {
    max_abs.off <- max(abs(S[upper.tri(S, diag = FALSE)]))
    S2 <- S*S
    sum.row <- rowsum(S2, groups)
    block.sum <- rowsum(t(sum.row), groups)
    block.norm <- sqrt(block.sum)
    max_block.norm.off <- max(block.norm[upper.tri(block.norm, diag = FALSE)])
    idx_list <- split(seq_along(groups), groups)
    block.off_list <- list()
    for (g in seq_along(idx_list)) {
      for (gp in seq_along(idx_list)) {
        if (g < gp) {
          block.off_list[[length(block.off_list)+1]] <- S[idx_list[[g]], idx_list[[gp]], drop = FALSE]
        }
      }
    }
    lambda_path <- lapply(alpha, function(a) {
      lambda.ind <- ifelse(a > 0, max_abs.off / a, 0)
      lambda.grp <- ifelse(a < 1, max_block.norm.off / (1-a), 0)
      lambda.safe <- max(lambda.ind, lambda.grp)
      if (a > 0 & a < 1) {
        lambda.max <- line_search_lambda_max(block.off_list, lambda.safe, a, growiter.lambda, tol.lambda, maxiter.lambda)
      } else {
        lambda.max <- lambda.safe
      }
      lambda.min <- lambda.min.ratio*lambda.max
      return(list(lambda = exp(seq(log(lambda.max), log(lambda.min), length = nlambda)),
                  lambda.safe = lambda.safe))
    })
    lambda <- unlist(lapply(lambda_path, `[[`, "lambda"))
    lambda.safe <- sapply(lambda_path, `[[`, "lambda.safe")
    names(lambda.safe) <- alpha
  }

  if (lambda_null & alpha_null) {
    parameter <- data.frame(alpha = rep(alpha, each = nlambda),
                            lambda = lambda)
  } else {
    parameter <- data.frame(alpha = alpha, lambda = lambda)
  }

  if (is.null(gamma)) {
    if (penalty == "adapt") {
      gamma <- 0.5
    } else if (penalty == "mcp") {
      gamma <- 3
    } else if (penalty == "scad") {
      gamma <- 3.7
    } else {
      gamma <- NA
    }
  }

  if (nrow(parameter) > 1) {

    if (crit == "CV") {

      if(is.null(X)) {
        stop("CV requires the n-by-p data matrix!")
      }

      if (kfold < 2 | kfold > n) {
        stop("'kfold' must be between 2 and the row dimension of X!")
      }

      CV <- ADMMsggm_CV(X = X, group_idx = group.idx, penalty = penalty,
                        diag_ind = diag.ind, diag_grp = diag.grp, diag_include = diag.include,
                        lambdas = parameter$lambda, alphas = parameter$alpha, gamma = gamma,
                        rho = rho, tau_incr = tau.incr, tau_decr = tau.decr, nu = nu,
                        tol_abs = tol.abs, tol_rel = tol.rel, maxiter = maxiter,
                        kfold = kfold)

      result <- ADMMsggm(S = S, group_idx = group.idx, penalty = penalty,
                         diag_ind = diag.ind, diag_grp = diag.grp, diag_include = diag.include,
                         lambda = CV$lambda.opt, alpha = CV$alpha.opt, gamma = gamma,
                         rho = rho, tau_incr = tau.incr, tau_decr = tau.decr, nu = nu,
                         tol_abs = tol.abs, tol_rel = tol.rel, maxiter = maxiter)
      result$lambda.grid <- parameter$lambda
      result$alpha.grid <- parameter$alpha
      if (lambda_null) {
        result$lambda.safe <- lambda.safe
      }
      result$loss <- CV$loss.min
      result$CV.loss <- CV$CV.loss

    } else {

      IC <- ADMMsggm_IC(S = S, group_idx = group.idx, penalty = penalty,
                        diag_ind = diag.ind, diag_grp = diag.grp, diag_include = diag.include,
                        lambdas = parameter$lambda, alphas = parameter$alpha, gamma = gamma,
                        rho = rho, tau_incr = tau.incr, tau_decr = tau.decr, nu = nu,
                        tol_abs = tol.abs, tol_rel = tol.rel, maxiter = maxiter,
                        crit = crit, n = n, ebic_tuning = ebic.tuning)

      result <- IC$result
      result$lambda.grid <- parameter$lambda
      result$alpha.grid <- parameter$alpha
      if (lambda_null) {
        result$lambda.safe <- lambda.safe
      }
      result$score <- IC$score.min
      result$IC.score <- as.vector(IC$IC.score)

    }

  } else {

    result <- ADMMsggm(S = S, group_idx = group.idx, penalty = penalty,
                       diag_ind = diag.ind, diag_grp = diag.grp, diag_include = diag.include,
                       lambda = lambda, alpha = alpha, gamma = gamma,
                       rho = rho, tau_incr = tau.incr, tau_decr = tau.decr, nu = nu,
                       tol_abs = tol.abs, tol_rel = tol.rel, maxiter = maxiter)

    if (lambda_null) {
      result$lambda.safe <- lambda.safe
    }

  }

  return(result)

}

