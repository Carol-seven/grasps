#' Groupwise Regularized Adaptive Sparse Precision Solution
#'
#' @description
#' Provide a collection of statistical methods that incorporate both
#' element-wise and group-wise penalties to estimate a precision matrix.
#'
#' @param X \enumerate{
#' \item An \eqn{n \times p} data matrix with sample size \eqn{n} and
#' dimension \eqn{p}.
#' \item A \eqn{p \times p} sample covariance matrix with dimension \eqn{p}.
#' }
#' Entries must be finite and real.
#' Covariance input must be positive semidefinite with positive diagonal entries;
#' data columns must have positive sample variances.
#' Symmetric matrices are interpreted as covariance input.
#'
#' @param n An integer (default = \code{nrow(X)}) specifying the sample size.
#' \enumerate{
#' \item For an \eqn{n \times p} data matrix, it defaults to \code{nrow(X)} and,
#' if supplied, must equal \code{nrow(X)}.
#' \item For a \eqn{p \times p} sample covariance matrix, it must be supplied
#' explicitly.
#' }
#' In either case, \code{n} must be at least 2.
#'
#' @param membership A numeric vector specifying the group membership.
#' Its length must equal the dimension \eqn{p}.
#' Labels must be finite, non-missing integer values; consecutive labels are not
#' required.
#'
#' @param penalty A character string specifying the penalty used to estimate
#' the precision matrix.
#' Available options include:
#' \enumerate{
#' \item "lasso": Least absolute shrinkage and selection operator
#' \insertCite{tibshirani1996regression,friedman2008sparse}{grasps}.
#' \item "adapt": Adaptive lasso
#' \insertCite{zou2006adaptive,fan2009network}{grasps}.
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
#'
#' @param diag.ind A logical value (default = TRUE) specifying whether to
#' penalize the diagonal elements.
#'
#' @param diag.grp A logical value (default = TRUE) specifying whether to
#' penalize the within-group blocks.
#'
#' @param diag.include A logical value (default = FALSE) specifying whether to
#' include the diagonal entries in the penalty for within-group blocks when
#' \code{diag.grp = TRUE}.
#'
#' @param lambda A positive numeric vector specifying the grid for
#' the regularization parameter.
#' The default is \code{NULL}, which generates its own \code{lambda} sequence
#' based on \code{nlambda} and \code{lambda.min.ratio}.
#'
#' @param alpha A numeric vector in [0,1] specifying the grid for the mixing
#' parameter balancing the element-wise individual L1 penalty and the block-wise
#' group L2 penalty.
#' An alpha of 1 corresponds to the individual penalty only; an alpha of 0
#' corresponds to the group penalty only.
#' The default value is a sequence from 0.1 to 0.9 with increments of 0.1.
#' When both \code{alpha} and \code{lambda} are supplied,
#' all Cartesian combinations of their values are evaluated.
#'
#' @param gamma A finite numeric scalar specifying the additional parameter for
#' the chosen \code{penalty}.
#' The default value and valid range depend on the penalty:
#' \enumerate{
#' \item "adapt": default 0.5; must be positive.
#' \item "atan": default 0.005; must be positive.
#' \item "exp": default 0.01; must be positive.
#' \item "lq": default 0.5; must be in \eqn{(0, 1)}.
#' \item "lsp": default 0.1; must be positive.
#' \item "mcp": default 3; must be greater than 1.
#' \item "scad": default 3.7; must be greater than 2.
#' }
#' This argument is ignored when \code{penalty = "lasso"}.
#'
#' @param initial A character string (default = "lasso") or a finite symmetric
#' numeric \eqn{p \times p} matrix specifying the initial estimator used to
#' construct penalty weights.
#' Available character options include:
#' \enumerate{
#' \item "lasso": Use the sparse-group graphical lasso estimate.
#' \item "invS": Use the inverse of the empirical covariance matrix.
#' The relevant empirical covariance matrix must be invertible.
#' }
#' For \code{crit = "CV"}, character-specified initial estimators are recomputed
#' separately from each training fold.
#' A user-supplied matrix cannot be used when \code{crit = "CV"}.
#' This argument is ignored when \code{penalty = "lasso"}.
#'
#' @param nlambda An integer (default = 10) specifying the number of
#' \code{lambda} values to generate when \code{lambda = NULL}.
#'
#' @param lambda.min.ratio A numeric value > 0 (default = 0.01) specifying
#' the fraction of the maximum \code{lambda} value \eqn{\lambda_{\max}} to
#' generate the minimum \code{lambda} \eqn{\lambda_{\min}}.
#' If \code{lambda = NULL}, a \code{lambda} grid of length \code{nlambda} is
#' automatically generated on a log scale, ranging from \eqn{\lambda_{\max}}
#' down to \eqn{\lambda_{\min}}.
#'
#' @param growiter.lambda An integer (default = 30) specifying the maximum
#' number of exponential growth steps during the initial search for an
#' admissible upper bound \eqn{\lambda_{\max}}.
#'
#' @param tol.lambda A numeric value > 0 (default = 1e-03) specifying
#' the relative tolerance for the bisection stopping rule on the interval width.
#'
#' @param maxiter.lambda An integer (default = 50) specifying the maximum number
#' of bisection iterations in the line search for  \eqn{\lambda_{\max}}.
#'
#' @param lla A logical value (default = FALSE) specifying whether to use
#' iterative local linear approximation (LLA) \insertCite{zou2008one}{grasps}
#' for non-convex penalties.
#' \enumerate{
#' \item If \code{FALSE}, a one-step LLA estimator is computed.
#' \item If \code{TRUE}, penalty weights are updated iteratively until the LLA
#' stopping criterion is satisfied or \code{maxiter.lla} is reached.
#' }
#' This argument is ignored for \code{penalty = "lasso"} and
#' \code{penalty = "adapt"}.
#'
#' @param tol.lla A numeric value > 0 (default = 1e-04) specifying the stopping
#' tolerance for the relative change between successive precision matrix
#' estimates in iterative LLA.
#' This argument is used only when \code{lla = TRUE}.
#'
#' @param maxiter.lla An integer (default = 10) specifying the maximum number of
#' weighted subproblems solved by iterative LLA.
#' Construction of the initial estimator is not included in this count.
#' This argument is used only when \code{lla = TRUE}.
#'
#' @param rho A numeric value > 0 (default = 2) specifying the ADMM
#' augmented-Lagrangian penalty parameter (often called the ADMM step size).
#' Larger values typically put more weight on enforcing the consensus
#' constraints at each iteration; smaller values yield more conservative updates.
#'
#' @param tau.incr A numeric value > 1 (default = 2) specifying
#' the multiplicative factor used to increase \code{rho} when the primal
#' residual dominates the dual residual in ADMM.
#'
#' @param tau.decr A numeric value > 1 (default = 2) specifying
#' the multiplicative factor used to decrease \code{rho} when the dual residual
#' dominates the primal residual in ADMM.
#'
#' @param nu A numeric value > 1 (default = 10) controlling how aggressively
#' \code{rho} is rescaled in the adaptive-\code{rho} scheme (residual balancing).
#'
#' @param tol.abs A numeric value > 0 (default = 1e-04) specifying the absolute
#' tolerance for ADMM stopping (applied to primal/dual residual norms).
#'
#' @param tol.rel A numeric value > 0 (default = 1e-04) specifying the relative
#' tolerance for ADMM stopping (applied to primal/dual residual norms).
#'
#' @param maxiter An integer (default = 1e+04) specifying the maximum number of
#' ADMM iterations.
#'
#' @param crit A character string (default = "BIC") specifying the parameter
#' selection criterion to use. Available options include:
#' \enumerate{
#' \item "AIC": Akaike information criterion
#' \insertCite{akaike1973information}{grasps}.
#' \item "BIC": Bayesian information criterion
#' \insertCite{schwarz1978estimating}{grasps}.
#' \item "EBIC": extended Bayesian information criterion
#' \insertCite{chen2008extended,foygel2010extended}{grasps}.
#' \item "HBIC": high dimensional Bayesian information criterion
#' \insertCite{wang2013calibrating,fan2017high}{grasps}.
#' \item "CV": k-fold cross validation with negative log-likelihood loss.
#' }
#' This argument is used only when the parameter grid contains more than one
#' candidate.
#'
#' @param kfold An integer (default = 5) specifying the number of folds used for
#' \code{crit = "CV"}.
#' It must be at least 2 and leave at least 2 observations per validation fold.
#' This argument is used only when \code{crit = "CV"} and the parameter grid
#' contains more than one candidate.
#'
#' @param ebic.tuning A numeric value in [0,1] (default = 0.5) specifying
#' the tuning parameter to calculate for \code{crit = "EBIC"}.
#' This argument is used only when \code{crit = "EBIC"} and the parameter grid
#' contains more than one candidate.
#'
#' @return
#' An object with S3 class \code{"grasps"} containing the following components:
#' \describe{
#' \item{hatOmega}{The estimated precision matrix.}
#' \item{lambda}{The selected regularization parameter.}
#' \item{alpha}{The selected mixing parameter.}
#' \item{initial}{The initial estimate used to construct penalty weights
#' when \code{penalty != "lasso"}.}
#' \item{gamma}{The additional parameter used for the selected penalty when
#' \code{penalty != "lasso"}.}
#' \item{iterations}{For \code{penalty = "lasso"}, the number of ADMM iterations.
#' For other penalties, the total number of ADMM iterations across all LLA
#' subproblems, excluding any ADMM iterations used to construct the initial
#' estimate.}
#' \item{converged}{A logical value indicating whether all ADMM subproblems and,
#' when \code{initial = "lasso"}, the initial lasso fit satisfied their
#' convergence criteria before reaching \code{maxiter}.}
#' \item{initial.converged}{For \code{penalty != "lasso"}, a logical value
#' indicating whether the initial lasso fit converged.
#' It is \code{NA} when \code{initial = "invS"} or a matrix.}
#' \item{lla.iterations}{The number of LLA steps performed for non-convex
#' penalties. It is \code{1} for one-step LLA and \code{NA} for \code{"lasso"}
#' and \code{"adapt"}.}
#' \item{lla.converged}{For iterative LLA, a logical value indicating whether
#' the stopping criterion was satisfied before reaching \code{maxiter.lla}.
#' It is \code{NA} for one-step LLA, \code{penalty = "lasso"}, and
#' \code{penalty = "adapt"}.}
#' \item{lambda.safe}{The conservative upper bound used to initialize
#' the line search for \eqn{\lambda_{\max}}, corresponding to each distinct
#' value in \code{alpha.grid}, when \code{lambda = NULL}.}
#' \item{loss}{The minimum mean k-fold loss when \code{crit = "CV"}.}
#' \item{CV.loss}{Matrix of CV losses, with rows for parameter combinations and
#' columns for CV folds, when \code{crit = "CV"}.}
#' \item{score}{The optimal information criterion score when \code{crit} is set
#' to \code{"AIC"}, \code{"BIC"}, \code{"EBIC"}, or \code{"HBIC"}.}
#' \item{IC.score}{The information criterion score for each parameter
#' combination when \code{crit} is set to \code{"AIC"}, \code{"BIC"},
#' \code{"EBIC"}, or \code{"HBIC"}.}
#' \item{converged.grid}{The convergence status of the candidate fits.
#' For information criterion selection, this is a logical vector corresponding
#' to \code{lambda.grid} and \code{alpha.grid}.
#' For cross-validation, this is a logical matrix with parameter combinations
#' in rows and folds in columns.
#' All candidates are retained for model selection regardless of convergence
#' status.
#' Non-finite scores and estimates that are not positive definite are assigned
#' infinite scores.
#' Selection stops with an error if no candidate has a finite score.
#' A warning is issued if the selected fit or its selected CV fits did not
#' satisfy the ADMM convergence criteria.}
#' \item{lla.converged.grid}{The iterative LLA convergence status of
#' the candidate fits.
#' For information criterion selection, this is a logical vector corresponding
#' to \code{lambda.grid} and \code{alpha.grid}.
#' For cross-validation, this is a logical matrix with parameter combinations
#' in rows and folds in columns.
#' Entries are \code{NA} when iterative LLA is not used.
#' A warning is issued if the selected fit or its selected CV fits did not
#' satisfy the iterative LLA convergence criterion.}
#' \item{lambda.grid}{The actual lambda grid used in the program.}
#' \item{alpha.grid}{The actual alpha grid used in the program.}
#' \item{membership}{The group membership.}
#' }
#'
#' @details
#' For non-convex penalties, the default \code{lla = FALSE} computes a one-step
#' LLA estimator.
#' Penalty weights are constructed from the estimator specified by \code{initial}
#' and held fixed while the resulting weighted subproblem is solved.
#'
#' When \code{lla = TRUE}, the penalty weights are updated after each weighted
#' subproblem until the relative change between successive precision matrix
#' estimates is no greater than \code{tol.lla}, or \code{maxiter.lla} is reached.
#' This iterative procedure is intended to approach a stationary point of
#' the original non-convex objective but does not guarantee a global optimum.
#'
#' For \code{penalty = "adapt"}, weights are constructed from the estimator
#' specified by \code{initial} and remain fixed; therefore, \code{lla} is
#' ignored.
#' The \code{converged} component concerns convergence of the fixed-weight ADMM
#' subproblems and, when \code{initial = "lasso"}, the initial lasso fit.
#' It does not indicate stationarity of the original non-convex objective.
#' The \code{lla.converged} component concerns the outer iterative LLA procedure.
#'
#' @references
#' \insertAllCited{}
#'
#' @example
#' inst/example/ex-grasps.R
#'
#' @importFrom stats cov
#' @importFrom Rdpack reprompt
#'
#' @export

grasps <- function(X, n = nrow(X), membership, penalty,
                   diag.ind = TRUE, diag.grp = TRUE, diag.include = FALSE,
                   lambda = NULL, alpha = NULL, gamma = NULL, initial = "lasso",
                   nlambda = 10, lambda.min.ratio = 0.01,
                   growiter.lambda = 30, tol.lambda = 1e-03, maxiter.lambda = 50,
                   lla = FALSE, tol.lla = 1e-04, maxiter.lla = 10,
                   rho = 2, tau.incr = 2, tau.decr = 2, nu = 10,
                   tol.abs = 1e-04, tol.rel = 1e-04, maxiter = 1e+04,
                   crit = "BIC", kfold = 5, ebic.tuning = 0.5) {

  ## Validate the input matrix and sample size ---------------------------------

  if (!is.matrix(X) || !is.numeric(X) || is.complex(X) ||
      nrow(X) < 1L || ncol(X) < 1L || any(!is.finite(X))) {
    stop('`X` must be a nonempty finite numeric matrix.', call. = FALSE)
  }

  cov.input <- isSymmetric(X)

  if (cov.input && missing(n)) {
    stop('`n` must be supplied when `X` is a covariance matrix.', call. = FALSE)
  }

  validate_scalar(n, "n", lower = 2, integer = TRUE, lower.inclusive = TRUE)

  if (cov.input) {
    S <- X
    X <- NULL
  } else {
    if (nrow(X) < 2L || n != nrow(X)) {
      stop('For data input, `n` must equal `nrow(X)`, which must be at least 2.',
           call. = FALSE)
    }
    S <- (n-1)/n*cov(X)
  }

  if (any(!is.finite(S)) || any(diag(S) <= 0)) {
    stop('The covariance matrix must be finite with positive diagonal entries.',
         call. = FALSE)
  }

  if (cov.input) {
    eigvals <- eigen(S, symmetric = TRUE, only.values = TRUE)$values
    tol <- sqrt(.Machine$double.eps) * max(abs(eigvals))
    if (min(eigvals) < -tol) {
      stop('The covariance matrix must be positive semidefinite.', call. = FALSE)
    }
  }

  p <- ncol(S)

  ## Validate group membership -------------------------------------------------

  if (!is.numeric(membership) || is.complex(membership) ||
      !is.null(dim(membership)) || any(!is.finite(membership)) ||
      any(membership != floor(membership))) {
    stop('`membership` must be a vector of finite integer group labels.',
         call. = FALSE)
  }

  if (length(membership) != p) {
    stop('The length of `membership` must equal `p`.', call. = FALSE)
  }

  ## Validate the penalty specification ----------------------------------------

  penalties <- c("lasso", "adapt", "atan", "exp", "lq", "lsp", "mcp", "scad")

  if (!is.character(penalty) || length(penalty) != 1L ||
      is.na(penalty) || !(penalty %in% penalties)) {
    stop(sprintf("Invalid `penalty`.\nAvailable options: %s.",
                 paste(sprintf('"%s"', penalties), collapse = ", ")),
         call. = FALSE)
  }

  ## Validate diagonal-penalty controls ----------------------------------------

  diags <- list(diag.ind = diag.ind,
                diag.grp = diag.grp,
                diag.include = diag.include)

  for (name in names(diags)) {
    value <- diags[[name]]
    if (!is.logical(value) || length(value) != 1L || is.na(value)) {
      stop(sprintf('`%s` must be TRUE or FALSE.', name), call. = FALSE)
    }
  }

  ## Validate tuning parameters ------------------------------------------------

  lambda_null <- is.null(lambda)

  if (!lambda_null &&
      (!is.numeric(lambda) || is.complex(lambda) || length(lambda) == 0 ||
       any(!is.finite(lambda)) || any(lambda <= 0))) {
    stop('`lambda` must contain positive finite numeric values.', call. = FALSE)
  }

  if (lambda_null) {
    validate_scalar(nlambda, "nlambda", integer = TRUE)
    validate_scalar(lambda.min.ratio, "lambda.min.ratio", upper = 1)
    validate_scalar(growiter.lambda, "growiter.lambda", integer = TRUE,
                    lower.inclusive = TRUE)
    validate_scalar(tol.lambda, "tol.lambda")
    validate_scalar(maxiter.lambda, "maxiter.lambda", integer = TRUE)
  }

  alpha_null <- is.null(alpha)

  if (!alpha_null &&
      (!is.numeric(alpha) || is.complex(alpha) || length(alpha) == 0 ||
       any(!is.finite(alpha)) || any(alpha < 0 | alpha > 1))) {
    stop('`alpha` must contain finite numeric values in [0,1].', call. = FALSE)
  }

  if (penalty == "lasso") {
    gamma <- NA_real_
  } else {
    if (is.null(gamma)) {
      gamma <- switch(penalty,
                      "adapt" = 0.5, "atan" = 0.005, "exp"  = 0.01, "lq" = 0.5,
                      "lsp"  = 0.1, "mcp"  = 3, "scad" = 3.7)
    }
    lower <- switch(penalty, "mcp" = 1, "scad" = 2, 0)
    validate_scalar(gamma, "gamma", lower = lower)
    if (penalty == "lq" && gamma >= 1) {
      stop('For "lq", `gamma` must be in (0,1).', call. = FALSE)
    }

    valid.initial <- if (is.character(initial)) {
      length(initial) == 1L && !is.na(initial) &&
        initial %in% c("lasso", "invS")
    } else {
      is.matrix(initial) && is.numeric(initial) && !is.complex(initial) &&
        identical(dim(initial), c(p, p)) && all(is.finite(initial)) &&
        isSymmetric(initial)
    }
    if (!valid.initial) {
      stop('`initial` must be "lasso", "invS", or a finite symmetric p-by-p ',
           'numeric matrix.', call. = FALSE)
    }
  }

  if (penalty != "lasso" && identical(initial, "invS") && n <= p) {
    stop('`initial = "invS"` requires n > p.', call. = FALSE)
  }

  ## Validate LLA controls -----------------------------------------------------

  if (!is.logical(lla) || length(lla) != 1L || is.na(lla)) {
    stop('`lla` must be TRUE or FALSE.', call. = FALSE)
  }

  if (lla && !(penalty %in% c("lasso", "adapt"))) {
    validate_scalar(tol.lla, "tol.lla")
    validate_scalar(maxiter.lla, "maxiter.lla", integer = TRUE)
  }

  ## Validate ADMM controls ----------------------------------------------------

  validate_scalar(rho, "rho")
  validate_scalar(tau.incr, "tau.incr", lower = 1)
  validate_scalar(tau.decr, "tau.decr", lower = 1)
  validate_scalar(nu, "nu", lower = 1)
  validate_scalar(tol.abs, "tol.abs")
  validate_scalar(tol.rel, "tol.rel")
  validate_scalar(maxiter, "maxiter", integer = TRUE)

  ## Construct the alpha grid --------------------------------------------------

  if (alpha_null) {
    alpha <- seq(0.1, 0.9, 0.1)
  }

  ## Validate model-selection controls -----------------------------------------

  ncandidate <- length(alpha) *
    if (lambda_null) nlambda else length(lambda)

  if (ncandidate > 1L) {

    criteria <- c("AIC", "BIC", "EBIC", "HBIC", "CV")

    if (!is.character(crit) || length(crit) != 1L ||
        is.na(crit) || !(crit %in% criteria)) {
      stop(sprintf("Invalid `crit`.\nAvailable options: %s.",
                   paste(sprintf('"%s"', criteria), collapse = ", ")),
           call. = FALSE)
    }

    if (crit == "EBIC") {
      validate_scalar(ebic.tuning, "ebic.tuning", upper = 1,
                      lower.inclusive = TRUE)
    }

    if (crit == "CV") {
      if (cov.input) {
        stop('CV requires the n-by-p data matrix.', call. = FALSE)
      }
      validate_scalar(kfold, "kfold", lower = 2, integer = TRUE,
                      lower.inclusive = TRUE)
      if (kfold > floor(n / 2)) {
        stop('`kfold` must leave at least 2 observations in every validation fold.',
             call. = FALSE)
      }
      if (penalty != "lasso" && identical(initial, "invS") &&
          n - ceiling(n / kfold) <= p) {
        stop('`initial = "invS"` requires more training observations than ',
             'variables in every cross-validation fold.', call. = FALSE)
      }
      if (penalty != "lasso" && is.matrix(initial)) {
        stop('A user-supplied `initial` matrix cannot be used when ',
             '`crit = "CV"`.', call. = FALSE)
      }
    }
  }

  ## Construct group indices ---------------------------------------------------

  grp <- sort(unique(membership))
  group.idx <- lapply(grp, function(g) which(membership == g) - 1)
  names(group.idx) <- grp

  ## Construct the lambda grid -------------------------------------------------

  if (lambda_null) {
    max_abs.off <- max(abs(S[upper.tri(S, diag = FALSE)]), 0)
    S2 <- S*S
    sum.row <- rowsum(S2, membership)
    block.sum <- rowsum(t(sum.row), membership)
    block.norm <- sqrt(block.sum)
    max_block.norm.off <- max(block.norm[upper.tri(block.norm, diag = FALSE)], 0)
    idx_list <- split(seq_along(membership), membership)
    block.off_list <- list()
    if (length(idx_list) >= 2L) {
      for (g in seq_along(idx_list)) {
        for (gp in seq_along(idx_list)) {
          if (g < gp) {
            block.off_list[[length(block.off_list) + 1L]] <-
              S[idx_list[[g]], idx_list[[gp]], drop = FALSE]
          }
        }
      }
    }
    lambda_path <- lapply(alpha, function(a) {
      lambda.ind <- ifelse(a > 0, max_abs.off / a, 0)
      lambda.grp <- ifelse(a < 1, max_block.norm.off / (1 - a), 0)
      lambda.safe <- max(lambda.ind, lambda.grp)
      if (a > 0 && a < 1 && length(block.off_list) > 0L) {
        lambda.max <- line_search_lambda_max(
          blockList = block.off_list, lambda.safe = lambda.safe, alpha = a,
          growiter = growiter.lambda, tol = tol.lambda, maxiter = maxiter.lambda)
      } else {
        lambda.max <- lambda.safe
      }
      if (!is.finite(lambda.max) || lambda.max <= 0) {
        stop("Cannot generate a positive `lambda` sequence for alpha = ", a,
             ": the computed `lambda.max` is ", lambda.max,
             ".\nPlease supply a positive `lambda` explicitly.", call. = FALSE)
      }
      lambda.min <- lambda.min.ratio * lambda.max
      return(list(lambda = exp(seq(log(lambda.max), log(lambda.min),
                                   length = nlambda)),
                  lambda.safe = lambda.safe))
    })
    lambda.safe <- sapply(lambda_path, `[[`, "lambda.safe")
    names(lambda.safe) <- alpha
  }

  ## Construct the parameter grid ----------------------------------------------

  if (lambda_null) {
    parameter <- do.call(rbind, lapply(seq_along(alpha), function(k) {
      data.frame(alpha = rep(alpha[k], nlambda),
                 lambda = lambda_path[[k]]$lambda)
    }))
    rownames(parameter) <- NULL
  } else {
    parameter <- data.frame(alpha = rep(alpha, each = length(lambda)),
                            lambda = rep(lambda, times = length(alpha)))
  }

  ## Fit the model -------------------------------------------------------------

  if (ncandidate > 1L) {
    if (crit == "CV") {
      CV <- ADMMsggm_CV(
        X = X, group_idx = group.idx, penalty = penalty,
        diag_ind = diag.ind, diag_grp = diag.grp, diag_include = diag.include,
        lambdas = parameter$lambda, alphas = parameter$alpha, gamma = gamma,
        initial = initial,
        lla = lla, tol_lla = tol.lla, maxiter_lla = maxiter.lla,
        rho = rho, tau_incr = tau.incr, tau_decr = tau.decr, nu = nu,
        tol_abs = tol.abs, tol_rel = tol.rel, maxiter = maxiter,
        kfold = kfold
      )
      result <- ADMMsggm(
        S = S, group_idx = group.idx, penalty = penalty,
        diag_ind = diag.ind, diag_grp = diag.grp, diag_include = diag.include,
        lambda = CV$lambda.opt, alpha = CV$alpha.opt, gamma = gamma,
        initial = initial,
        lla = lla, tol_lla = tol.lla, maxiter_lla = maxiter.lla,
        rho = rho, tau_incr = tau.incr, tau_decr = tau.decr, nu = nu,
        tol_abs = tol.abs, tol_rel = tol.rel, maxiter = maxiter
      )
      if (lambda_null) {
        result$lambda.safe <- lambda.safe
      }
      result$loss <- CV$loss.min
      result$CV.loss <- CV$CV.loss
      result$converged.grid <- CV$converged.grid == 1
      result$lla.converged.grid <- CV$lla.converged.grid
    } else {
      IC <- ADMMsggm_IC(
        S = S, group_idx = group.idx, penalty = penalty,
        diag_ind = diag.ind, diag_grp = diag.grp, diag_include = diag.include,
        lambdas = parameter$lambda, alphas = parameter$alpha, gamma = gamma,
        initial = initial,
        lla = lla, tol_lla = tol.lla, maxiter_lla = maxiter.lla,
        rho = rho, tau_incr = tau.incr, tau_decr = tau.decr, nu = nu,
        tol_abs = tol.abs, tol_rel = tol.rel, maxiter = maxiter,
        crit = crit, n = n, ebic_tuning = ebic.tuning
      )
      result <- IC$result
      if (lambda_null) {
        result$lambda.safe <- lambda.safe
      }
      result$score <- IC$score.min
      result$IC.score <- as.vector(IC$IC.score)
      result$converged.grid <- as.vector(IC$converged.grid == 1)
      result$lla.converged.grid <- as.vector(IC$lla.converged.grid)
    }
  } else {
    result <- ADMMsggm(
      S = S, group_idx = group.idx, penalty = penalty,
      diag_ind = diag.ind, diag_grp = diag.grp, diag_include = diag.include,
      lambda = parameter$lambda, alpha = parameter$alpha, gamma = gamma,
      initial = initial,
      lla = lla, tol_lla = tol.lla, maxiter_lla = maxiter.lla,
      rho = rho, tau_incr = tau.incr, tau_decr = tau.decr, nu = nu,
      tol_abs = tol.abs, tol_rel = tol.rel, maxiter = maxiter
    )
    if (lambda_null) {
      result$lambda.safe <- lambda.safe
    }
  }

  ## Store the parameter grid --------------------------------------------------

  result$lambda.grid <- parameter$lambda
  result$alpha.grid <- parameter$alpha

  ## Validate the selected estimate --------------------------------------------

  positive.definite <- isTRUE(tryCatch({
    chol(result$hatOmega)
    TRUE
  }, error = function(e) FALSE))

  if (any(!is.finite(result$hatOmega)) || !positive.definite) {
    stop('The selected fit did not produce a finite positive-definite ',
         'precision matrix.\n',
         'Try increasing `maxiter` or changing the tuning parameters.',
         call. = FALSE)
  }

  cv.failed <- FALSE
  cv.lla.failed <- FALSE

  if (!is.null(result$CV.loss)) {
    selected <- which.min(rowMeans(result$CV.loss))
    cv.failed <- !all(result$converged.grid[selected, ])
    cv.lla.failed <- any(result$lla.converged.grid[selected, ] == FALSE,
                         na.rm = TRUE)
  }

  if (isFALSE(result$converged) || cv.failed) {
    warning('The selected fit or its cross-validation fits did not satisfy ',
            'the ADMM convergence criteria.\n',
            'Try increasing `maxiter` and inspect `converged` and ',
            '`converged.grid`.', call. = FALSE)
  }

  if (isFALSE(result$lla.converged) || cv.lla.failed) {
    warning('The selected fit or its cross-validation fits did not satisfy ',
            'the iterative LLA convergence criterion.\n',
            'Try increasing `maxiter.lla` and inspect `lla.converged` and ',
            '`lla.converged.grid`.', call. = FALSE)
  }

  ## Finalize the result -------------------------------------------------------

  result$membership <- membership
  class(result) <- c("grasps", "blkmat")

  return(result)

}
