#' SBM-based Gaussian Data with a Block-Structured Precision Matrix
#'
#' @description
#' Generate multivariate Gaussian data \eqn{X \in \mathbb{R}^{n \times p}} whose
#' precision matrix exhibits a block (group) structure induced by a Stochastic
#' Block Model (SBM).
#'
#' @param n An integer specifying the number of sample size.
#'
#' @param p An integer specifying the number of variables (dimensions).
#'
#' @param seed An integer (default = 1) specifying the random seed for
#' reproducibility.
#'
#' @param block.sizes An integer vector (default = NULL) specifying the size of
#' each group. If \code{NULL}, the \code{p} variables are divided as evenly as
#' possible across \code{K} groups.
#'
#' @param K An integer (default = 3) specifying the number of groups.
#' Ignored if \code{block.sizes} is provided; then \code{K <- length(block.sizes)}.
#'
#' @param prob.mat A \code{K}-by-\code{K} symmetric matrix (default = NULL)
#' specifying the Bernoulli rates. Element (i,j) gives the probability of
#' creating an edge between vertices from groups i and j. If \code{NULL},
#' a matrix with \code{within.prob} on the diagonal and \code{between.prob}
#' off-diagonal is used.
#'
#' @param within.prob A scalar in [0,1] (default = 0.25) specifying the
#' probability of creating an edge between vertices within the same group.
#'
#' @param between.prob A scalar in [0,1] (default = 0.05) specifying the
#' probability of creating an edge between vertices from different groups.
#'
#' @param weight.mat A \code{p}-by-\code{p} symmetric matrix (default = NULL)
#' specifying the edge weights. If \code{NULL}, weights are generated block-wise
#' according to \code{weight.dists} and \code{weight.paras}.
#'
#' @param weight.dists A character vector (default = c("gamma", "unif"))
#' specifying the distribution family for each block of weights.
#' Its length can be: \itemize{
#' \item length = 1: same distribution for all blocks.
#' \item length = 2: first for within-group blocks, second for between-group
#' blocks.
#' \item length = \code{K + K(K-1)/2}: full specification for each block (see
#' Details).
#' }
#' Accepted distributions (base R functions in parentheses) include:
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
#' @param weight.paras A list (default = list(c(shape = 1e4, rate = 1e2),
#' c(min = 0, max = 5)) specifying the parameters associated with
#' \code{weight.dists}. It must follow the same length rules as
#' \code{weight.dists}. Each element should be a named vector or list
#' suitable for the corresponding base R generator.
#'
#' @param min.eig A scalar (default = 0.1) specifying the minimum eigenvalue
#' target for the precision matrix. A diagonal shift ensures positive
#' definiteness.
#'
#' @details
#' Within-group and between-group edges are generated independently according to
#' Bernoulli distributions specified by \code{prob.mat}, or by \code{within.prob}
#' and \code{between.prob} if \code{prob.mat} is not supplied.
#'
#' Conditional on the adjacency structure, edge weights are sampled block-wise
#' from distributions specified in \code{weight.dists} and \code{weight.paras}.
#' The length of \code{weight.dists} (and \code{weight.paras}) determines how
#' weight distributions are assigned: \itemize{
#' \item length = 1: same distribution for all blocks;
#' \item length = 2: first for within-group blocks, second for between-group blocks;
#' \item length = \eqn{K + K(K - 1)/2}: full specification for each block.
#' }
#'
#' Block indexing for weight specifications uses \code{K} within-group blocks
#' with indices \code{1, \dots, K}, followed by \code{K(K-1)/2} between-group
#' blocks ordered as \code{(1,2), (1,3), \dots, (1,K), (2,3), \dots, (K-1,K)}.
#'
#' The weighted adjacency matrix is then symmetrized and used as the precision
#' matrix \eqn{\Omega}. Since arbitrary block-structured weights may not be
#' positive definite, a diagonal adjustment is applied to ensure the minimum
#' eigenvalue of \eqn{\Omega} is at least \code{min.eig}. The covariance matrix
#' is then computed as \eqn{\Sigma = \Omega^{-1}}, and the data matrix \eqn{X}
#' is generated from a multivariate normal distribution
#' \eqn{X \sim \mathcal{N}(0, \Sigma)}.
#'
#' @importFrom igraph as_adjacency_matrix sample_sbm
#' @importFrom MASS mvrnorm
#'
#' @return
#' A list with the following components:
#' \describe{
#' \item{Omega}{The precision matrix with SBM block structure.}
#' \item{Sigma}{The covariance matrix, i.e., the inverse of \code{Omega}.}
#' \item{sparsity}{Proportion of zero entries in \code{Omega}.}
#' \item{X}{An n-by-p matrix of Gaussian observations sampled from
#' \eqn{\mathcal{N}(0, \Sigma)}.}
#' \item{membership}{An integer vector specifying the group membership.}
#' }
#'
#' @export

gen_sbm_data <- function(n, p, seed = 1,
                         block.sizes = NULL, K = 3,
                         prob.mat = NULL, within.prob = 0.25, between.prob = 0.05,
                         weight.mat = NULL,
                         weight.dists = c("gamma", "unif"),
                         weight.paras = list(c(shape = 1e4, rate = 1e2),
                                             c(min = 0, max = 5)),
                         min.eig = 0.1) {

  ## reproducibility for everything
  set.seed(seed)

  ## block sizes (allow p not divisible by K)
  if (is.null(block.sizes)) {
    block.sizes <- rep(floor(p / K), K)
    rem <- p - sum(block.sizes)
    if (rem > 0) {
      block.sizes[seq_len(rem)] <- block.sizes[seq_len(rem)] + 1
    }
  }
  stopifnot(sum(block.sizes) == p)
  K <- length(block.sizes)

  ## membership and indices
  membership <- rep(1:K, times = block.sizes)
  idx_list <- split(1:p, membership)

  ## probability matrix
  if (is.null(prob.mat)) {
    prob.mat <- matrix(between.prob, K, K)
    diag(prob.mat) <- within.prob
  } else {
    stopifnot(all(dim(prob.mat) == c(K, K)))
  }

  ## SBM adjacency (undirected, no loops)
  SBM_sim <- igraph::sample_sbm(n = p, pref.matrix = prob.mat,
                                block.sizes = block.sizes,
                                directed = FALSE, loops = FALSE)
  SBM_adj <- as.matrix(igraph::as_adjacency_matrix(SBM_sim, sparse = FALSE))
  # ## generate SBM adjacency (upper triangle once, then mirror)
  # # SBM_adj <- matrix(0, p, p)
  # for (i in 1:K) {
  #   rows <- idx_list[[i]]
  #   for (j in i:K) {
  #     cols <- idx_list[[j]]
  #     block <- matrix(rbinom(length(rows) * length(cols), 1, prob.mat[i,j]),
  #                     length(rows), length(cols))
  #     SBM_adj[rows, cols] <- block
  #     if (j != i) {
  #       SBM_adj[cols, rows] <- t(block)
  #     }
  #   }
  # }
  # diag(SBM_adj) <- 0

  ## build weights
  ## fill upper-triangular; mirror to keep symmetric.
  if (is.null(weight.mat)) {
    dists_expand <- expand_spec(weight.dists, K)
    paras_expand <- expand_spec(weight.paras, K)
    weight.mat <- matrix(0, p, p)
    for (i in 1:K) {
      rows <- idx_list[[i]]
      nr <- length(rows)
      weight.mat[rows, rows] <- draw_dist(dists_expand[[i]], paras_expand[[i]], nr*nr)
      if (i < K) {
        for (j in (i+1):K) {
          cols <- idx_list[[j]]
          nc <- length(cols)
          k <- K + ((i-1)*(2*K-i))/2 + (j-i)
          weight.mat[rows, cols] <- draw_dist(dists_expand[[k]], paras_expand[[k]], nr*nc)
          weight.mat[cols, rows] <- t(weight.mat[rows, cols])
        }
      }
    }
  }

  ## precision matrix
  Omega <- SBM_adj * weight.mat
  Omega <- (Omega + t(Omega)) / 2 ## symmetric; diag still 0

  ## ensure positive definiteness by shifting diagonal
  # add a scalar tau to the diagonal so that lambda_min(Omega + tau*I) = min.eig
  ## (shift all eigenvalues by the same tau)
  eigmin <- min(eigen(Omega, only.values = TRUE)$values)
  tau <- ifelse(eigmin < min.eig, min.eig - eigmin, 0)
  diag(Omega) <- diag(Omega) + tau

  ## covariance matrix
  Sigma <- solve(Omega)
  ## sample
  X <- MASS::mvrnorm(n = n, mu = rep(0, ncol(Omega)), Sigma = Sigma)

  return(list(Omega = Omega, Sigma = Sigma, sparsity = sum(Omega == 0) / length(Omega),
              X = X, membership = membership))
}

