#' Block-Band Matrix with Group Structure
#'
#' @description
#' Construct a masked precision matrix where entries are kept only within
#' specified group neighborhoods. Optionally generate Gaussian samples from
#' the resulting precision matrix.
#'
#' @param mat A p-by-p precision-like matrix specifying the base matrix to be
#' masked.
#'
#' @param membership An integer vector specifying the group membership.
#' The length of \code{membership} must be consistent with the dimension p.
#'
#' @param within.diag A boolean (default = TRUE) specifying whether within-group
#' blocks keep only diagonal entries (i.e., identity-like structure within each
#' group).
#'
#' @param neighbor.range An integer (default = 1) specifying the neighbor range,
#' where groups whose labels differ by at most \code{neighbor.range} are
#' considered neighbors and kept in the mask.
#'
#' @param group.diag An integer vector (default = NULL) specifying which
#' within-group blocks should keep diagonal-only structure when
#' \code{within.diag = FALSE}.
#'
#' @param n An integer (default = NULL) specifying the sample size for
#' generating Gaussian data from the resulting precision matrix. If \code{NULL},
#' no data are generated.
#'
#' @param seed An integer (default = 1) specifying the random seed for
#' reproducibility.
#'
#' @importFrom MASS mvrnorm
#'
#' @return A list containing:
#' \describe{
#' \item{Omega}{The masked precision matrix.}
#' \item{Sigma}{The covariance matrix, i.e., the inverse of \code{Omega}.}
#' \item{sparsity}{Proportion of zero entries in \code{Omega}.}
#' \item{X}{If \code{n} is not \code{NULL}, an \code{n}-by-\code{p} matrix of
#' Gaussian observations sampled from \eqn{\mathcal{N}(0, \Sigma)}.}
#' \item{membership}{An integer vector specifying the group membership.}
#' }
#'
#' @export

blockband <- function(mat, membership,
                      within.diag = TRUE, neighbor.range = 1, group.diag = NULL,
                      n, seed = 1) {

  stopifnot(length(membership) == ncol(mat))

  ## determine which entries to keep: membership within 'neighbor.range'
  mask <- abs(outer(membership, membership, `-`)) <= neighbor.range

  ## adjust within-group blocks
  if (within.diag) {
    ## within-group blocks -> diagonal matrix
    block_idx <- split(seq_along(membership), membership)
    for (idx in block_idx) {
      mask[idx, idx] <- FALSE
      diag(mask[idx, idx]) <- TRUE
    }
  } else {
    if (!is.null(group.diag)) {
      for (group_idx in group.diag) {
        idx <- which(membership == group_idx)
        mask[idx, idx] <- FALSE
        diag(mask[idx, idx]) <- TRUE
      }
    }
  }

  ## apply mask to matrix
  Omega <- mat * mask
  ## compute covariance
  Sigma <- solve(Omega)

  result <- list(Omega = Omega, Sigma = Sigma,
                 sparsity = sum(Omega == 0) / length(Omega),
                 membership = membership)

  if (!is.null(n)) {
    set.seed(seed)
    ## sample
    result$X <- MASS::mvrnorm(n = n, mu = rep(0, ncol(Omega)), Sigma = Sigma)
  }

  return(result)
}

