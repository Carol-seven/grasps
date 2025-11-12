#' Groupwise Block-Banded Sparsifier
#'
#' @description
#' Make a precision-like matrix block-banded according to group membership,
#' keeping only entries within specified group neighborhoods.
#'
#' @param mat A p-by-p precision-like matrix specifying the base matrix to be
#' masked.
#'
#' @param membership An integer vector specifying the group membership.
#' The length of \code{membership} must be consistent with the dimension p.
#'
#' @param neighbor.range An integer (default = 1) specifying the neighbor range,
#' where groups whose labels differ by at most \code{neighbor.range} are
#' considered neighbors and kept in the mask.
#'
#' @return A list containing:
#' \describe{
#' \item{Omega}{The masked precision matrix.}
#' \item{Sigma}{The covariance matrix, i.e., the inverse of \code{Omega}.}
#' \item{sparsity}{Proportion of zero entries in \code{Omega}.}
#' \item{membership}{An integer vector specifying the group membership.}
#' }
#'
#' @examples
#' ## reproducibility for everything
#' set.seed(1234)
#'
#' ## precision matrix estimation
#' X <- matrix(rnorm(200), 10, 20)
#' membership <- c(rep(1,5), rep(2,5), rep(3,4), rep(4,6))
#' est <- grasps(X, membership = membership, penalty = "lasso", crit = "BIC")
#'
#' ## default: keep blocks within ±1 of each group
#' res1 <- sparsify_block_banded(est$hatOmega, membership, neighbor.range = 1)
#' ## visualization
#' visualize(res1$Omega, res1$membership)
#'
#' ## wider band: keep blocks within ±2 of each group
#' res2 <- sparsify_block_banded(est$hatOmega, membership, neighbor.range = 2)
#' ## visualization
#' visualize(res2$Omega, res2$membership)
#'
#' @export

sparsify_block_banded <- function(mat, membership, neighbor.range = 1) {

  p <- ncol(mat)
  if (length(membership) != p) {
    stop(sprintf("Length of 'membership' (%d) must equal the matrix dimension p (%d).",
                 length(membership), p))
  }

  ## determine which entries to keep: membership within 'neighbor.range'
  mask <- abs(outer(membership, membership, `-`)) <= neighbor.range

  ## apply mask to matrix
  Omega <- mat * mask
  ## compute covariance
  Sigma <- solve(Omega)

  return(list(Omega = Omega, Sigma = Sigma,
              sparsity = sum(Omega == 0) / length(Omega),
              membership = membership))
}
