#' Groupwise Block-Diagonal Sparsifier
#'
#' @description
#' Make a precision-like matrix block-diagonal according to group membership,
#' keeping only within-group blocks. Optionally, further reduce selected
#' within-group blocks to diagonal-only (identity-like) structure.
#'
#' @param mat A p-by-p precision-like matrix specifying the base matrix to be
#' masked.
#'
#' @param membership An integer vector specifying the group membership.
#' The length of \code{membership} must be consistent with the dimension p.
#'
#' @param group.diag (default = \code{NULL})
#' \enumerate{
#' \item \code{NULL}: Make \code{mat} block-diagonal, i.e., only keep all
#' within-group blocks.
#' \item An integer vector specifying which within-group blocks are further
#' reduced to diagonal-only structure (based on the block-diagonal form).
#' }
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
#' ## block-diagonalizel; within-group block for group 3 made diagonal-only.
#' res <- sparsify_block_diag(est$hatOmega, membership, group.diag = 3)
#'
#' @export

sparsify_block_diag <- function(mat, membership, group.diag = NULL) {

  p <- ncol(mat)
  if (length(membership) != p) {
    stop(sprintf("Length of 'membership' (%d) must equal the matrix dimension p (%d).",
                 length(membership), p))
  }

  ## block-diagonal mask: TRUE only when in same group
  mask <- outer(membership, membership, `==`)

  ## adjust within-group blocks
  if (!is.null(group.diag)) {
    block_idx <- split(seq_along(membership), membership)
    for (g in group.diag) {
      idx <- block_idx[[as.character(g)]]
      mask[idx, idx] <- FALSE
      diag(mask[idx, idx]) <- TRUE
    }
  }

  ## apply mask to matrix
  Omega <- mat * mask
  ## compute covariance
  Sigma <- solve(Omega)

  return(list(Omega = Omega, Sigma = Sigma,
              sparsity = sum(Omega == 0) / length(Omega),
              membership = membership))
}
