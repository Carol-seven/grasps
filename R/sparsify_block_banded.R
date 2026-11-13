#' Groupwise Block-Banded Sparsifier
#'
#' @description
#' Apply a groupwise block-banded mask to a symmetric matrix, keeping only
#' entries within specified group neighborhoods.
#'
#' @param mat A numeric symmetric \eqn{p \times p} matrix, with \eqn{p} > 1,
#' specifying the matrix to be masked.
#'
#' @param membership An integer vector specifying the group membership.
#' Its length must equal \eqn{p}.
#' The numerical differences between group labels determine which groups are
#' considered neighbors.
#'
#' @param neighbor.range A non-negative integer (default = 1) specifying
#' the neighbor range.
#' Groups whose labels differ by at most \code{neighbor.range} are retained
#' in the mask.
#'
#' @return An object with S3 class \code{"sparsify_block_banded"} containing
#' the following components:
#' \describe{
#' \item{Omega}{The block-banded matrix.}
#' \item{Sigma}{The inverse of \code{Omega} if \code{Omega} is positive definite;
#' otherwise, \code{NULL}.}
#' \item{sparsity}{Proportion of absent edges among all possible undirected
#' edges in \code{Omega}.}
#' \item{membership}{An integer vector specifying the group membership.}
#' }
#'
#' @details
#' Entry \eqn{(i, j)} is retained when the absolute difference between the group
#' labels of variables \eqn{i} and \eqn{j} does not exceed \code{neighbor.range}.
#'
#' Applying the block-banded mask does not necessarily preserve positive
#' definiteness.
#' If the masked matrix is not positive definite, the function issues a warning
#' and returns \code{Sigma = NULL}.
#'
#' @example
#' inst/example/ex-sparsify_block_banded.R
#'
#' @export

sparsify_block_banded <- function(mat, membership, neighbor.range = 1) {

  ## validate `mat`
  validate_symmat(mat, "mat", min.dim = 2)

  ## validate `membership`
  p <- ncol(mat)
  if (!is.numeric(membership) || is.complex(membership) ||
      length(membership) != p || anyNA(membership) ||
      any(!is.finite(membership)) || any(membership != floor(membership))) {
    stop(sprintf("`membership` must contain %d finite integer values.", p),
         call. = FALSE)
  }

  ## validate `neighbor.range`
  validate_scalar(neighbor.range, "neighbor.range", lower = 0,
                  integer = TRUE, lower.inclusive = TRUE)

  ## determine which entries to keep: membership within 'neighbor.range'
  mask <- abs(outer(membership, membership, `-`)) <= neighbor.range

  ## apply mask to matrix
  Omega <- mat * mask

  ## verify positive definiteness and compute the covariance matrix
  chol.Omega <- tryCatch(
    chol(Omega),
    error = function(e) NULL
  )
  if (!is.null(chol.Omega)) {
    Sigma <- chol2inv(chol.Omega)
  } else {
    Sigma <- NULL
    warning(paste("The block-banded matrix is not positive definite.",
                  "`Sigma` is returned as `NULL`."), call. = FALSE)
  }

  Omega_edge <- Omega[upper.tri(Omega)]
  sparsity <- sum(Omega_edge == 0) / length(Omega_edge)

  result <- list(Omega = Omega,
                 Sigma = Sigma,
                 sparsity = sparsity,
                 membership = membership)
  class(result) <- c("sparsify_block_banded", "blkmat")
  return(result)
}
