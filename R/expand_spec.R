#' Expand Block-Level Specifications for SBM Weight Generation
#'
#' @description
#' Expand a user-provided specification for within-group and between-group block
#' weight-generation rules into a full list of length \code{K + K(K-1)/2}, where:
#' - the first \code{K} elements correspond to within-group blocks, and
#' - the remaining \code{K(K-1)/2} elements correspond to between-group blocks.
#'
#' @param spec A list specifying weight parameters. Acceptable forms:
#' \enumerate{
#' \item length = 1: use the same specification for all blocks.
#' \item length = 2: first for within-group blocks, second for between-group
#' blocks.
#' \item length = \code{K + K(K-1)/2}: full specification for each block.
#' }
#'
#' @param K An integer specifying the number of groups.
#'
#' @return
#' A list of length \code{K + K(K-1)/2}, where the first \code{K} elements
#' correspond to within-group blocks, and the remaining \code{K(K-1)/2}
#' correspond to between-group blocks.
#'
#' @noRd

expand_spec <- function(spec, K) {

  ## number of within-group and between-group blocks
  n_within <- K
  n_between <- K*(K-1)/2
  n_total <- n_within + n_between

  ## ensure the input is a list
  ## if it is a vector, convert to a list
  if (!is.list(spec)) {
    spec <- as.list(spec)
  }

  if (length(spec) == 1) {
    ## Case 1: one rule -> use the same specification for all blocks
    return(rep(list(spec[[1]]), n_total))

  } else if (length(spec) == 2) {
    ## Case 2: two rules -> 1st rule for all within-group blocks
    ##                      2nd rule for all between-group blocks
    return(c(rep(list(spec[[1]]), n_within),
             rep(list(spec[[2]]), n_between)))

  } else if (length(spec) == n_total) {
    ## Case 3: full specification explicitly provided
    return(spec)
  }

  ## otherwise: invalid specification length
  stop(paste0(
    "`weight.dists`/`weight.paras` length must be 1, 2, or K + K*(K-1)/2.\n",
    "For K = ", K, ", expected lengths: 1, 2, or ", n_total, "."
  ))
}

