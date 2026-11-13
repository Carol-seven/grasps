# Groupwise Block-Banded Sparsifier

Apply a groupwise block-banded mask to a symmetric matrix, keeping only
entries within specified group neighborhoods.

## Usage

``` r
sparsify_block_banded(mat, membership, neighbor.range = 1)
```

## Arguments

- mat:

  A numeric symmetric \\p \times p\\ matrix, with \\p\\ \> 1, specifying
  the matrix to be masked.

- membership:

  An integer vector specifying the group membership. Its length must
  equal \\p\\. The numerical differences between group labels determine
  which groups are considered neighbors.

- neighbor.range:

  A non-negative integer (default = 1) specifying the neighbor range.
  Groups whose labels differ by at most `neighbor.range` are retained in
  the mask.

## Value

An object with S3 class `"sparsify_block_banded"` containing the
following components:

- Omega:

  The block-banded matrix.

- Sigma:

  The inverse of `Omega` if `Omega` is positive definite; otherwise,
  `NULL`.

- sparsity:

  Proportion of absent edges among all possible undirected edges in
  `Omega`.

- membership:

  An integer vector specifying the group membership.

## Details

Entry \\(i, j)\\ is retained when the absolute difference between the
group labels of variables \\i\\ and \\j\\ does not exceed
`neighbor.range`.

Applying the block-banded mask does not necessarily preserve positive
definiteness. If the masked matrix is not positive definite, the
function issues a warning and returns `Sigma = NULL`.

## Examples

``` r
library(grasps)

## reproducibility for everything
set.seed(1234)

## precision matrix estimation
X <- matrix(rnorm(200), 10, 20)
membership <- rep(1:5, each = 4)
est <- grasps(X, membership = membership, penalty = "lasso",
              lambda = 0.1, alpha = 0.5)

## default: keep blocks within ±1 of each group
res1 <- sparsify_block_banded(est$hatOmega, membership, neighbor.range = 1)
plot(res1)


## wider band: keep blocks within ±2 of each group
res2 <- sparsify_block_banded(est$hatOmega, membership, neighbor.range = 2)
plot(res2)


## special case: block-diagonal matrix
res3 <- sparsify_block_banded(est$hatOmega, membership, neighbor.range = 0)
plot(res3)
```
