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
