library(grasps)

## reproducibility for everything
set.seed(1234)

## block-structured precision matrix based on SBM
sim <- gen_prec_sbm(d = 60, K = 3,
                    within.prob = 0.5, between.prob = 0.05,
                    weight.dists = list("gamma", "unif"),
                    weight.paras = list(c(shape = 20, scale = 5), c(min = 0, max = 1)),
                    cond.target = 100)
## visualization
plot(sim)

## n-by-d data matrix
library(MASS)
X <- mvrnorm(n = 30, mu = rep(0, 60), Sigma = sim$Sigma)

## lasso, BIC
res <- grasps(X = X, membership = sim$membership, penalty = "lasso", crit = "BIC")

## visualization
plot(res)
