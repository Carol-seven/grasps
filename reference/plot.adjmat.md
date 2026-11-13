# Plot Function for S3 Class "adjmat"

Visualize an adjacency matrix as a heatmap. This function is shared by
objects returned from
[`prec_to_adj`](https://shiying-xiao.com/grasps/reference/prec_to_adj.md).

## Usage

``` r
# S3 method for class 'adjmat'
plot(x, ...)
```

## Arguments

- x:

  An object inheriting from S3 class `"adjmat"`, typically returned by
  [`prec_to_adj`](https://shiying-xiao.com/grasps/reference/prec_to_adj.md).

- ...:

  Additional arguments passed to
  [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html).

## Value

A heatmap of class `ggplot` showing the matrix entries. The plot title
also reports matrix dimension and sparsity.

## Examples

``` r
library(grasps)

## reproducibility for everything
set.seed(1234)

## block-structured precision matrix based on SBM
sim <- gen_prec_sbm(p = 30, K = 3,
                    within.prob = 0.25, between.prob = 0.05,
                    weight.dists = list("gamma", "unif"),
                    weight.paras = list(c(shape = 20, rate = 10),
                                        c(min = 0, max = 5)),
                    cond.target = 100)
## ground truth visualization
plot(sim)


## n-by-p data matrix
library(MASS)
X <- mvrnorm(n = 20, mu = rep(0, 30), Sigma = sim$Sigma)

## precision matrix: adaptive lasso; BIC
prec <- grasps(X = X, membership = sim$membership, penalty = "adapt", crit = "BIC")

## precision matrix visualization
plot(prec)


## performance
performance(hatOmega = prec$hatOmega, Omega = sim$Omega)
#>      measure    value
#> 1   sparsity   0.8253
#> 2  Frobenius  25.7397
#> 3         KL   9.3274
#> 4  quadratic 144.4687
#> 5   spectral  11.1326
#> 6         TP  25.0000
#> 7         TN 336.0000
#> 8         FP  51.0000
#> 9         FN  23.0000
#> 10       TPR   0.5208
#> 11       FPR   0.1318
#> 12        F1   0.4032
#> 13       MCC   0.3210

## adjacency matrix: diagonal = 0; raw partial correlations;
##                   no thresholding; weighted network
adj <- prec_to_adj(prec$hatOmega,
                   diag.zero = TRUE, absolute = FALSE,
                   threshold = NULL, weighted = TRUE)

## adjacency matrix visualization
plot(adj)
```
