
<!-- README.md is generated from README.Rmd. Please edit that file -->

# grasps <img src="man/figures/logo.png" align="right" alt="" width="150">

## Groupwise Regularized Adaptive Sparse Precision Solution

<!-- badges: start -->

[![GitHub R package version](https://img.shields.io/github/r-package/v/Carol-seven/grasps?label=R%20in%20dev&color=green)](https://github.com/Carol-seven/grasps/blob/main/DESCRIPTION)
[![GitHub last commit](https://img.shields.io/github/last-commit/Carol-seven/grasps)](https://github.com/Carol-seven/grasps/commits/main)
[![R-CMD-check](https://github.com/Carol-seven/grasps/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Carol-seven/grasps/actions/workflows/R-CMD-check.yaml)
[![GitHub License](https://img.shields.io/github/license/Carol-seven/grasps?color=blue)](https://github.com/Carol-seven/grasps/blob/main/LICENSE.md)
<!-- badges: end -->

The goal of **grasps** is to provide a collection of statistical methods that
incorporate both element-wise and group-wise penalties to estimate a precision
matrix, making them user-friendly and useful for researchers and practitioners.

## Penalties

The package **grasps** provides functions to estimate precision matrices using
the following penalties:

1.  Adaptive lasso (Zou 2006; Fan et al. 2009)

2.  Lasso (Tibshirani 1996; Friedman et al. 2008)

3.  Minimax concave penalty (MCP) (Zhang 2010)

4.  Smoothly clipped absolute deviation (SCAD) (Fan and Li 2001; Fan et al. 2009)

## Installation

You can install the development version of **grasps** from
[GitHub](https://github.com/) with:

    # install.packages("devtools")
    devtools::install_github("Carol-seven/grasps")

## Example

``` r
library(grasps)

## reproducibility for everything
set.seed(1234)

## block-structured precision matrix based on SBM
sim <- gen_prec_sbm(d = 60, K = 5,
                    within.prob = 0.5, between.prob = 0.05,
                    weight.dists = list("gamma", "unif"),
                    weight.paras = list(c(shape = 20, scale = 5), c(min = 0, max = 1)),
                    cond.target = 100)

## synthetic data
library(MASS)
X <- MASS::mvrnorm(n = 30, mu = rep(0, 60), Sigma = sim$Sigma)

## solution
res <- grasps(X = X, membership = sim$membership, penalty = "lasso", crit = "BIC")

## visualization
plot(res)
```

<img src="man/figures/README-unnamed-chunk-2-1.png" width="100%" />

## Reference

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-fan2009network" class="csl-entry">

Fan, Jianqing, Yang Feng, and Yichao Wu. 2009. “Network Exploration via the Adaptive LASSO and SCAD Penalties.” *The Annals of Applied Statistics* 3 (2): 521–41. <https://doi.org/10.1214/08-aoas215>.

</div>

<div id="ref-fan2001variable" class="csl-entry">

Fan, Jianqing, and Runze Li. 2001. “Variable Selection via Nonconcave Penalized Likelihood and Its Oracle Properties.” *Journal of the American Statistical Association* 96 (456): 1348–60. <https://doi.org/10.1198/016214501753382273>.

</div>

<div id="ref-friedman2008sparse" class="csl-entry">

Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. 2008. “Sparse Inverse Covariance Estimation with the Graphical Lasso.” *Biostatistics* 9 (3): 432–41. <https://doi.org/10.1093/biostatistics/kxm045>.

</div>

<div id="ref-tibshirani1996regression" class="csl-entry">

Tibshirani, Robert. 1996. “Regression Shrinkage and Selection via the Lasso.” *Journal of the Royal Statistical Society: Series B (Methodological)* 58 (1): 267–88. <https://doi.org/10.1111/j.2517-6161.1996.tb02080.x>.

</div>

<div id="ref-zhang2010nearly" class="csl-entry">

Zhang, Cun-Hui. 2010. “Nearly Unbiased Variable Selection Under Minimax Concave Penalty.” *The Annals of Statistics* 38 (2): 894–942. <https://doi.org/10.1214/09-AOS729>.

</div>

<div id="ref-zou2006adaptive" class="csl-entry">

Zou, Hui. 2006. “The Adaptive Lasso and Its Oracle Properties.” *Journal of the American Statistical Association* 101 (476): 1418–29. <https://doi.org/10.1198/016214506000000735>.

</div>

</div>
