# grasps

## Groupwise Regularized Adaptive Sparse Precision Solution

The goal of **grasps** is to provide a collection of statistical methods
that incorporate both element-wise and group-wise penalties to estimate
a precision matrix, making them user-friendly and useful for researchers
and practitioners.

## Penalties

The package **grasps** provides functions to estimate precision matrices
using the following penalties:

1.  Adaptive lasso (Zou 2006; Fan et al. 2009)

2.  Lasso (Tibshirani 1996; Friedman et al. 2008)

3.  Minimax concave penalty (MCP) (Zhang 2010)

4.  Smoothly clipped absolute deviation (SCAD) (Fan and Li 2001; Fan et
    al. 2009)

## Installation

You can install the development version of **grasps** from
[GitHub](https://github.com/) with:

``` R
# install.packages("devtools")
devtools::install_github("Carol-seven/grasps")
```

## Example

``` R
library(grasps)

X <- matrix(rnorm(200), 10, 20)
membership <- c(rep(1,5), rep(2,5), rep(3,4), rep(4,6))

grasps(X, membership = membership, penalty = "lasso", crit = "BIC")
```

## Reference

Fan, Jianqing, Yang Feng, and Yichao Wu. 2009. “Network Exploration via
the Adaptive LASSO and SCAD Penalties.” *The Annals of Applied
Statistics* 3 (2): 521–41. <https://doi.org/10.1214/08-aoas215>.

Fan, Jianqing, and Runze Li. 2001. “Variable Selection via Nonconcave
Penalized Likelihood and Its Oracle Properties.” *Journal of the
American Statistical Association* 96 (456): 1348–60.
<https://doi.org/10.1198/016214501753382273>.

Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. 2008. “Sparse
Inverse Covariance Estimation with the Graphical Lasso.” *Biostatistics*
9 (3): 432–41. <https://doi.org/10.1093/biostatistics/kxm045>.

Tibshirani, Robert. 1996. “Regression Shrinkage and Selection via the
Lasso.” *Journal of the Royal Statistical Society: Series B
(Methodological)* 58 (1): 267–88.
<https://doi.org/10.1111/j.2517-6161.1996.tb02080.x>.

Zhang, Cun-Hui. 2010. “Nearly Unbiased Variable Selection Under Minimax
Concave Penalty.” *The Annals of Statistics* 38 (2): 894–942.
<https://doi.org/10.1214/09-AOS729>.

Zou, Hui. 2006. “The Adaptive Lasso and Its Oracle Properties.” *Journal
of the American Statistical Association* 101 (476): 1418–29.
<https://doi.org/10.1198/016214506000000735>.
