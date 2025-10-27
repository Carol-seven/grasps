# grasps <img src="man/figure/grasps.png" align="right" alt="" width="150">


[![GitHub R package version](https://img.shields.io/github/r-package/v/Carol-seven/grasps?label=R%20in%20dev&color=green)](https://github.com/Carol-seven/grasps/blob/main/DESCRIPTION)
[![GitHub last commit](https://img.shields.io/github/last-commit/Carol-seven/grasps)](https://github.com/Carol-seven/grasps/commits/main)
[![R-CMD-check](https://github.com/Carol-seven/grasps/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Carol-seven/grasps/actions/workflows/R-CMD-check.yaml)
[![GitHub License](https://img.shields.io/github/license/Carol-seven/grasps)](https://github.com/Carol-seven/grasps/blob/main/LICENSE.md)


The goal of **grasps** is to provide a collection of statistical methods that
incorporate both element-wise and group-wise penalties to estimate a precision
matrix, making them user-friendly and useful for researchers and practitioners.


## Installation


You can install the development version of **grasps** from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("Carol-seven/grasps")
```


## Example


``` r
library(grasps)

X <- matrix(rnorm(200), 10, 20)
groups <- c(rep(1,5), rep(2,5), rep(3,4), rep(4,6))

sggm(X, groups = groups, penalty = "lasso", crit = "BIC")
```

