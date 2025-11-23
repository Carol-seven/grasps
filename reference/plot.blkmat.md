# Plot Function for Block-Structured Precision Matrices (Visualize a Matrix with Group Boundaries)

Visualize a precision matrix as a heatmap with dashed boundary lines
separating group blocks. This function is shared by objects returned
from [`grasps`](https://shiying-xiao.com/grasps/reference/grasps.md),
[`gen_prec_sbm`](https://shiying-xiao.com/grasps/reference/gen_prec_sbm.md),
and
[`sparsify_block_banded`](https://shiying-xiao.com/grasps/reference/sparsify_block_banded.md),
all of which inherit from the S3 class `"blkmat"`.

## Usage

``` r
# S3 method for class 'blkmat'
plot(x, colors = NULL, ...)
```

## Arguments

- x:

  An object inheriting from S3 class `"blkmat"`, typically returned by
  [`grasps`](https://shiying-xiao.com/grasps/reference/grasps.md),
  [`gen_prec_sbm`](https://shiying-xiao.com/grasps/reference/gen_prec_sbm.md)
  or
  [`sparsify_block_banded`](https://shiying-xiao.com/grasps/reference/sparsify_block_banded.md).

- colors:

  A vector of colors specifying an n-color gradient scale for the fill
  aesthetics.

- ...:

  Additional arguments passed to
  [`ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html).

## Value

A `ggplot2` heatmap showing the matrix entries. Dashed lines indicate
group boundaries. The plot title also reports matrix dimension and
sparsity.

## Examples

``` r
library(grasps)

## reproducibility for everything
set.seed(1234)

## block-structured precision matrix based on SBM
sim <- gen_prec_sbm(d = 100, K = 5,
                    within.prob = 0.5, between.prob = 0.05,
                    weight.dists = list("gamma", "unif"),
                    weight.paras = list(c(shape = 20, scale = 5), c(min = 0, max = 1)),
                    cond.target = 100)

## visualization
plot(sim)
```
