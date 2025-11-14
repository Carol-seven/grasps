library(grasps)

## reproducibility for everything
set.seed(1234)

## user-defined sampler
#### Case 1
my_gamma <- function(n) {
  rgamma(n, shape = 10, scale = 0.5)
}
draw_sample(my_gamma, n = 10)
## Case 2
my_unif <- function(n, min, max) {
  runif(n, min = 0, max = 1)
}
draw_sample(my_unif, para = list(min = 1, max = 5), n = 10)

## base R distribution
#### Case 1
draw_sample("gamma", para = list(shape = 10, scale = 0.5), n = 10)
#### Case 2
draw_sample("unif", para = list(min = 1, max = 5), n = 10)
