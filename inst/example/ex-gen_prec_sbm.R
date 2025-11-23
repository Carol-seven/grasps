library(grasps)

## reproducibility for everything
set.seed(1234)

## block-structured precision matrix based on SBM
#### case 1: base R distribution
sim1 <- gen_prec_sbm(d = 100, K = 5,
                     within.prob = 0.5, between.prob = 0.05,
                     weight.dists = list("gamma", "unif"),
                     weight.paras = list(c(shape = 20, scale = 5),
                                         c(min = 0, max = 1)),
                     cond.target = 100)
#### visualization
plot(sim1)

#### case 2: user-defined sampler
my_gamma <- function(n) {
  rgamma(n, shape = 10, scale = 5)
}
sim2 <- gen_prec_sbm(d = 100, K = 5,
                     within.prob = 0.5, between.prob = 0.05,
                     weight.dists = list(my_gamma, "unif"),
                     weight.paras = list(NULL,
                                         c(min = 0, max = 1)),
                     cond.target = 100)
#### visualization
plot(sim2)
