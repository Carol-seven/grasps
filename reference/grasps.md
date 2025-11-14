# Groupwise Regularized Adaptive Sparse Precision Solution

Provide a collection of statistical methods that incorporate both
element-wise and group-wise penalties to estimate a precision matrix.

## Usage

``` r
grasps(
  X,
  n = nrow(X),
  membership,
  penalty,
  diag.ind = TRUE,
  diag.grp = TRUE,
  diag.include = FALSE,
  lambda = NULL,
  alpha = NULL,
  gamma = NULL,
  nlambda = 10,
  lambda.min.ratio = 0.01,
  growiter.lambda = 30,
  tol.lambda = 0.001,
  maxiter.lambda = 50,
  rho = 2,
  tau.incr = 2,
  tau.decr = 2,
  nu = 10,
  tol.abs = 1e-04,
  tol.rel = 1e-04,
  maxiter = 10000,
  crit = "BIC",
  kfold = 5,
  ebic.tuning = 0.5
)
```

## Arguments

- X:

  1.  An n-by-d data matrix with sample size n and dimension d.

  2.  A d-by-d sample covariance/correlation matrix with dimension d.

- n:

  An integer (default = `nrow(X)`) specifying the sample size. This is
  only required when the input matrix `X` is a d-by-d sample
  covariance/correlation matrix with dimension d.

- membership:

  An integer vector specifying the group membership. The length of
  `membership` must be consistent with the dimension d.

- penalty:

  A character string specifying the penalty for estimating precision
  matrix. Available options include:

  1.  "adapt": adaptive lasso (Zou 2006; Fan et al. 2009) .

  2.  "lasso": lasso (Tibshirani 1996; Friedman et al. 2008) .

  3.  "mcp": minimax concave penalty (Zhang 2010) .

  4.  "scad": smoothly clipped absolute deviation (Fan and Li 2001; Fan
      et al. 2009) .

- diag.ind:

  A boolean (default = TRUE) specifying whether to penalize the diagonal
  elements.

- diag.grp:

  A boolean (default = TRUE) specifying whether to penalize the
  within-group blocks.

- diag.include:

  A boolean (default = FALSE) specifying whether to include the diagonal
  entries in the penalty for within-group blocks when `diag.grp = TRUE`.

- lambda:

  A grid of non-negative scalars for the regularization parameter. The
  default is `NULL`, which generates its own `lambda` sequence based on
  `nlambda` and `lambda.min.ratio`.

- alpha:

  A grid of scalars in \[0,1\] specifying the parameter leveraging the
  element-wise individual L1 penalty and the block-wise group L2
  penalty. An alpha of 1 corresponds to the element penalty only; an
  alpha of 0 corresponds to the group penalty only. The default values
  is a sequence from 0.05 to 0.95 with increments of 0.05.

- gamma:

  A scalar specifying the hyperparameter for the chosen `penalty`.
  Default values:

  1.  "adapt": 0.5

  2.  "mcp": 3

  3.  "scad": 3.7

- nlambda:

  An integer (default = 10) specifying the number of `lambda` values to
  be generated when `lambda = NULL`.

- lambda.min.ratio:

  A scalar (default = 0.01) specifying the fraction of the maximum
  `lambda` value \\\lambda\_{max}\\ to generate the minimum `lambda`
  \\\lambda\_{min}\\. If `lambda = NULL`, the program automatically
  generates a `lambda` grid as a sequence of length `nlambda` in log
  scale, starting from \\\lambda\_{min}\\ to \\\lambda\_{max}\\.

- growiter.lambda:

  An integer (default = 30) specifying the maximum number of exponential
  growth steps during the initial search for an admissible upper bound
  \\\lambda\_{\max}\\.

- tol.lambda:

  A scalar (default = 1e-03) specifying the relative tolerance for the
  bisection stopping rule on the interval width.

- maxiter.lambda:

  An integer (default = 50) specifying the maximum number of bisection
  iterations in the line search for \\\lambda\_{\max}\\.

- rho:

  A scalar \> 0 (default = 2) specifying the ADMM augmented-Lagrangian
  penalty parameter (often called the ADMM step size). Larger values
  typically put more weight on enforcing the consensus constraints each
  iteration; smaller values yield more conservative updates.

- tau.incr:

  A scalar \> 1 (default = 2) specifying the multiplicative factor used
  to increase `rho` when the primal residual dominates the dual residual
  in ADMM.

- tau.decr:

  A scalar \> 1 (default = 2) specifying the multiplicative factor used
  to decrease `rho` when the dual residual dominates the primal residual
  in ADMM.

- nu:

  A scalar \> 1 (default = 10) controlling how aggressively `rho` is
  rescaled in the adaptive-`rho` scheme (residual balancing).

- tol.abs:

  A scalar \> 0 (default = 1e-04) specifying the absolute tolerance for
  ADMM stopping (applied to primal/dual residual norms).

- tol.rel:

  A scalar \> 0 (default = 1e-04) specifying the relative tolerance for
  ADMM stopping (applied to primal/dual residual norms).

- maxiter:

  An integer (default = 1e+04) specifying the maximum number of ADMM
  iterations.

- crit:

  A string (default = "BIC") specifying the parameter selection method
  to use. Available options include:

  1.  "AIC": Akaike information criterion (Akaike 1973) .

  2.  "BIC": Bayesian information criterion (Schwarz 1978) .

  3.  "EBIC": extended Bayesian information criterion (Foygel and
      Drton 2010) .

  4.  "HBIC": high dimensional Bayesian information criterion (Wang et
      al. 2013; Fan et al. 2017) .

  5.  "CV": k-fold cross validation with negative log-likelihood loss.

- kfold:

  An integer (default = 5) specifying the number of folds used for
  `crit = "CV"`.

- ebic.tuning:

  A scalar (default = 0.5) specifying the tuning parameter to calculate
  for `crit = "EBIC"`.

## Value

An object with S3 class "grasps" containing the following components:

- hatOmega:

  The estimated precision matrix.

- lambda:

  The optimal regularization parameter.

- alpha:

  The optimal penalty balancing parameter.

- initial:

  The initial estimate of `hatOmega` when `penalty` is set to `"adapt"`,
  `"mcp"`, or `"scad"`.

- gamma:

  The optimal hyperparameter when `penalty` is set to `"adapt"`,
  `"mcp"`, or `"scad"`.

- iterations:

  The number of ADMM iterations.

- lambda.grid:

  The actual lambda grid used in the program.

- alpha.grid:

  The actual alpha grid used in the program.

- lambda.safe:

  The bisection-refined upper bound \\\lambda\_{\max}\\, corresponding
  to `alpha.grid`, when `lambda = NULL`.

- loss:

  The optimal k-fold loss when `crit = "CV"`.

- CV.loss:

  Matrix of CV losses, with rows for parameter combinations and columns
  for CV folds, when `crit = "CV"`.

- score:

  The optimal information criterion score when `crit` is set to `"AIC"`,
  `"BIC"`, `"EBIC"`, or `"HBIC"`.

- IC.score:

  The information criterion score for each parameter combination when
  `crit` is set to `"AIC"`, `"BIC"`, `"EBIC"`, or `"HBIC"`.

- membership:

  The group membership.

## References

Akaike H (1973). “Information Theory and an Extension of the Maximum
Likelihood Principle.” In Petrov BN, Csáki F (eds.), *Second
International Symposium on Information Theory*, 267–281. Akad\\emiai
Kiad\\o, Budapest, Hungary.  
  
Fan J, Feng Y, Wu Y (2009). “Network Exploration via the Adaptive LASSO
and SCAD Penalties.” *The Annals of Applied Statistics*, **3**(2),
521–541. [doi:10.1214/08-aoas215](https://doi.org/10.1214/08-aoas215)
.  
  
Fan J, Li R (2001). “Variable Selection via Nonconcave Penalized
Likelihood and its Oracle Properties.” *Journal of the American
Statistical Association*, **96**(456), 1348–1360.
[doi:10.1198/016214501753382273](https://doi.org/10.1198/016214501753382273)
.  
  
Fan J, Liu H, Ning Y, Zou H (2017). “High Dimensional Semiparametric
Latent Graphical Model for Mixed Data.” *Journal of the Royal
Statistical Society Series B: Statistical Methodology*, **79**(2),
405–421. [doi:10.1111/rssb.12168](https://doi.org/10.1111/rssb.12168)
.  
  
Foygel R, Drton M (2010). “Extended Bayesian Information Criteria for
Gaussian Graphical Models.” In Lafferty J, Williams C, Shawe-Taylor J,
Zemel R, Culotta A (eds.), *Advances in Neural Information Processing
Systems 23 (NIPS 2010)*, 604–612.  
  
Friedman J, Hastie T, Tibshirani R (2008). “Sparse Inverse Covariance
Estimation with the Graphical Lasso.” *Biostatistics*, **9**(3),
432–441.
[doi:10.1093/biostatistics/kxm045](https://doi.org/10.1093/biostatistics/kxm045)
.  
  
Schwarz G (1978). “Estimating the Dimension of a Model.” *The Annals of
Statistics*, **6**(2), 461–464.
[doi:10.1214/aos/1176344136](https://doi.org/10.1214/aos/1176344136) .  
  
Tibshirani R (1996). “Regression Shrinkage and Selection via the Lasso.”
*Journal of the Royal Statistical Society: Series B (Methodological)*,
**58**(1), 267–288.
[doi:10.1111/j.2517-6161.1996.tb02080.x](https://doi.org/10.1111/j.2517-6161.1996.tb02080.x)
.  
  
Wang L, Kim Y, Li R (2013). “Calibrating Nonconvex Penalized Regression
in Ultra-High Dimension.” *The Annals of Statistics*, **41**(5),
2505–2536. [doi:10.1214/13-AOS1159](https://doi.org/10.1214/13-AOS1159)
.  
  
Zhang C (2010). “Nearly Unbiased Variable Selection under Minimax
Concave Penalty.” *The Annals of Statistics*, **38**(2), 894–942.
[doi:10.1214/09-AOS729](https://doi.org/10.1214/09-AOS729) .  
  
Zou H (2006). “The Adaptive Lasso and Its Oracle Properties.” *Journal
of the American Statistical Association*, **101**(476), 1418–1429.
[doi:10.1198/016214506000000735](https://doi.org/10.1198/016214506000000735)
.

## Examples

``` r
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
```
