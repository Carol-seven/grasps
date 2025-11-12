# Penalized Precision Matrix Estimation

## Preliminary

Consider the following setting:

- Gaussian graphical model (GGM) assumption:  
  The data $X_{n \times p}$ consists of independent and identically
  distributed samples $X_{1},\ldots,X_{n} \sim N_{p}(0,\Sigma)$.

- Disjoint group structure:  
  The $p$ variables can be partitioned into disjoint groups.

- Goal:  
  Estimate the precision matrix
  $\Omega = \Sigma^{- 1} = \left( \omega_{ij} \right)_{p \times p}$.

## Bi-level penalty

$$\widehat{\Omega} = \operatorname{arg\,min}\limits_{\Omega \succ 0}\left\{ - \log\det(\Omega) + \operatorname{tr}(S\Omega) + \alpha P_{\text{individual}} + (1 - \alpha)P_{\text{group}} \right\},$$
where:

- $\alpha \in \lbrack 0,1\rbrack$ controls the balance between
  element-wise and block-wise penalties.
- $P_{\text{individual}}$ denotes the element-wise individual penalty
  term.
- $P_{\text{group}}$ denotes the block-wise group penalty term.

## Penalties

The package **grasps** estimates precision matrices using the following
penalties:

1.  Adaptive lasso (Zou 2006; Fan, Feng, and Wu 2009)

$$P_{\text{individual}} = \lambda\sum\limits_{i,j}\frac{|\omega_{ij}|}{|v_{ij}|}\quad\text{and}\quad P_{\text{group}} = \lambda\sum\limits_{g,g^{\prime}}\frac{\|\Omega_{gg^{\prime}}\|_{2}}{\| V_{gg^{\prime}}\|_{2}}$$

2.  Lasso (Tibshirani 1996; Friedman, Hastie, and Tibshirani 2008)

$$P_{\text{individual}} = \lambda\|\Omega\|_{1}\quad\text{and}\quad P_{\text{group}} = \lambda\sum\limits_{g,g^{\prime}}\|\Omega_{gg^{\prime}}\|_{2}$$

3.  Minimax concave penalty (MCP) (Zhang 2010)

$$P_{\text{individual}} = \sum\limits_{i,j}\xi_{\lambda,\gamma}\left( |\omega_{ij}| \right)\quad\text{and}\quad P_{\text{group}} = \sum\limits_{g,g^{\prime}}\xi_{\lambda,\gamma}\left( \|\Omega_{gg^{\prime}}\|_{2} \right)$$

4.  Smoothly clipped absolute deviation (SCAD) (Fan and Li 2001; Fan,
    Feng, and Wu 2009)

$$P_{\text{individual}} = \sum\limits_{i,j}\psi_{\lambda,\gamma}\left( |\omega_{ij}| \right)\quad\text{and}\quad P_{\text{group}} = \sum\limits_{g,g^{\prime}}\psi_{\lambda,\gamma}\left( \|\Omega_{gg^{\prime}}\|_{2} \right)$$

where:

- $\Omega_{gg^{\prime}}$ denotes the submatrix of $\Omega$ with the rows
  from group $g$ and columns from group $g^{\prime}$.

- The norms are defined as
  $$\|\Omega\|_{1} = \sum\limits_{i,j}|\omega_{ij}|\quad\text{and}\quad\|\Omega\|_{2} = \|\Omega\|_{F} = \sqrt{\sum\limits_{i,j}|\omega_{ij}|^{2}} = \sqrt{\operatorname{tr}\left( \Omega^{\top}\Omega \right)}.$$

- $\lambda > 0$ is a regularization parameter.

- $V = \left( v_{ij} \right)_{p \times p}$ is a matrix of adaptive
  weights, which is the estimate from `penalty = "lasso"`.

- $\xi_{\lambda,\gamma}$ is the penalty function of MCP.

- $\psi_{\lambda,\gamma}$ is the penalty function of SCAD.

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
