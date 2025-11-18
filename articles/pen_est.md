# Penalized Precision Matrix Estimation in grasps

## Preliminary

Consider the following setting:

- Gaussian graphical model (GGM) assumption:  
  The data $X_{n \times d}$ consists of independent and identically
  distributed samples $X_{1},\ldots,X_{n} \sim N_{d}(\mu,\Sigma)$.

- Disjoint group structure:  
  The $d$ variables can be partitioned into disjoint groups.

- Goal:  
  Estimate the precision matrix
  $\Omega = \Sigma^{- 1} = \left( \omega_{ij} \right)_{d \times d}$.

## Sparse-Group Estimator

\$\$\begin{align} \hat{\Omega}(\lambda,\alpha,\gamma) =
\operatorname\*{arg\\min}\_{\Omega \succ 0} \Bigl\\ - \log\det(\Omega) +
\operatorname{tr}(S\Omega) + P\_{\lambda,\alpha,\gamma}(\Omega) \Bigr\\,
\\ P\_{\lambda,\alpha,\gamma}(\Omega) = \alpha
P^\text{individual}\_{\lambda,\gamma}(\Omega) + (1-\alpha)
P^\text{group}\_{\lambda,\gamma}(\Omega), \\
P^\text{individual}\_{\lambda,\gamma}(\Omega) = \sum\_{i,j}
p\_{\lambda,\gamma}(\vert\omega\_{ij}\vert), \\
P^\text{group}\_{\lambda,\gamma}(\Omega) = \sum\_{g,g^\prime}
p\_{\lambda,\gamma}(\Vert\Omega\_{gg^\prime}\Vert_F), \end{align}\$\$

where:

- $S = n^{- 1}\sum_{i = 1}^{n}\left( X_{i} - \bar{X} \right)\left( X_{i} - \bar{X} \right)^{\top}$
  is the empirical covariance matrix.

- $\lambda \geq 0$ is the global regularization parameter controlling
  overall shrinkage.

- $\alpha \in \lbrack 0,1\rbrack$ is the mixing parameter controlling
  the balance between element-wise and block-wise penalties.

- $\gamma$ is the additional parameter controlling the curvature and
  effective degree of nonconvexity of the penalty.

- $P_{\lambda,\alpha,\gamma}(\Omega)$ is a generic bi-level penalty
  template that can incorporate convex or non-convex regularizers while
  preserving the intrinsic group structure among variables.

- $P_{\lambda,\gamma}^{\text{individual}}(\Omega)$ is the element-wise
  individual penalty term.

- $P_{\lambda,\gamma}^{\text{group}}(\Omega)$ is the block-wise group
  penalty term.

- $p_{\lambda,\gamma}( \cdot )$ is a penalty function parameterized by
  $\lambda$ and $\gamma$.

- $\Omega_{gg^{\prime}}$ is the submatrix of $\Omega$ with the rows from
  group $g$ and columns from group $g^{\prime}$.

- The Frobenius norm $\|\Omega\|_{F}$ is defined as
  $\|\Omega\|_{F} = \left( \sum_{i,j}|\omega_{ij}|^{2} \right)^{1/2} = \left\lbrack \operatorname{tr}\left( \Omega^{\top}\Omega \right) \right\rbrack^{1/2}$.

**Note**: For convex penalties, the parameter $\gamma$ is not required,
and the penalty function $p_{\lambda,\gamma}( \cdot )$ simplifies to
$p_{\lambda}( \cdot )$.

## Penalties

1.  Lasso: Least absolute shrinkage and selection operator ([Tibshirani
    1996](#ref-tibshirani1996regression); [Friedman, Hastie, and
    Tibshirani 2008](#ref-friedman2008sparse))

$$p_{\lambda}\left( \omega_{ij} \right) = \lambda|\omega_{ij}|.$$

2.  Adaptive lasso ([Zou 2006](#ref-zou2006adaptive); [Fan, Feng, and Wu
    2009](#ref-fan2009network))

$$p_{\lambda,\gamma}\left( \omega_{ij} \right) = \lambda\frac{|\omega_{ij}|}{v_{ij}},$$
where
$V = \left( v_{ij} \right)_{d \times d} = \left( |{\widetilde{\omega}}_{ij}|^{\gamma} \right)_{d \times d}$
is a matrix of adaptive weights, and ${\widetilde{\omega}}_{ij}$ is the
initial estimate obtained using `penalty = "lasso"`.

3.  Atan: Arctangent type penalty ([Wang and Zhu
    2016](#ref-wang2016variable))

$$p_{\lambda,\gamma}\left( \omega_{ij} \right) = \lambda\left( \gamma + \frac{2}{\pi} \right)\arctan\left( \frac{|\omega_{ij}|}{\gamma} \right),\quad\gamma > 0.$$

4.  Exp: Exponential type penalty ([Wang, Fan, and Zhu
    2018](#ref-wang2018variable))

$$p_{\lambda,\gamma}\left( \omega_{ij} \right) = \lambda\left\lbrack 1 - \exp\left( - \frac{|\omega_{ij}|}{\gamma} \right) \right\rbrack,\quad\gamma > 0.$$

5.  Lq ([Frank and Friedman 1993](#ref-frank1993statistical); [Fu
    1998](#ref-fu1998penalized); [Fan and Li
    2001](#ref-fan2001variable))

$$p_{\lambda,\gamma}\left( \omega_{ij} \right) = \lambda|\omega_{ij}|^{\gamma},\quad 0 < \gamma < 1.$$

6.  LSP: Log-sum penalty ([Candès, Wakin, and Boyd
    2008](#ref-candes2008enhancing))

$$p_{\lambda,\gamma}\left( \omega_{ij} \right) = \lambda\log\left( 1 + \frac{|\omega_{ij}|}{\gamma} \right),\quad\gamma > 0.$$

7.  MCP: Minimax concave penalty ([Zhang 2010](#ref-zhang2010nearly))

$$p_{\lambda,\gamma}\left( \omega_{ij} \right) = \begin{cases}
{\lambda|\omega_{ij}| - \frac{\omega_{ij}^{2}}{2\gamma},} & {{\text{if}\mspace{6mu}}|\omega_{ij}| \leq \gamma\lambda,} \\
{\frac{1}{2}\gamma\lambda^{2},} & {{\text{if}\mspace{6mu}}|\omega_{ij}| > \gamma\lambda.}
\end{cases}\quad\gamma > 1.$$

8.  SCAD: Smoothly clipped absolute deviation ([Fan and Li
    2001](#ref-fan2001variable); [Fan, Feng, and Wu
    2009](#ref-fan2009network))

$$p_{\lambda,\gamma}\left( \omega_{ij} \right) = \begin{cases}
{\lambda|\omega_{ij}|} & {{\text{if}\mspace{6mu}}|\omega_{ij}| \leq \lambda,} \\
\frac{2\gamma\lambda|\omega_{ij}| - \omega_{ij}^{2} - \lambda^{2}}{2(\gamma - 1)} & {{\text{if}\mspace{6mu}}\lambda < |\omega_{ij}| < \gamma\lambda,} \\
\frac{\lambda^{2}(\gamma + 1)}{2} & {{\text{if}\mspace{6mu}}|\omega_{ij}| \geq \gamma\lambda.}
\end{cases}\quad\gamma > 2.$$

## Illustrative Visualization

Figure 1 illustrates a comparison of various penalty functions
$p(\omega)$ evaluated over a range of $\omega$ values. The main panel
(right) provides a wider view of the penalty functions’ behavior for
larger $|\omega|$, while the inset panel (left) magnifies the region
near zero $\lbrack - 1,1\rbrack$.

![Figure 1: Illustrative penalty
functions.](pen_est_files/figure-html/pen-1.png)

Figure 1: Illustrative penalty functions.

Figure 2 displays the derivative function $p^{\prime}(\omega)$
associated with a range of penalty types. The Lasso exhibits a constant
derivative, corresponding to uniform shrinkage. For MCP and SCAD, the
derivatives are piecewise: initially equal to the Lasso derivative, then
decreasing over an intermediate region, and eventually dropping to zero,
indicating that large $|\omega|$ receive no shrinkage. Other non-convex
penalties show smoothly diminishing derivatives as $|\omega|$ increases,
reflecting their tendency to shrink small $|\omega|$ strongly while
exerting little to no shrinkage on large ones.

![Figure 2: Illustrative penalty
derivatives.](pen_est_files/figure-html/unnamed-chunk-2-1.png)

Figure 2: Illustrative penalty derivatives.

## Reference

Candès, Emmanuel J., Michael B. Wakin, and Stephen P. Boyd. 2008.
“Enhancing Sparsity by Reweighted $\ell_{1}$ Minimization.” *Journal of
Fourier Analysis and Applications* 14 (5): 877–905.
<https://doi.org/10.1007/s00041-008-9045-x>.

Fan, Jianqing, Yang Feng, and Yichao Wu. 2009. “Network Exploration via
the Adaptive LASSO and SCAD Penalties.” *The Annals of Applied
Statistics* 3 (2): 521–41. <https://doi.org/10.1214/08-aoas215>.

Fan, Jianqing, and Runze Li. 2001. “Variable Selection via Nonconcave
Penalized Likelihood and Its Oracle Properties.” *Journal of the
American Statistical Association* 96 (456): 1348–60.
<https://doi.org/10.1198/016214501753382273>.

Frank, Lldiko E., and Jerome H. Friedman. 1993. “A Statistical View of
Some Chemometrics Regression Tools.” *Technometrics* 35 (2): 109–35.
<https://doi.org/10.1080/00401706.1993.10485033>.

Friedman, Jerome, Trevor Hastie, and Robert Tibshirani. 2008. “Sparse
Inverse Covariance Estimation with the Graphical Lasso.” *Biostatistics*
9 (3): 432–41. <https://doi.org/10.1093/biostatistics/kxm045>.

Fu, Wenjiang J. 1998. “Penalized Regressions: The Bridge Versus the
Lasso.” *Journal of Computational and Graphical Statistics* 7 (3):
397–416. <https://doi.org/10.1080/10618600.1998.10474784>.

Tibshirani, Robert. 1996. “Regression Shrinkage and Selection via the
Lasso.” *Journal of the Royal Statistical Society: Series B
(Methodological)* 58 (1): 267–88.
<https://doi.org/10.1111/j.2517-6161.1996.tb02080.x>.

Wang, Yanxin, Qibin Fan, and Li Zhu. 2018. “Variable Selection and
Estimation Using a Continuous Approximation to the $L_{0}$ Penalty.”
*Annals of the Institute of Statistical Mathematics* 70 (1): 191–214.
<https://doi.org/10.1007/s10463-016-0588-3>.

Wang, Yanxin, and Li Zhu. 2016. “Variable Selection and Parameter
Estimation with the Atan Regularization Method.” *Journal of Probability
and Statistics* 2016: 6495417. <https://doi.org/10.1155/2016/6495417>.

Zhang, Cun-Hui. 2010. “Nearly Unbiased Variable Selection Under Minimax
Concave Penalty.” *The Annals of Statistics* 38 (2): 894–942.
<https://doi.org/10.1214/09-AOS729>.

Zou, Hui. 2006. “The Adaptive Lasso and Its Oracle Properties.” *Journal
of the American Statistical Association* 101 (476): 1418–29.
<https://doi.org/10.1198/016214506000000735>.
