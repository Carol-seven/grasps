# grasps 0.1.3

+ Fixed within-group local linear approximation (LLA) weights so that diagonal
  entries are excluded when `diag.include = FALSE`.
+ Revised adaptive group weights so that `gamma` is applied consistently
  when constructing Frobenius-norm-based weights.
+ Added optional iterative LLA for non-convex penalties through `lla`, `tol.lla`,
  and `maxiter.lla`; one-step LLA remains the default.
+ Added configurable initial estimators for adaptive and non-convex penalties
  through `initial`, including the lasso estimate (default), the inverse
  empirical covariance matrix, and a user-supplied symmetric matrix.
+ Added `iterations`, `converged`, `initial.converged`, `lla.iterations`, and
  `lla.converged` to report inner ADMM and outer LLA convergence information.
+ Added `converged.grid` and `lla.converged.grid` to retain convergence status
  for every tuning-parameter candidate and cross-validation fold.
+ Added warnings when the selected estimate, its selected cross-validation fits,
  or the iterative LLA procedure does not satisfy its convergence criterion.
+ Strengthened validation of data, covariance matrices, group labels, penalty
  parameters, initial estimators, and solver controls before calling the native
  solvers.
+ Required an explicit sample size for covariance-matrix input and strengthened
  validation of cross-validation folds and training-fold initial estimators.
+ Limited validation of model-selection and automatic-grid controls to cases
  where those controls are active.
+ Added positive-definiteness and finite-score checks during model selection.
  Invalid candidates receive infinite scores, and selection stops when no
  candidate has a finite score.
+ Strengthened validation in `gen_prec_sbm()` and corrected its reported
  sparsity to use undirected off-diagonal edges.
+ Strengthened validation in `sparsify_block_banded()`, corrected its reported
  sparsity, and added explicit handling of non-positive-definite results.
+ Strengthened validation in `performance()`, retained computable metrics for
  non-positive-definite inputs, and returned `NA` for undefined metrics.
+ Strengthened validation in `prec_to_adj()` and clarified the handling of
  diagonal entries, thresholding, and non-positive-definite inputs.
+ Added vectorized input validation and penalty-specific `gamma` checks to
  `compute_penalty()` and `compute_derivative()`.
+ Expanded standalone regression tests for input validation, convergence
  reporting, invalid candidate handling, simulation utilities, performance
  measures, adjacency conversion, and penalty calculations.

# grasps 0.1.2

+ Fixed the tuning-parameter grid so that all combinations of user-supplied
  `alpha` and `lambda` values are evaluated.
+ Fixed automatic tuning-parameter generation so that a separate lambda sequence
  is constructed for each value of `alpha` using its corresponding `lambda.max`.
+ Validation of `nlambda` and `lambda.min.ratio` is now performed only when
  `lambda = NULL`.
+ Added a check for a non-positive or non-finite `lambda.max` to prevent
  the generation of invalid lambda sequences.
+ Fixed within-group weight generation in `gen_prec_sbm()` so that each edge
  weight is sampled once and assigned symmetrically to the corresponding entries
  of the precision matrix.
+ Added input validation for `cond.target` and explicit handling of empty graphs
  in `gen_prec_sbm()`.
+ Clarified the roles of `lambda.safe` and `lambda.max` in the documentation of
  the automatic lambda-grid construction.
+ Expanded the documentation of the positive-definiteness and condition-number
  adjustment in `gen_prec_sbm()` and added relevant references.

# grasps 0.1.1

+ Added function `plot.adjmat()`.
+ Added function `prec_to_adj()`.
+ Use `d` instead of `p` to denote the dimension.
+ Revised the estimator expression in the vignette.
+ Added NEWS file to record the changelog.

# grasps 0.1.0

+ Initial release of the grasps package.
+ The grasps is a toolbox for precision matrix estimation with group structure.
+ Key features:
  - Unified regularization framework for sparse network learning.
  - Supports element-wise sparsity and group-wise shrinkage.
  - Includes both convex and non-convex penalties (e.g., adaptive lasso, SCAD, MCP).
  - Provides model selection tools (e.g., CV, EBIC, HBIC).
+ Designed for applications where variables exhibit grouped or modular
  relationships (e.g., brain networks, biological pathways).
