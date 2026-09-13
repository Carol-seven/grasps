# Changelog

## grasps 0.1.2

CRAN release: 2026-09-12

- Fixed the tuning-parameter grid so that all combinations of
  user-supplied `alpha` and `lambda` values are evaluated.
- Fixed automatic tuning-parameter generation so that a separate lambda
  sequence is constructed for each value of `alpha` using its
  corresponding `lambda.max`.
- Validation of `nlambda` and `lambda.min.ratio` is now performed only
  when `lambda = NULL`.
- Added a check for a non-positive or non-finite `lambda.max` to prevent
  the generation of invalid lambda sequences.
- Fixed within-group weight generation in
  [`gen_prec_sbm()`](https://shiying-xiao.com/grasps/reference/gen_prec_sbm.md)
  so that each edge weight is sampled once and assigned symmetrically to
  the corresponding entries of the precision matrix.
- Added input validation for `cond.target` and support for empty graphs
  in
  [`gen_prec_sbm()`](https://shiying-xiao.com/grasps/reference/gen_prec_sbm.md).
- Clarified the roles of `lambda.safe` and `lambda.max` in the
  documentation of the automatic lambda-grid construction.
- Expanded the documentation of the positive-definiteness and
  condition-number adjustment in
  [`gen_prec_sbm()`](https://shiying-xiao.com/grasps/reference/gen_prec_sbm.md)
  and added relevant references.

## grasps 0.1.1

CRAN release: 2026-05-02

- Added function
  [`plot.adjmat()`](https://shiying-xiao.com/grasps/reference/plot.adjmat.md).
- Added function
  [`prec_to_adj()`](https://shiying-xiao.com/grasps/reference/prec_to_adj.md).
- Use `d` instead of `p` to denote the dimension.
- Revised the estimator expression in the vignette.
- Added NEWS file to record the changelog.

## grasps 0.1.0

CRAN release: 2025-11-27

- Initial release of the grasps package.
- The grasps is a toolbox for precision matrix estimation with group
  structure.
- Key features:
  - Unified regularization framework for sparse network learning.
  - Supports element-wise sparsity and group-wise shrinkage.
  - Includes both convex and non-convex penalties (e.g., adaptive lasso,
    SCAD, MCP).
  - Provides model selection tools (e.g., CV, EBIC, HBIC).
- Designed for applications where variables exhibit grouped or modular
  relationships (e.g., brain networks, biological pathways).
