#' Performance Measures for Precision Matrix Estimation
#'
#' @description
#' Compute a collection of loss-based and structure-based measures to evaluate
#' the performance of an estimated precision matrix.
#'
#' @param hatOmega A numeric symmetric \eqn{p \times p} matrix, with \eqn{p} > 1,
#' specifying the estimated precision matrix.
#'
#' @param Omega A numeric symmetric \eqn{p \times p} matrix, with \eqn{p} > 1,
#' specifying the reference (typically true) precision matrix.
#' It must have the same dimensions as \code{hatOmega}.
#'
#' @return
#' A data frame of S3 class \code{"performance"}, with one row per performance
#' metric and two columns:
#' \describe{
#' \item{measure}{The name of each performance metric. The reported metrics
#' include: sparsity, Frobenius norm loss, Kullback-Leibler divergence,
#' quadratic norm loss, spectral norm loss, true positive, true negative,
#' false positive, false negative, true positive rate, false positive rate,
#' F1 score, and Matthews correlation coefficient.}
#' \item{value}{The corresponding numeric value.
#' Undefined measures are returned as \code{NA}.}
#' }
#'
#' @details
#' Let \eqn{\Omega_{p \times p}} and \eqn{\hat{\Omega}_{p \times p}} be
#' the reference (true) and estimated precision matrices, respectively.
#' When \eqn{\Omega} is positive definite, \eqn{\Sigma = \Omega^{-1}} is
#' the corresponding covariance matrix.
#' Edges are defined by nonzero off-diagonal entries in the upper triangle of
#' the precision matrices.
#'
#' The input matrices are not required to be positive definite.
#' If either matrix is not positive definite, the function issues a warning
#' but continues computing the performance measures.
#' In this case, covariance- and likelihood-based measures may not have their
#' usual statistical interpretation.
#' If \eqn{\Omega} is singular, the KL and quadratic losses are returned as
#' \code{NA}, while the remaining measures are still computed.
#' For non-positive-definite but invertible matrices, the reported KL value uses
#' the log absolute determinant and need not be non-negative.
#' A classification-based measure is returned as \code{NA} when its denominator
#' is zero.
#'
#' Sparsity is treated as a structural summary, while the remaining measures
#' are grouped into loss-based measures, raw confusion-matrix counts, and
#' classification-based (structure-recovery) measures.
#'
#' "sparsity": Sparsity is computed as the proportion of zero entries among
#' the off-diagonal elements in the upper triangle of \eqn{\hat{\Omega}}.
#'
#' \strong{Loss-based measures:}
#' \itemize{
#' \item "Frobenius": Frobenius (Hilbert-Schmidt) norm loss
#' \eqn{= \Vert \Omega - \hat{\Omega} \Vert_F}.
#' \item "KL": Kullback-Leibler divergence
#' \eqn{= \mathrm{tr}(\Sigma \hat{\Omega}) - \log\det(\Sigma \hat{\Omega}) - p}.
#' \item "quadratic": Quadratic norm loss
#' \eqn{= \Vert \Sigma \hat{\Omega} - I_p \Vert_F^2}.
#' \item "spectral": Spectral (operator) norm loss
#' \eqn{= \Vert \Omega - \hat{\Omega} \Vert_{2,2} = e_1},
#' where \eqn{e_1^2} is the largest eigenvalue of \eqn{(\Omega - \hat{\Omega})^2}.
#' }
#'
#' \strong{Confusion-matrix counts:}
#' \itemize{
#' \item "TP": True positive \eqn{=} number of edges
#' in both \eqn{\Omega} and \eqn{\hat{\Omega}}.
#' \item "TN": True negative \eqn{=} number of edges
#' in neither \eqn{\Omega} nor \eqn{\hat{\Omega}}.
#' \item "FP": False positive \eqn{=} number of edges
#' in \eqn{\hat{\Omega}} but not in \eqn{\Omega}.
#' \item "FN": False negative \eqn{=} number of edges
#' in \eqn{\Omega} but not in \eqn{\hat{\Omega}}.
#' }
#'
#' \strong{Classification-based (structure-recovery) measures:}
#' \itemize{
#' \item "TPR": True positive rate (TPR), recall, sensitivity
#' \eqn{= \mathrm{TP} / (\mathrm{TP} + \mathrm{FN})}.
#' \item "FPR": False positive rate (FPR)
#' \eqn{= \mathrm{FP} / (\mathrm{FP} + \mathrm{TN})}.
#' \item "F1": \eqn{F_1} score
#' \eqn{= 2\,\mathrm{TP} / (2\,\mathrm{TP} + \mathrm{FN} + \mathrm{FP})}
#' \item "MCC": Matthews correlation coefficient (MCC)
#' \eqn{= (\mathrm{TP}\times\mathrm{TN} - \mathrm{FP}\times\mathrm{FN}) /
#' \sqrt{(\mathrm{TP}+\mathrm{FP})(\mathrm{TP}+\mathrm{FN})
#' (\mathrm{TN}+\mathrm{FP})(\mathrm{TN}+\mathrm{FN})}}
#' }
#'
#' The following table summarizes the confusion matrix and related rates:
#' \tabular{lllll}{
#' \tab \strong{Predicted Positive} \tab \strong{Predicted Negative} \tab \tab \cr
#' \strong{Real Positive} (P) \tab True positive (TP) \tab False negative (FN)
#' \tab True positive rate (TPR), recall, sensitivity = TP / P = 1 - FNR
#' \tab False negative rate (FNR) = FN / P = 1 - TPR \cr
#' \strong{Real Negative} (N) \tab False positive (FP) \tab True negative (TN)
#' \tab False positive rate (FPR) = FP / N = 1 - TNR
#' \tab True negative rate (TNR), specificity = TN / N = 1 - FPR \cr
#' \tab Positive predictive value (PPV), precision = TP / (TP + FP) = 1 - FDR
#' \tab False omission rate (FOR) = FN / (TN + FN) = 1 - NPV \tab \tab \cr
#' \tab False discovery rate (FDR) = FP / (TP + FP) = 1 - PPV
#' \tab Negative predictive value (NPV) = TN / (TN + FN) = 1 - FOR \tab \tab \cr
#' }
#'
#' @example
#' inst/example/ex-grasps.R
#'
#' @export

performance <- function(hatOmega, Omega) {

  ## validate `hatOmega`, `Omega`
  validate_symmat(hatOmega, "hatOmega", min.dim = 2)
  validate_symmat(Omega, "Omega", min.dim = 2)
  if (!identical(dim(hatOmega), dim(Omega))) {
    stop("`hatOmega` and `Omega` must have the same dimensions.", call. = FALSE)
  }

  ## check positive definiteness
  chol.hatOmega <- tryCatch(chol(hatOmega), error = function(e) NULL)
  if (is.null(chol.hatOmega)) {
    warning(paste(
      "`hatOmega` is not positive definite.",
      "Likelihood-based measures may not have their usual interpretation."
    ), call. = FALSE)
  }
  chol.Omega <- tryCatch(chol(Omega), error = function(e) NULL)

  ## true covariance
  Sigma <- tryCatch(solve(Omega), error = function(e) NULL)
  if (is.null(Sigma)) {
    warning(paste("`Omega` is singular.",
                  "KL and quadratic losses are returned as `NA`."),
            call. = FALSE)
  } else if (is.null(chol.Omega)) {
    warning(paste(
      "`Omega` is not positive definite.",
      "Covariance- and likelihood-based measures may not have their",
      "usual interpretation."
    ), call. = FALSE)
  }

  ## dimension
  p <- ncol(hatOmega)

  ## Frobenius norm loss
  FL <- norm(Omega - hatOmega, "F")

  ## Kullback-Leibler divergence
  ## quadratic loss
  if (is.null(Sigma)) {
    KL <- NA_real_
    QL <- NA_real_
  } else {
    Sigma.hatOmega <- Sigma %*% hatOmega
    KL <- sum(diag(Sigma.hatOmega)) -
      determinant(Sigma.hatOmega, logarithm = TRUE)$modulus[1] - p
    QL <- norm(Sigma.hatOmega - diag(p), "F")^2
  }

  ## spectral norm loss
  SL <- svd(Omega - hatOmega)$d[1]
  ## Equivalent version:
  ## temp <- (Omega - hatOmega) %*% (Omega - hatOmega)
  ## e1 <- max(Re(eigen(temp)$values))
  ## SL <- sqrt(e1)

  ## off-diagonal upper triangle entries define edges
  Omega_edge <- Omega[upper.tri(Omega)]
  hatOmega_edge <- hatOmega[upper.tri(hatOmega)]

  ## sparsity
  sparsity <- sum(hatOmega_edge == 0) / length(hatOmega_edge)

  ## confusion-matrix counts
  TP <- sum(Omega_edge != 0 & hatOmega_edge != 0)
  TN <- sum(Omega_edge == 0 & hatOmega_edge == 0)
  FP <- sum(Omega_edge == 0 & hatOmega_edge != 0)
  FN <- sum(Omega_edge != 0 & hatOmega_edge == 0)

  ## classification-based (structure-recovery) measures
  TPR <- if (TP + FN > 0) TP / (TP + FN) else NA_real_
  FPR <- if (TN + FP > 0) FP / (TN + FP) else NA_real_
  F1 <- if (2 * TP + FN + FP > 0) 2 * TP / (2 * TP + FN + FP) else NA_real_
  MCC <- if (all(c(TP + FP, TP + FN, TN + FP, TN + FN) > 0)) {
    (TP * TN - FP * FN) /
      (sqrt(TP + FP) * sqrt(TP + FN) * sqrt(TN + FP) * sqrt(TN + FN))
  } else {
    NA_real_
  }

  result <- data.frame(
    measure = c("sparsity", "Frobenius", "KL", "quadratic", "spectral",
                "TP", "TN", "FP", "FN", "TPR", "FPR", "F1", "MCC"),
    value = c(sparsity, FL, KL, QL, SL,
              TP, TN, FP, FN, TPR, FPR, F1, MCC)
  )
  class(result) <- c("performance", class(result))
  return(result)
}


#' @export
#' @noRd

print.performance <- function(x, ...) {
  x$value <- round(x$value, 4)
  NextMethod("print", x)
}
