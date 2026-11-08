#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


// Element-Wise Soft-Thresholding on Matrix
// sign(M_{ij}) * max(\vert M_{ij} \vert - Thres_{ij}, 0)
// [[Rcpp::export]]
arma::mat soft_matrix(const arma::mat& M, const arma::mat& Thres) {
  return arma::sign(M) % arma::max(arma::abs(M) - Thres, arma::zeros<arma::mat>(M.n_rows, M.n_cols));
}


// Index of K-Fold
arma::uvec index_kfold(const int n, const int kfold) {
  arma::uvec folds(n);
  for (int i = 0; i < n; ++i) {
    folds[i] = i % kfold;
  }
  return arma::shuffle(folds);
}


// criterion
double criterion(const arma::mat& hatOmega, const arma::mat& S, int n,
                 const std::string& crit, double ebic_tuning) {

  // dimensionality
  int p = S.n_cols;
  // Gaussian log-likelihood
  arma::mat factor;
  if (!hatOmega.is_finite() || !arma::chol(factor, hatOmega)) {
    return std::numeric_limits<double>::infinity();
  }
  double logdet = 2.0 * arma::accu(arma::log(factor.diag()));
  double loglik = (n/2.0) * (logdet - arma::accu(hatOmega % S));
  if (!std::isfinite(loglik)) {
    return std::numeric_limits<double>::infinity();
  }

  if (crit == "negloglik") {
    return -loglik;
  }

  // cardinality of the edge set
  arma::uword edges = arma::accu(arma::trimatu(hatOmega, 1) != 0);

  double score;
  if (crit == "AIC") {
    score = -2*loglik + 2*edges;
  }
  else if (crit == "BIC") {
    score = -2*loglik + std::log(n)*edges;
  }
  else if (crit == "EBIC") {
    score = -2*loglik + std::log(n)*edges + 4*ebic_tuning*std::log(p)*edges;
  }
  else if (crit == "HBIC") {
    score = -2*loglik + std::log(std::log(n))*std::log(p)*edges;
  }
  else {
    Rcpp::stop("Unknown criterion!");
  }
  return std::isfinite(score) ? score : std::numeric_limits<double>::infinity();
}

