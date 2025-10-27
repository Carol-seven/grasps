#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;


static const double eps = 1e-10;


// Element-Wise Tuning Parameter Matrix of Individual Penalty
// "lasso": lambda
// "adapt": lambda / \vert Omega_{ij} \vert^gamma
// "mcp": lambda - \vert Omega_{ij} \vert / gamma  if \vert Omega_{ij} \vert \leq gamma*lambda
// "scad": lambda                                                 if \vert Omega_{ij} \vert \leq lambda
//         (gamma*lambda - \vert Omega_{ij} \vert) / (gamma - 1)  if lambda < \vert Omega_{ij} \vert < gamma*lambda
arma::mat lambda_individual(const arma::mat& Omega, std::string penalty,
                            double lambda, double gamma) {

  arma::mat W = arma::zeros(Omega.n_rows, Omega.n_cols);

  if (penalty == "lasso") {
    W.fill(lambda);

  } else if (penalty == "adapt") {
    W = lambda / (arma::pow(arma::abs(Omega), gamma) + eps);

  } else if (penalty == "mcp") {
    arma::mat abs_Omega = arma::abs(Omega);
    W = (lambda - abs_Omega / gamma) % (abs_Omega <= gamma * lambda);

  } else if (penalty == "scad") {
    arma::mat abs_Omega = arma::abs(Omega);
    W.fill(lambda);
    W = W % (abs_Omega <= lambda);
    W += ((gamma * lambda - abs_Omega) / (gamma - 1)) %
      ((abs_Omega > lambda) % (abs_Omega <= gamma * lambda));

  } else {
    stop("Unsupported penalty type!");
  }

  return W;
}


// Tuning Parameter Scalar of Group Penalty
// "lasso": lambda
// "adapt": lambda / \Vert V \Vert_2, where V_{ij} = \vert Omega_{ij} \vert^gamma
// "mcp": lambda - \Vert Omega \Vert_2 / gamma  if \Vert Omega \Vert_2 \leq gamma*lambda
// "scad": lambda                                              if \Vert Omega \Vert_2 \leq lambda
//         (gamma*lambda - \Vert Omega \Vert_2) / (gamma - 1)  if lambda < \Vert Omega \Vert_2 < gamma*lambda
double lambda_group(const arma::mat& Omega, std::string penalty,
                    double lambda, double gamma) {

  if (penalty == "lasso") {
    return lambda;

  } else if (penalty == "adapt") {
    return lambda / (arma::norm(arma::pow(arma::abs(Omega), gamma), "fro") + eps);

  } else if (penalty == "mcp") {
    double Omega_norm = arma::norm(Omega, "fro");
    if (Omega_norm <= gamma * lambda) {
      return lambda - Omega_norm / gamma;
    } else {
      return 0.0;
    }

  } else if (penalty == "scad") {
    double Omega_norm = arma::norm(Omega, "fro");
    if (Omega_norm <= lambda) {
      return lambda;
    } else if (Omega_norm <= gamma * lambda) {
      return (gamma * lambda - Omega_norm) / (gamma - 1);
    } else {
      return 0.0;
    }

  } else {
    stop("Unsupported penalty type!");
  }
}

