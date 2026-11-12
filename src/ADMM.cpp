#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <algorithm> // for std::max
#include <cmath> // for std::isfinite()
#include "internals.h"
using namespace Rcpp;


// Sparse-Group Graphical Lasso via ADMM
List ADMMsggl(
    const arma::mat& S, const List& group_idx,
    bool diag_ind, bool diag_grp, bool diag_include,
    const double& lambda, const double& alpha,
    double rho, const double tau_incr, const double tau_decr, const double nu,
    const double tol_abs, const double tol_rel, const int maxiter) {

  int p = S.n_cols;
  int iterations = 0;
  bool converged = false;
  arma::mat Omega = arma::diagmat(1 / arma::diagvec(S));
  arma::mat Z = arma::zeros(p, p);
  arma::mat U = arma::zeros(p, p);
  arma::mat Z_old = Z;

  arma::mat I = arma::ones(p, p);
  if (!diag_ind) {
    I -= arma::diagmat(I);
  }

  int grp_length = group_idx.length();

  // loop until convergence
  for (int iter = 0; iter < maxiter; ++iter) {

    Z_old = Z;

    // update Omega
    Omega = update_Omega(rho * (Z - U) - S, rho);

    // update Z
    // step 1: element-wise individual L1 penalty
    Z = soft_matrix(Omega + U, lambda * I * alpha / rho);
    arma::vec Z_diag = arma::diagvec(Z);
    // step 2: block-wise group L2 penalty
    // update Z_{gh}
    for (int g = 0; g < grp_length; ++g) {
      // rows corresponding to group g
      arma::uvec row_idx = group_idx[g];
      for (int h = g; h < grp_length; ++h) {
        if (g == h && !diag_grp) {
          continue;
        }
        // columns corresponding to group h
        arma::uvec col_idx = group_idx[h];
        // get subview for in-place update
        auto Z_block = Z.submat(row_idx, col_idx);
        if (g != h) {
          double block_norm = arma::norm(Z_block, "fro");
          double scale_grp = std::max(
            1.0 - lambda * (1 - alpha) / (rho * (block_norm + eps)),
            0.0
          );
          Z_block *= scale_grp;
          Z.submat(col_idx, row_idx) = Z_block.t();
          continue;
        }
        if (diag_include) {
          double block_norm = arma::norm(Z_block, "fro");
          double scale_grp = std::max(
            1.0 - lambda * (1 - alpha) / (rho * (block_norm + eps)),
            0.0
          );
          Z_block *= scale_grp;
        } else {
          double sq_all = arma::accu(arma::square(Z_block));
          arma::vec de = Z_diag.elem(row_idx);
          double sq_diag = arma::dot(de, de);
          double offblock_norm = std::sqrt(std::max(sq_all - sq_diag, 0.0));
          double scale_grp = std::max(
            1.0 - lambda * (1 - alpha) / (rho * (offblock_norm + eps)),
            0.0
          );
          Z_block *= scale_grp;
          for (arma::uword k = 0; k < row_idx.n_elem; ++k) {
            Z(row_idx[k], row_idx[k]) = de[k];
          }
        }
      }
    }

    // update U
    U += Omega - Z;

    // residual
    double r = arma::norm(Omega - Z, "fro");
    double s = arma::norm(rho * (Z - Z_old), "fro");

    // tolerance
    double tol_pri = p * tol_abs +
      tol_rel * std::max(arma::norm(Omega, "fro"), arma::norm(Z, "fro"));
    double tol_dual = p * tol_abs +
      tol_rel * rho * arma::norm(U, "fro");

    iterations = iter + 1;

    if (r <= tol_pri && s <= tol_dual) {
      converged = true;
      break;
    }

    // update rho
    if (r > nu * s) {
      rho *= tau_incr;
      U /= tau_incr;
    }
    if (s > nu * r) {
      rho /= tau_decr;
      U *= tau_decr;
    }

    // R_CheckUserInterrupt
    if (iter % 1000 == 0) {
      R_CheckUserInterrupt();
    }
  }

  return List::create(
    Named("hatOmega") = Z,
    Named("lambda") = lambda,
    Named("alpha") = alpha,
    Named("iterations") = iterations,
    Named("converged") = converged
  );
}


// Fixed-Weight LLA Subproblem via ADMM
List ADMMlla(
    const arma::mat& S, const List& group_idx, const std::string& penalty,
    bool diag_ind, bool diag_grp, bool diag_include,
    const double& lambda, const double& alpha, const double& gamma,
    const arma::mat& weight_ref,
    double rho, const double tau_incr, const double tau_decr, const double nu,
    const double tol_abs, const double tol_rel, const int maxiter) {

  int p = S.n_cols;
  int iterations = 0;
  bool converged = false;

  arma::mat Omega = arma::diagmat(1 / arma::diagvec(S));
  arma::mat Z = arma::zeros(p, p);
  arma::mat U = arma::zeros(p, p);
  arma::mat Z_old = Z;

  arma::mat I = arma::ones(p, p);
  if (!diag_ind) {
    I -= arma::diagmat(I);
  }

  // penalty weights remain fixed within this ADMM subproblem
  arma::mat lambda_ind =
    lambda_individual(weight_ref, penalty, lambda, gamma) % I;

  // loop until convergence
  for (int iter = 0; iter < maxiter; ++iter) {

    Z_old = Z;

    // update Omega
    Omega = update_Omega(rho * (Z - U) - S, rho);

    // update Z
    // step 1: element-wise individual L1 penalty
    Z = soft_matrix(Omega + U, lambda_ind * alpha / rho);

    // step 2: block-wise group L2 penalty
    Z = update_Zblock(Z, penalty, group_idx,
                      lambda, gamma, weight_ref, alpha,
                      diag_grp, diag_include, rho);

    // update U
    U += Omega - Z;

    // residuals
    double r = arma::norm(Omega - Z, "fro");
    double s = arma::norm(rho * (Z - Z_old), "fro");

    // tolerances
    double tol_pri = p * tol_abs +
      tol_rel * std::max(arma::norm(Omega, "fro"), arma::norm(Z, "fro"));
    double tol_dual = p * tol_abs +
      tol_rel * rho * arma::norm(U, "fro");

    iterations = iter + 1;

    if (r <= tol_pri && s <= tol_dual) {
      converged = true;
      break;
    }

    // update rho
    if (r > nu * s) {
      rho *= tau_incr;
      U /= tau_incr;
    }

    if (s > nu * r) {
      rho /= tau_decr;
      U *= tau_decr;
    }

    // R_CheckUserInterrupt
    if (iter % 1000 == 0) {
      R_CheckUserInterrupt();
    }
  }

  return List::create(
    Named("hatOmega") = Z,
    Named("iterations") = iterations,
    Named("converged") = converged
  );
}


// Sparse-Group Graphical Non-Convex Penalties via LLA
List ADMMsggn(
    const arma::mat& S, const List& group_idx, const std::string& penalty,
    bool diag_ind, bool diag_grp, bool diag_include,
    const double& lambda, const double& alpha, const double& gamma,
    const Rcpp::RObject& initial,
    bool lla, const double tol_lla, const int maxiter_lla,
    double rho, const double tau_incr, const double tau_decr, const double nu,
    const double tol_abs, const double tol_rel, const int maxiter) {

  // initial estimate used to construct penalty weights
  arma::mat initial_est;
  Rcpp::LogicalVector initial_converged(1, NA_LOGICAL);

  if (Rf_isString(initial)) {
    List sgglres = ADMMsggl(S, group_idx,
                            diag_ind, diag_grp, diag_include,
                            lambda, alpha,
                            rho, tau_incr, tau_decr, nu,
                            tol_abs, tol_rel, maxiter);
    initial_est = Rcpp::as<arma::mat>(sgglres["hatOmega"]);
    initial_converged[0] = Rcpp::as<bool>(sgglres["converged"]);
  } else {
    initial_est = Rcpp::as<arma::mat>(initial);
  }

  // adaptive lasso weights remain fixed at the initial estimate
  bool iterative = lla && penalty != "adapt";
  int steps = iterative ? maxiter_lla : 1;

  arma::mat weight_ref = initial_est;
  arma::mat hatOmega = initial_est;

  int iterations = 0;
  int lla_iterations = 0;
  bool all_admm_converged = true;
  bool lla_converged = false;

  for (int lla_iter = 0; lla_iter < steps; ++lla_iter) {

    List step = ADMMlla(S, group_idx, penalty,
                        diag_ind, diag_grp, diag_include,
                        lambda, alpha, gamma,
                        weight_ref,
                        rho, tau_incr, tau_decr, nu,
                        tol_abs, tol_rel, maxiter);

    hatOmega = Rcpp::as<arma::mat>(step["hatOmega"]);

    bool step_converged = Rcpp::as<bool>(step["converged"]);

    iterations += Rcpp::as<int>(step["iterations"]);
    ++lla_iterations;

    if (!step_converged) {
      all_admm_converged = false;
      break;
    }

    if (iterative) {
      double relative_change =
        arma::norm(hatOmega - weight_ref, "fro") /
          std::max(1.0, arma::norm(weight_ref, "fro"));

      if (relative_change <= tol_lla) {
        lla_converged = true;
        break;
      }

      weight_ref = hatOmega;
    }
  }

  return List::create(
    Named("hatOmega") = hatOmega,
    Named("lambda") = lambda,
    Named("alpha") = alpha,
    Named("initial") = initial_est,
    Named("gamma") = gamma,
    Named("iterations") = iterations,
    Named("converged") = initial_converged[0] != FALSE && all_admm_converged,
    Named("initial.converged") = initial_converged,
    Named("lla.iterations") = penalty == "adapt" ? NA_INTEGER : lla_iterations,
    Named("lla.converged") =
      Rcpp::LogicalVector::create(iterative ? lla_converged : NA_LOGICAL)
  );
}


// Sparse-Group Graphical Model via ADMM
// [[Rcpp::export]]
List ADMMsggm(
    const arma::mat& S, const List& group_idx, const std::string& penalty,
    bool diag_ind, bool diag_grp, bool diag_include,
    const double& lambda, const double& alpha, const double& gamma,
    const Rcpp::RObject& initial,
    bool lla, const double tol_lla, const int maxiter_lla,
    double rho, const double tau_incr, const double tau_decr, const double nu,
    const double tol_abs, const double tol_rel, const int maxiter) {

  if (penalty == "lasso") {
    List result = ADMMsggl(S, group_idx,
                           diag_ind, diag_grp, diag_include,
                           lambda, alpha,
                           rho, tau_incr, tau_decr, nu,
                           tol_abs, tol_rel, maxiter);
    result["lla.iterations"] = NA_INTEGER;
    result["lla.converged"] = Rcpp::LogicalVector::create(NA_LOGICAL);
    return result;
  } else {
    Rcpp::RObject initial_est = gen_initial(initial, S);
    return ADMMsggn(S, group_idx, penalty,
                    diag_ind, diag_grp, diag_include,
                    lambda, alpha, gamma, initial_est,
                    lla, tol_lla, maxiter_lla,
                    rho, tau_incr, tau_decr, nu,
                    tol_abs, tol_rel, maxiter);
  }
}


// ADMMsggm_CV
// [[Rcpp::export]]
List ADMMsggm_CV(
    const arma::mat& X, const List& group_idx, const std::string& penalty,
    bool diag_ind, bool diag_grp, bool diag_include,
    const arma::colvec& lambdas, const arma::colvec& alphas, const double& gamma,
    const Rcpp::RObject& initial,
    bool lla, const double tol_lla, const int maxiter_lla,
    double rho, const double tau_incr, const double tau_decr, const double nu,
    const double tol_abs, const double tol_rel, const int maxiter,
    int kfold) {

  int n = X.n_rows;
  int nlambda = lambdas.n_elem;

  arma::mat CV_loss(nlambda, kfold, arma::fill::zeros);
  arma::umat converged_grid(nlambda, kfold, arma::fill::zeros);
  Rcpp::LogicalMatrix lla_converged_grid(nlambda, kfold);
  std::fill(lla_converged_grid.begin(), lla_converged_grid.end(), NA_LOGICAL);
  arma::uvec index = index_kfold(n, kfold);

  for (int k = 0; k < kfold; ++k) {

    // indices for test and training sets; test and training sets
    arma::uvec idx_test = arma::find(index == (unsigned)k);
    arma::mat X_test  = X.rows(idx_test);
    int n_test = X_test.n_rows;
    arma::uvec idx_train = arma::find(index != (unsigned)k);
    arma::mat X_train = X.rows(idx_train);

    // empirical covariance matrix
    arma::mat S_test = arma::cov(X_test, 1);
    arma::mat S_train = arma::cov(X_train, 1);

    Rcpp::RObject initial_cv = initial;
    if (penalty != "lasso") {
      initial_cv = gen_initial(initial, S_train);
    }

    for (int i = 0; i < nlambda; ++i) {

      double lambda = lambdas[i];
      double alpha = alphas[i];

      List result = ADMMsggm(S_train, group_idx, penalty,
                             diag_ind, diag_grp, diag_include,
                             lambda, alpha, gamma, initial_cv,
                             lla, tol_lla, maxiter_lla,
                             rho, tau_incr, tau_decr, nu,
                             tol_abs, tol_rel, maxiter);

      bool converged = Rcpp::as<bool>(result["converged"]);
      converged_grid(i, k) = converged;

      Rcpp::LogicalVector lla_status = result["lla.converged"];
      lla_converged_grid(i, k) = lla_status[0];

      // negative Gaussian log-likelihood
      arma::mat hatOmega = result["hatOmega"];
      CV_loss(i,k) = criterion(hatOmega, S_test, n_test, "negloglik", 0.0);

    }
  }

  // the mean of the k-fold loss for each parameter grid value
  arma::mat loss_avg = arma::mean(CV_loss, 1);
  if (!std::isfinite(loss_avg.min())) {
    stop("No candidate has a finite cross-validation score. "
         "Try increasing `maxiter` or changing the tuning grid.");
  }
  arma::uword opt_idx = loss_avg.index_min();
  double loss_min = loss_avg(opt_idx);
  double lambda_opt = lambdas[opt_idx];
  double alpha_opt = alphas[opt_idx];

  return List::create(
    Named("lambda.opt") = lambda_opt,
    Named("alpha.opt") = alpha_opt,
    Named("CV.loss") = CV_loss,
    Named("loss.min") = loss_min,
    Named("converged.grid") = converged_grid,
    Named("lla.converged.grid") = lla_converged_grid
  );
}


// ADMMsggm_IC
// [[Rcpp::export]]
List ADMMsggm_IC(
    const arma::mat& S, const List& group_idx, const std::string& penalty,
    bool diag_ind, bool diag_grp, bool diag_include,
    const arma::colvec& lambdas, const arma::colvec& alphas, const double& gamma,
    const Rcpp::RObject& initial,
    bool lla, const double tol_lla, const int maxiter_lla,
    double rho, const double tau_incr, const double tau_decr, const double nu,
    const double tol_abs, const double tol_rel, const int maxiter,
    const std::string& crit, int n, double ebic_tuning) {

  int nlambda = lambdas.n_elem;
  arma::vec IC_score(nlambda, arma::fill::zeros);
  arma::uvec converged_grid(nlambda, arma::fill::zeros);
  Rcpp::LogicalVector lla_converged_grid(nlambda, NA_LOGICAL);
  List result_list(nlambda);

  Rcpp::RObject initial_ic = initial;
  if (penalty != "lasso") {
    initial_ic = gen_initial(initial, S);
  }

  for (int i = 0; i < nlambda; ++i) {

    double lambda = lambdas[i];
    double alpha = alphas[i];
    List result = ADMMsggm(S, group_idx, penalty,
                           diag_ind, diag_grp, diag_include,
                           lambda, alpha, gamma, initial_ic,
                           lla, tol_lla, maxiter_lla,
                           rho, tau_incr, tau_decr, nu,
                           tol_abs, tol_rel, maxiter);
    result_list[i] = result;

    bool converged = Rcpp::as<bool>(result["converged"]);
    converged_grid[i] = converged;

    Rcpp::LogicalVector lla_status = result["lla.converged"];
    lla_converged_grid[i] = lla_status[0];

    arma::mat hatOmega = result["hatOmega"];
    IC_score[i] = criterion(hatOmega, S, n, crit, ebic_tuning);

  }

  // the mean of the score for each parameter grid value
  if (!std::isfinite(IC_score.min())) {
    stop("No candidate has a finite information-criterion score. "
         "Try increasing `maxiter` or changing the tuning grid.");
  }
  arma::uword opt_idx = IC_score.index_min();
  double score_min = IC_score(opt_idx);
  double lambda_opt = lambdas[opt_idx];
  double alpha_opt = alphas[opt_idx];

  return List::create(
    Named("lambda.opt") = lambda_opt,
    Named("alpha.opt") = alpha_opt,
    Named("IC.score") = IC_score,
    Named("score.min") = score_min,
    Named("converged.grid") = converged_grid,
    Named("lla.converged.grid") = lla_converged_grid,
    Named("result") = result_list[opt_idx]
  );
}

