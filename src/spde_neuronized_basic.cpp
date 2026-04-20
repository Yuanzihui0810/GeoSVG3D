// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// ============================================================
// 1. NPrior-style activation function
//    directly adapted from NPrior/src/Neuro_Linear.cpp
//
// type = 1 : BL      -> T(t) = t
// type = 2 : HS      -> T(t) = exp(0.5*sign(t)*t^2 + 0.733*t)
// type = 3 : SpSL-L  -> T(t) = max(0, t)
// type = 4 : SpSL-C  -> T(t) = exp(0.5*t^2 - 1.27*t + 0.29) * 1(t>0)
// type = 5 : custom continuous (abs-quadratic)
//            T(t) = exp(lam1 * t * |t| + lam2 * t + lam3)
// type = 6 : custom continuous (quadratic)
//            T(t) = exp(lam1 * t^2 + lam2 * t + lam3)
// type = 7 : custom thresholded (abs-quadratic)
//            T(t) = exp(lam1 * t * |t| + lam2 * t + lam3) * I(t > 0)
// type = 8 : custom thresholded (quadratic)
//            T(t) = exp(lam1 * t^2 + lam2 * t + lam3) * I(t > 0)
// ============================================================
inline double T_fun_np_scalar(const double t, const int type,
                              const double lam1 = 0.0,
                              const double lam2 = 0.0,
                              const double lam3 = 0.0) {
  if (type == 1) {
    return t;
  } else if (type == 2) {
    double z = 0.5 * t * std::abs(t) + 0.733 * t;
    return std::exp(z);
  } else if (type == 3) {
    return (t > 0.0) ? t : 0.0;
  } else if (type == 4) {
    if (t > 0.0) {
      double z = 0.5 * t * t - 1.27 * t + 0.29;
      return std::exp(z);
    } else {
      return 0.0;
    }
  } else if (type == 5) {   // custom continuous: exp(lam1 * t|t| + lam2 * t + lam3)
    double z = lam1 * t * std::abs(t) + lam2 * t + lam3;
    return std::exp(z);
  }else if (type == 6) {   // custom continuous: exp(lam1 * t^2 + lam2 * t + lam3)
    double z = lam1 * t * t + lam2 * t + lam3;
    return std::exp(z);
  }else if (type == 7) {   // custom thresholded: exp(lam1 * t|t| + lam2 * t + lam3) * I(t>0)
    if (t <= 0.0) return 0.0;
    double z = lam1 * t * std::abs(t) + lam2 * t + lam3;
    return std::exp(z);
  }else if (type == 8) {   // custom thresholded: exp(lam1 * t^2 + lam2 * t + lam3) * I(t>0)
    if (t <= 0.0) return 0.0;
    double z = lam1 * t * t + lam2 * t + lam3;
    return std::exp(z);
  }
  Rcpp::stop("unsupported type in T_fun_np_scalar");
}

// vectorized activation via Armadillo transform
arma::vec T_fun_np_vec(const arma::vec& t, const int type,
                       const double lam1 = 0.0,
                       const double lam2 = 0.0,
                       const double lam3 = 0.0) {
  arma::vec out = t;
  if (type == 1) {
    return out;
  } else if (type == 2) {
    out = arma::exp(0.5 * t % arma::abs(t) + 0.733 * t);
    return out;
  } else if (type == 3) {
    out = arma::clamp(t, 0.0, arma::datum::inf);
    return out;
  } else if (type == 4) {
    out.zeros();
    arma::uvec ind = arma::find(t > 0.0);
    if (!ind.is_empty()) {
      arma::vec tp = t.elem(ind);
      out.elem(ind) = arma::exp(0.5 * tp % tp - 1.27 * tp + 0.29);
    }
    return out;
  } else if (type == 5) {   // custom continuous abs-quadratic
    return arma::exp(lam1 * (t % arma::abs(t)) + lam2 * t + lam3);
  }else if (type == 6) {   // custom continuous quadratic
    return arma::exp(lam1 * (t % t) + lam2 * t + lam3);
  }else if (type == 7) {   // custom thresholded abs-quadratic
    out.zeros();
    arma::uvec ind = arma::find(t > 0.0);
    if (!ind.is_empty()) {
      arma::vec tp = t.elem(ind);
      out.elem(ind) = arma::exp(lam1 * (tp % arma::abs(tp)) + lam2 * tp + lam3);
    }
    return out;
  }else if (type == 8) {   // custom thresholded quadratic
    out.zeros();
    arma::uvec ind = arma::find(t > 0.0);
    if (!ind.is_empty()) {
      arma::vec tp = t.elem(ind);
      out.elem(ind) = arma::exp(lam1 * (tp % tp) + lam2 * tp + lam3);
    }
    return out;
  }
  Rcpp::stop("unsupported type in T_fun_np_vec");
}

// ============================================================
// 2. one-gene MCMC under current workflow
//
// y_g = mu_g 1_n + B theta_g + eps_g
// theta_gk = s_gk * T(xi_gk - xi_g0) * omega_gk
// s_gk = (kappa_g^2 + lambda_k)^(-alpha_g/2)
//
// MCMC steps:
//   1) update mu_g
//   2) update sigma_g^2
//   3) for k=1..K:
//        marginalized MH for xi_gk
//        Gibbs for omega_gk
//   4) update tau_g^2
//
// Only the activation/prior family is borrowed from NPrior.
// ============================================================

// [[Rcpp::export]]
Rcpp::List spde_neuronized_basic_one_gene_cpp(
    const arma::vec& y_g,      // n x 1
    const arma::mat& B,        // n x K
    const arma::vec& lambda,   // K x 1
    const int N,
    const int BURN,
    const int thin,
    const int mh_per_k,
    const double kappa_g,
    const double alpha_g,
    const double xi_g0,
    const double nu_g0,
    const double eta_g0_sq,
    const double m_g0,
    const double gamma_g0,
    const double a_g0,
    const double b_g0,
    const int prior_type,
    const double lam1,
    const double lam2,
    const double lam3,
    const double xi_prop_sd,
    const int tau2_g_update,
    double mu_g,
    double sigma2_g,
    double tau2_g,
    arma::vec xi_g,
    arma::vec omega_g
) {
  Rcpp::RNGScope scope;
  const int n = B.n_rows;
  const int K = B.n_cols;
  const int N_total = N + BURN;
  
  // ----------------------------------------------------------
  // vectorized SPDE scaling:
  // s_g = (kappa_g^2 + lambda)^(-alpha_g/2)
  // ----------------------------------------------------------
  arma::vec s_g = arma::pow(kappa_g * kappa_g + lambda, -0.5 * alpha_g);
  
  // If B columns are orthonormal Laplacian eigenvectors,
  // || s_gk * u_k ||^2 = s_gk^2.
  // This avoids recomputing column norms every iteration.
  arma::vec x_norm2 = arma::square(s_g);
  
  // ----------------------------------------------------------
  // vectorized initialization
  // ----------------------------------------------------------
  arma::vec t_g = xi_g - xi_g0;
  arma::vec act_g = T_fun_np_vec(t_g, prior_type, lam1, lam2, lam3);
  arma::vec theta_g = s_g % act_g % omega_g;
  
  arma::vec one_n(n, fill::ones);
  arma::vec fit_g = B * theta_g;
  arma::vec res_g = y_g - mu_g * one_n - fit_g;
  
  // storage size
  const int n_keep = (N <= 0) ? 0 : ((N - 1) / thin + 1);
  
  arma::vec MU_G(n_keep, fill::zeros);
  arma::vec SIGMA2_G(n_keep, fill::zeros);
  arma::vec TAU2_G(n_keep, fill::zeros);
  arma::mat XI_G(K, n_keep, fill::zeros);
  arma::mat OMEGA_G(K, n_keep, fill::zeros);
  arma::mat THETA_G(K, n_keep, fill::zeros);
  arma::vec ACC_XI_G(K, fill::zeros);
  
  int save_idx = 0;
  
  for (int iter = 0; iter < N_total; ++iter) {
    
    // ========================================================
    // Step 1. update mu_g
    // mu_g | rest ~ N(v_mu_g * [1^T(y_g - B theta_g)/sigma2_g + nu_g0/eta_g0_sq], v_mu_g)
    // Since y_g - B theta_g = res_g + mu_g * 1_n
    // ========================================================
    {
      double v_mu_g = 1.0 / ( (double)n / sigma2_g + 1.0 / eta_g0_sq );
      double sum_y_minus_Btheta = arma::sum(res_g) + n * mu_g;
      double m_mu_g = v_mu_g * ( sum_y_minus_Btheta / sigma2_g + nu_g0 / eta_g0_sq );
      
      double mu_g_new = R::rnorm(m_mu_g, std::sqrt(v_mu_g));
      
      // res_new = y - mu_new 1 - B theta = res_old - (mu_new - mu_old) 1
      res_g -= (mu_g_new - mu_g);
      mu_g = mu_g_new;
    }
    
    // ========================================================
    // Step 2. update sigma_g^2
    // sigma_g^2 | rest ~ IG(m_g0 + n/2, gamma_g0 + ||res_g||^2 / 2)
    // ========================================================
    {
      double shape_sig = m_g0 + 0.5 * n;
      double rate_sig = gamma_g0 + 0.5 * arma::dot(res_g, res_g);
      sigma2_g = 1.0 / R::rgamma(shape_sig, 1.0 / rate_sig);
    }
    
    // ========================================================
    // Step 3. for k = 1,...,K:
    //         (a) marginalized MH for xi_gk
    //         (b) Gibbs for omega_gk
    // ========================================================
    for (int k = 0; k < K; ++k) {
      
      const arma::vec u_k = B.col(k);
      
      // remove current k contribution:
      // res_g becomes y_g - mu_g 1_n - sum_{l != k} u_l theta_gl
      res_g += u_k * theta_g(k);
      
      // x_{gk} = s_{gk} u_k
      // rx_gk = x_{gk}^T res_g = s_{gk} u_k^T res_g
      double rx_gk = s_g(k) * arma::dot(u_k, res_g);
      
      // ------------------------------------------------------
      // Step 3(a). marginalized MH for xi_gk
      //
      // This is the same NPrior-style target structure:
      //   -0.5 log( ||x_gk||^2 T^2 + sigma2_g/tau2_g )
      //   -0.5 xi_gk^2
      //   +0.5 rx_gk^2 T^2 / ( ||x_gk||^2 T^2 + sigma2_g/tau2_g ) / sigma2_g
      //
      // with x_gk = s_gk u_k.
      // ------------------------------------------------------
      for (int rep = 0; rep < mh_per_k; ++rep) {
        double xi_curr = xi_g(k);
        double t_curr = xi_curr - xi_g0;
        double act_curr = T_fun_np_scalar(t_curr, prior_type, lam1, lam2, lam3);
        
        double xi_cand = R::rnorm(xi_curr, xi_prop_sd);
        double t_cand = xi_cand - xi_g0;
        double act_cand = T_fun_np_scalar(t_cand, prior_type, lam1, lam2, lam3);
        
        double denom_curr = x_norm2(k) * act_curr * act_curr + sigma2_g / tau2_g;
        double denom_cand = x_norm2(k) * act_cand * act_cand + sigma2_g / tau2_g;
        
        double log_curr =
          -0.5 * std::log(denom_curr)
          -0.5 * xi_curr * xi_curr
          +0.5 * rx_gk * rx_gk * act_curr * act_curr / denom_curr / sigma2_g;
          
          double log_cand =
          -0.5 * std::log(denom_cand)
          -0.5 * xi_cand * xi_cand
          +0.5 * rx_gk * rx_gk * act_cand * act_cand / denom_cand / sigma2_g;
          
          double log_ratio = log_cand - log_curr;
          if (std::log(R::runif(0.0, 1.0)) < log_ratio) {
            xi_g(k) = xi_cand;
            act_g(k) = act_cand;
            ACC_XI_G(k) += 1.0;
          }
      }
      
      // ------------------------------------------------------
      // Step 3(b). Gibbs for omega_gk
      //
      // a_gk = s_gk * T(xi_gk - xi_g0)
      //
      // omega_gk | rest ~ N(nu_wgk, Sigma_wgk)
      // Sigma_wgk = [ a_gk^2 ||u_k||^2 / sigma_g^2 + 1/tau_g^2 ]^{-1}
      // Because ||u_k||^2 = 1 for orthonormal eigenvectors, use that directly.
      // ------------------------------------------------------
      {
        double a_gk = s_g(k) * act_g(k);
        double Sigma_wgk = 1.0 / ( (a_gk * a_gk) / sigma2_g + 1.0 / tau2_g );
        double nu_wgk = Sigma_wgk * ( a_gk * arma::dot(u_k, res_g) / sigma2_g );
        
        omega_g(k) = R::rnorm(nu_wgk, std::sqrt(Sigma_wgk));
        theta_g(k) = a_gk * omega_g(k);
      }
      
      // add updated k contribution back
      res_g -= u_k * theta_g(k);
    }
    
    // ========================================================
    // Step 4. update tau_g^2
    // tau_g^2 | omega_g ~ IG(a_g0 + K/2, b_g0 + ||omega_g||^2/2)
    // ========================================================
    if (tau2_g_update == 1) {
      double shape_tau = a_g0 + 0.5 * K;
      double rate_tau  = b_g0 + 0.5 * arma::dot(omega_g, omega_g);
      tau2_g = 1.0 / R::rgamma(shape_tau, 1.0 / rate_tau);
    }
    
    // save
    if (iter >= BURN && ((iter - BURN) % thin == 0)) {
      MU_G(save_idx) = mu_g;
      SIGMA2_G(save_idx) = sigma2_g;
      TAU2_G(save_idx) = tau2_g;
      XI_G.col(save_idx) = xi_g;
      OMEGA_G.col(save_idx) = omega_g;
      THETA_G.col(save_idx) = theta_g;
      save_idx++;
    }
  }
  
  ACC_XI_G /= (double)(N_total * mh_per_k);
  
  return Rcpp::List::create(
    Rcpp::Named("MU_G") = MU_G,
    Rcpp::Named("SIGMA2_G") = SIGMA2_G,
    Rcpp::Named("TAU2_G") = TAU2_G,
    Rcpp::Named("XI_G") = XI_G,
    Rcpp::Named("OMEGA_G") = OMEGA_G,
    Rcpp::Named("THETA_G") = THETA_G,
    Rcpp::Named("ACC_XI_G") = ACC_XI_G,
    Rcpp::Named("s_g") = s_g
  );
}