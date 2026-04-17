prior_to_type <- function(prior,
                          act_type = NULL,
                          SpSL_type = NULL) {
  if (prior == "BL") return(1L)
  if (prior == "HS") return(2L)
  if (prior == "SpSL-L") return(3L)
  if (prior == "SpSL-C") return(4L)
  
  if (prior == "custom") {
    if (is.null(act_type) || is.null(SpSL_type)) {
      stop("for prior = 'custom', act_type and SpSL_type must be provided")
    }
    
    if (act_type == 1L && SpSL_type == 0L) return(5L)
    if (act_type == 2L && SpSL_type == 0L) return(6L)
    if (act_type == 1L && SpSL_type == 1L) return(7L)
    if (act_type == 2L && SpSL_type == 1L) return(8L)
    
    stop("unsupported combination of act_type and SpSL_type")
  }
  
  stop("unsupported prior")
}


fit_one_gene <- function(
    y_g, B, lambda, N = 2000, BURN = 1000, thin = 1, mh_per_k = 1,
    
    kappa_g = 118, alpha_g = 2.5, xi_g0 = 0,
    
    nu_g0 = 0, eta_g0_sq = 100, m_g0 = 2, gamma_g0 = 1, a_g0 = 1, b_g0 = 1,
    
    prior = c("BL", "HS", "SpSL-L", "SpSL-C", "custom"),
    act_type = NULL, SpSL_type = NULL, lam = c(0, 0, 0),
    
    xi_prop_sd = 0.15,
    tau2_g_update = FALSE, tau2_g_init = 1,
    mu_g_init = mean(y_g), sigma2_g_init = max(var(y_g), 1e-6),
    xi_g_init = NULL, omega_g_init = NULL
) {
  prior <- match.arg(prior)
  
  prior_type <- prior_to_type(prior, act_type, SpSL_type)
  
  if (length(lam) != 3L) {
    stop("lam must be a numeric vector of length 3")
  }
  lam <- as.numeric(lam)
  lam1 <- lam[1]
  lam2 <- lam[2]
  lam3 <- lam[3]
  
  K <- ncol(B)
  
  if (prior %in% c("BL", "HS")) {
    xi_g0 <- 0
  }
  
  if (is.null(xi_g_init)) {
    xi_g_init <- rep(0, K)
  }
  if (is.null(omega_g_init)) {
    omega_g_init <- rep(0, K)
  }
  
  spde_neuronized_basic_one_gene_cpp(
    y_g = y_g, B = as.matrix(B), lambda = as.numeric(lambda),
    N = N, BURN = BURN, thin = thin, mh_per_k = mh_per_k,
    
    kappa_g = kappa_g, alpha_g = alpha_g, xi_g0 = xi_g0,
    
    nu_g0 = nu_g0, eta_g0_sq = eta_g0_sq, m_g0 = m_g0, gamma_g0 = gamma_g0,
    a_g0 = a_g0, b_g0 = b_g0,
    
    prior_type = as.integer(prior_type),
    lam1 = lam1, lam2 = lam2, lam3 = lam3, xi_prop_sd = xi_prop_sd,
    
    tau2_g_update = as.integer(tau2_g_update),
    mu_g = mu_g_init, sigma2_g = sigma2_g_init,
    tau2_g = tau2_g_init, xi_g = xi_g_init, omega_g = omega_g_init
  )
}

fit_all_genes <- function(
    Y, S, K = 50, knn = 10, h = NULL, 
    symmetric = TRUE, normalized_laplacian = FALSE,
    N = 2000, BURN = 1000, thin = 1, mh_per_k = 1,
    kappa_g = 118, alpha_g = 2.5, xi_g0 = 0,
    nu_g0 = 0, eta_g0_sq = 100, m_g0 = 2, gamma_g0 = 1, a_g0 = 1, b_g0 = 1,
    prior = c("BL", "HS", "SpSL-L", "SpSL-C", "custom"),
    act_type = NULL, SpSL_type = NULL, lam = c(0,0,0),
    xi_prop_sd = 0.15, tau2_g_update = FALSE, tau2_g_init = 1,
    ncores = 1L
) {
  stopifnot(nrow(Y) == nrow(S))
  prior <- match.arg(prior)
  
  W <- build_W_knn(S, k = knn, h = h, symmetric = symmetric)
  L <- compute_laplacian(W, normalized = normalized_laplacian)
  eig <- compute_eigen(L, K = K)
  
  B_basis <- eig$B
  lambda <- eig$lambda
  G <- ncol(Y)
  
  expand_to_G <- function(x) {
    if (length(x) == 1L) rep(x, G) else x
  }
  
  kappa_g_vec <- expand_to_G(kappa_g)
  alpha_g_vec <- expand_to_G(alpha_g)
  xi_g0_vec <- expand_to_G(xi_g0)
  tau2_g_init_vec <- expand_to_G(tau2_g_init)
  
  fit_one <- function(g) {
    fit_one_gene(
      y_g = Y[, g], B = B_basis, lambda = lambda,
      N = N, BURN = BURN, thin = thin, mh_per_k = mh_per_k,
      kappa_g = kappa_g_vec[g], alpha_g = alpha_g_vec[g], xi_g0 = xi_g0_vec[g],
      nu_g0 = nu_g0, eta_g0_sq = eta_g0_sq, m_g0 = m_g0, gamma_g0 = gamma_g0,
      a_g0 = a_g0, b_g0 = b_g0,
      prior = prior, act_type = act_type, SpSL_type = SpSL_type, lam = lam,
      xi_prop_sd = xi_prop_sd,
      tau2_g_update = tau2_g_update, tau2_g_init = tau2_g_init_vec[g]
    )
  }
  
  if (ncores <= 1L) {
    fits <- lapply(seq_len(G), fit_one)
  } else {
    fits <- mclapply(seq_len(G), fit_one, mc.cores = ncores)
  }
  
  list(
    W = W,
    L = L,
    lambda = lambda,
    B = B_basis,
    fits = fits,
    prior = prior
  )
}