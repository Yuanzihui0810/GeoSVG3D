#' GeoSVG3D package
#'
#' Geometry-aware Bayesian detection of spatially variable genes in
#' three-dimensional spatial transcriptomics data.
#'
#' @useDynLib GeoSVG3D, .registration = TRUE
#' @importFrom Rcpp evalCpp
#' @importFrom utils setTxtProgressBar
#' @keywords internal
"_PACKAGE"

#' Convert a prior specification to its internal integer code
#'
#' @param prior Character string. One of "BL", "HS", "SpSL-L", "SpSL-C",
#'   "SpSL-G", or "custom".
#' @param act_type Integer, only used when `prior = "custom"`. Activation
#'   type code (1 or 2).
#' @param SpSL_type Integer, only used when `prior = "custom"`. Spike-and-slab
#'   type code (0 or 1).
#'
#' @return An integer prior type code used internally by the C++ sampler.
#' @noRd
prior_to_type <- function(prior,
                          act_type = NULL,
                          SpSL_type = NULL) {
  if (prior == "BL") return(1L)
  if (prior == "HS") return(2L)
  if (prior == "SpSL-L") return(3L)
  if (prior == "SpSL-C") return(4L)
  if (prior == "SpSL-G") return(5L)
  
  if (prior == "custom") {
    if (is.null(act_type) || is.null(SpSL_type)) {
      stop("for prior = 'custom', act_type and SpSL_type must be provided")
    }
    
    if (act_type == 1L && SpSL_type == 0L) return(6L)
    if (act_type == 2L && SpSL_type == 0L) return(7L)
    if (act_type == 1L && SpSL_type == 1L) return(8L)
    if (act_type == 2L && SpSL_type == 1L) return(9L)
    
    stop("unsupported combination of act_type and SpSL_type")
  }
  
  stop("unsupported prior")
}
#' Scale three-dimensional spatial coordinates to the unit cube

#'

#' Each coordinate column is independently transformed to the interval from 0 to 1.
#'
#' @param S A numeric matrix with observations in rows and three spatial
#'   coordinates in columns.
#'
#' @return A numeric matrix with the same dimensions, row names, and column
#'   names as `S`, with every coordinate column scaled to `[0, 1]`.
#'
#' @export
scale_spatial_coords_01 <- function(S) {
  S <- as.matrix(S)
  storage.mode(S) <- "double"
  
  if (ncol(S) != 3) {
    stop("S must have 3 columns.")
  }
  
  if (is.null(colnames(S))) {
    colnames(S) <- c("x", "y", "z")
  }
  
  scale01 <- function(v) {
    r <- range(v, na.rm = TRUE)
    if (r[2] == r[1]) {
      return(rep(0, length(v)))
    }
    (v - r[1]) / (r[2] - r[1])
  }
  
  S_scaled <- apply(S, 2, scale01)
  S_scaled <- as.matrix(S_scaled)
  colnames(S_scaled) <- colnames(S)
  rownames(S_scaled) <- rownames(S)
  
  return(S_scaled)
}
#' Fit the GeoSVG3D model to one gene
#'
#' @param y_g Numeric expression vector for one gene.
#' @param B Numeric graph-basis matrix.
#' @param lambda Numeric vector of graph Laplacian eigenvalues.
#' @param N Number of retained MCMC iterations after burn-in.
#' @param BURN Number of burn-in iterations.
#' @param thin Thinning interval.
#' @param mh_per_k Number of Metropolis-Hastings updates per basis coefficient
#'   in each MCMC iteration.
#' @param kappa_g Spectral range parameter.
#' @param alpha_g Spectral smoothness parameter.
#' @param xi_g0 Threshold parameter in the neuronized prior.
#' @param nu_g0 Prior mean for the gene-specific intercept.
#' @param eta_g0_sq Prior variance for the gene-specific intercept.
#' @param m_g0 Shape parameter of the inverse-gamma prior for residual
#'   variance.
#' @param gamma_g0 Scale parameter of the inverse-gamma prior for residual
#'   variance.
#' @param a_g0 Shape parameter of the inverse-gamma prior for `tau2_g`.
#' @param b_g0 Scale parameter of the inverse-gamma prior for `tau2_g`.
#' @param prior Character string specifying the neuronized prior.
#' @param act_type Optional custom activation type.
#' @param SpSL_type Optional custom continuous or thresholded prior type.
#' @param lam Numeric vector of length three containing custom activation
#'   parameters.
#' @param xi_prop_sd Standard deviation of the random-walk proposal for
#'   `xi_g`.
#' @param tau2_g_update Logical; whether to update `tau2_g`.
#' @param tau2_g_init Initial value of `tau2_g`.
#' @param mu_g_init Optional initial value of the intercept.
#' @param sigma2_g_init Optional initial value of the residual variance.
#' @param xi_g_init Optional initial activation vector.
#' @param omega_g_init Optional initial coefficient-magnitude vector.
#'
#' @return A list of posterior samples returned by the C++ MCMC sampler.
#' @noRd
fit_one_gene <- function(
    y_g, B, lambda, N = 2000, BURN = 2000, thin = 1, mh_per_k = 1,
    
    kappa_g = 0, alpha_g = 10, xi_g0 = 0,
    
    nu_g0 = 0, eta_g0_sq = 100, m_g0 = 2, gamma_g0 = 1,
    a_g0 = 1, b_g0 = 1,
    
    prior = c("BL", "HS", "SpSL-L", "SpSL-C", "SpSL-G", "custom"),
    act_type = NULL, SpSL_type = NULL, lam = c(0, 0, 0),
    
    xi_prop_sd = 0.15,
    tau2_g_update = FALSE,
    tau2_g_init = 1,
    
    mu_g_init = NULL,
    sigma2_g_init = NULL,
    xi_g_init = NULL,
    omega_g_init = NULL
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
  
  if (is.null(mu_g_init)) {
    mu_g_init <- mean(y_g)
  }
  
  if (is.null(sigma2_g_init)) {
    sigma2_g_init <- max(stats::var(y_g), 1e-6)
  }
  
  if (is.null(xi_g_init)) {
    xi_g_init <- numeric(K)
  }
  
  if (is.null(omega_g_init)) {
    omega_g_init <- numeric(K)
  }
  
  if (prior %in% c("BL", "HS")) {
    xi_g0 <- 0
  }
  
  spde_neuronized_basic_one_gene_cpp(
    y_g = y_g, B = B, lambda = lambda,
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

#' Fit the GeoSVG3D model to all genes
#'
#' Constructs a spatial k-nearest-neighbor graph, computes graph Laplacian
#' basis functions, and fits the gene-specific Bayesian model independently
#' to every column of the expression matrix.
#'
#' @param Y A numeric or matrix-like expression object with observations in
#'   rows and genes in columns.
#' @param S A numeric matrix with observations in rows and three spatial
#'   coordinates in columns.
#' @param K Number of graph Laplacian basis functions.
#' @param knn Number of nearest neighbors used to construct the spatial graph.
#' @param h Optional Gaussian-kernel bandwidth passed to `build_W_knn()`.
#' @param normalized_laplacian Logical; whether to use the normalized graph
#'   Laplacian.
#' @param binary Logical; whether to use binary instead of Gaussian-kernel
#'   graph weights.
#' @param N Number of retained MCMC iterations after burn-in.
#' @param BURN Number of burn-in iterations.
#' @param thin Thinning interval.
#' @param mh_per_k Number of Metropolis-Hastings updates per basis coefficient
#'   in each MCMC iteration.
#' @param kappa_g Scalar or gene-specific vector of spectral range parameters.
#' @param alpha_g Scalar or gene-specific vector of spectral smoothness
#'   parameters.
#' @param xi_g0 Scalar or gene-specific vector of neuronized-prior thresholds.
#' @param nu_g0 Prior mean for each gene-specific intercept.
#' @param eta_g0_sq Prior variance for each gene-specific intercept.
#' @param m_g0 Shape parameter of the inverse-gamma prior for residual
#'   variance.
#' @param gamma_g0 Scale parameter of the inverse-gamma prior for residual
#'   variance.
#' @param a_g0 Shape parameter of the inverse-gamma prior for `tau2_g`.
#' @param b_g0 Scale parameter of the inverse-gamma prior for `tau2_g`.
#' @param prior Character string specifying the neuronized prior.
#' @param act_type Optional custom activation type.
#' @param SpSL_type Optional custom continuous or thresholded prior type.
#' @param lam Numeric vector of length three containing custom activation
#'   parameters.
#' @param xi_prop_sd Standard deviation of the random-walk proposal for
#'   `xi_g`.
#' @param tau2_g_update Logical; whether to update `tau2_g`.
#' @param tau2_g_init Scalar or gene-specific vector of initial `tau2_g`
#'   values.
#' @param ncores Number of parallel workers.
#' @param seed Integer seed used to initialize gene-specific random-number
#'   streams in parallel computation.
#'
#' @return A list containing the adjacency matrix, graph Laplacian,
#'   eigenvalues, basis matrix, gene-specific spectral parameters, posterior
#'   fits, and prior specification.
#'
#' @export
fit_all_genes <- function(
    Y, S, K = 5, knn = 100, h = NULL, 
    normalized_laplacian = FALSE, binary = FALSE,
    N = 2000, BURN = 2000, thin = 1, mh_per_k = 1,
    
    kappa_g = 0,
    alpha_g = 10, xi_g0 = 0,
    nu_g0 = 0, eta_g0_sq = 100, m_g0 = 2, gamma_g0 = 1, a_g0 = 1, b_g0 = 1,
    prior = c("BL", "HS", "SpSL-L", "SpSL-C", "SpSL-G", "custom"),
    act_type = NULL, SpSL_type = NULL, lam = c(0,0,0),
    xi_prop_sd = 0.15, tau2_g_update = FALSE, tau2_g_init = 1,
    ncores = 1L, seed = 1
) {
  stopifnot(nrow(Y) == nrow(S))
  S <- scale_spatial_coords_01(S)
  
  prior <- match.arg(prior)
  W <- build_W_knn(S, k = knn, h = h, binary = binary)
  L <- compute_laplacian(W, normalized = normalized_laplacian)
  d <- Matrix::rowSums(W)
  # cat(paste("Check graph degree...", sum(d == 0), "node with no edge\n"))
  eig <- compute_eigen(L, K = K, tol = 1e-10)
  # cat("Finish building graph\n")
  
  B_basis <- eig$U
  lambda <- eig$lambda
  
  Y <- as.matrix(Y)
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
  
  cat("Start analysing...\n")
  
  if (ncores <= 1L) {
    # fits <- lapply(seq_len(G), fit_one)
    fits <- vector("list", G)
    pb <- utils::txtProgressBar(min = 0, max = G, style = 3)
    on.exit(close(pb), add = TRUE)
    
    for (g in seq_len(G)) {
      fits[[g]] <- fit_one(g)
      setTxtProgressBar(pb, g)
    }
  } else {
    # t1 = Sys.time()
    cl <- parallel::makeCluster(ncores)
    # parallel::clusterSetRNGStream(cl, iseed = seed)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    # t2 = Sys.time()
    # print(t2-t1)
    
    # ***** Later change to the following code *****
    pkg <- "GeoSVG3D"  # Package name

    parallel::clusterExport(cl, "pkg", envir = environment())

    parallel::clusterEvalQ(cl, {
      library(pkg, character.only = TRUE)
      NULL
    })
    
    parallel::clusterExport(
      cl,
      varlist = c(
        #"GeoSVG3D",#"pkg_name", #"fit_one_gene", "prior_to_type",
        "Y", "B_basis", "lambda",
        "N", "BURN", "thin", "mh_per_k",
        "kappa_g_vec", "alpha_g_vec", "xi_g0_vec",
        "nu_g0", "eta_g0_sq", "m_g0", "gamma_g0", "a_g0", "b_g0",
        "prior", "act_type", "SpSL_type", "lam",
        "xi_prop_sd", "tau2_g_update", "tau2_g_init_vec",
        "seed"
      ),
      envir = environment()
    )
    
    fits <- parallel::parLapply(
      cl,
      seq_len(G),
      function(g) {
        set.seed(seed + g)
        
        fit_fun <- get(
          "fit_one_gene",
          envir = asNamespace(pkg),
          inherits = FALSE
        )
        
        fit_fun(
          y_g = Y[, g],
          B = B_basis,
          lambda = lambda,
          N = N, BURN = BURN, thin = thin, mh_per_k = mh_per_k,
          kappa_g = kappa_g_vec[g],
          alpha_g = alpha_g_vec[g],
          xi_g0 = xi_g0_vec[g],
          nu_g0 = nu_g0, eta_g0_sq = eta_g0_sq,
          m_g0 = m_g0, gamma_g0 = gamma_g0,
          a_g0 = a_g0, b_g0 = b_g0,
          prior = prior, act_type = act_type, SpSL_type = SpSL_type, lam = lam,
          xi_prop_sd = xi_prop_sd,
          tau2_g_update = tau2_g_update,
          tau2_g_init = tau2_g_init_vec[g]
        )
      }
    )
    
    # fits <- parallel::parLapply(cl, seq_len(G), function(g) {
    #   set.seed(seed+g)
    #   fit_one_gene(
    #     y_g = Y[, g],
    #     B = B_basis,
    #     lambda = lambda,
    #     N = N, BURN = BURN, thin = thin, mh_per_k = mh_per_k,
    #     kappa_g = kappa_g_vec[g],
    #     alpha_g = alpha_g_vec[g],
    #     xi_g0 = xi_g0_vec[g],
    #     nu_g0 = nu_g0, eta_g0_sq = eta_g0_sq,
    #     m_g0 = m_g0, gamma_g0 = gamma_g0,
    #     a_g0 = a_g0, b_g0 = b_g0,
    #     prior = prior, act_type = act_type, SpSL_type = SpSL_type, lam = lam,
    #     xi_prop_sd = xi_prop_sd,
    #     tau2_g_update = tau2_g_update,
    #     tau2_g_init = tau2_g_init_vec[g]
    #   )
    # })
  }
  
  list(
    W = W,
    L = L,
    lambda = lambda,
    kappa_g = kappa_g_vec,
    B = B_basis,
    fits = fits,
    prior = prior
  )
}

#' Compute posterior probabilities of no spatial effect
#'
#' For each gene, calculates the posterior proportion of MCMC samples in
#' which all spatial-basis coefficients are zero, or within a specified
#' tolerance of zero.
#'
#' @param fit_res An object returned by `fit_all_genes()`.
#' @param tol Nonnegative numerical tolerance used to define a zero
#'   coefficient.
#'
#' @return A numeric vector containing one posterior null probability per gene.
#'
#' @export
compute_p0g <- function(fit_res, tol = 0) {
  G <- length(fit_res$fits)
  p0g <- numeric(G)
  
  for (g in seq_len(G)) {
    Theta <- fit_res$fits[[g]]$THETA_G   # K x M
    
    if (tol == 0) {
      is_H0_each_iter <- apply(Theta == 0, 2, all)
    } else {
      is_H0_each_iter <- apply(abs(Theta) <= tol, 2, all)
    }
    
    p0g[g] <- mean(is_H0_each_iter)
  }
  
  return(p0g)
}

#' Select spatially variable genes using Bayesian FDR control
#'
#' @param p0g Numeric vector containing posterior probabilities of no spatial
#'   effect.
#' @param alpha Target Bayesian false discovery rate.
#'
#' @return A list containing the binary selection vector, selected gene
#'   indices, ordering information, cumulative Bayesian FDR values, and the
#'   selected cutoff.
#'
#' @export
select_svg_bfdr <- function(p0g, alpha = 0.05) {
  ord <- order(p0g)
  p_sorted <- p0g[ord]
  cum_bfdr <- cumsum(p_sorted) / seq_along(p_sorted)
  
  idx <- which(cum_bfdr <= alpha)
  k <- if (length(idx) == 0) 0L else max(idx)
  
  pred <- rep(0L, length(p0g))
  selected <- integer(0)
  
  if (k > 0) {
    selected <- ord[seq_len(k)]
    pred[selected] <- 1L
  }
  
  list(
    pred = pred,
    selected = selected,
    k = k,
    ord = ord,
    p_sorted = p_sorted,
    cum_bfdr = cum_bfdr,
    nu_star = if (k > 0) p_sorted[k] else NA_real_
  )
}

#' Evaluate spatially variable gene detection
#'
#' @param truth Binary vector of true SVG labels.
#' @param pred Binary vector of predicted SVG labels.
#'
#' @return A list containing a one-row summary of detection metrics and a
#'   two-by-two confusion matrix.
#'
#' @export
evaluate_svg_detection <- function(truth, pred) {
  truth <- as.integer(truth)
  pred <- as.integer(pred)
  
  TP <- sum(truth == 1 & pred == 1)
  TN <- sum(truth == 0 & pred == 0)
  FP <- sum(truth == 0 & pred == 1)
  FN <- sum(truth == 1 & pred == 0)
  
  correct_n <- TP + TN
  wrong_n <- FP + FN
  
  accuracy <- (TP + TN) / length(truth)
  power <- if ((TP + FN) == 0) NA else TP / (TP + FN)
  precision <- if ((TP + FP) == 0) NA else TP / (TP + FP)
  
  specificity <- if ((TN + FP) == 0) NA else TN / (TN + FP)
  
  recall <- power
  F1 <- if (is.na(precision) || is.na(recall) || (precision + recall) == 0) {
    NA
  } else {
    2 * precision * recall / (precision + recall)
  }
  
  confusion <- matrix(
    c(TP, FN, FP, TN),
    nrow = 2,
    byrow = TRUE,
    dimnames = list(
      Truth = c("SVG", "non-SVG"),
      Pred = c("SVG", "non-SVG")
    )
  )
  
  stats <- data.frame(
    correct_n = correct_n,
    wrong_n = wrong_n,
    accuracy = accuracy,
    power = power,
    specificity = specificity,
    precision = precision,
    recall = recall,
    F1 = F1
  )
  
  list(
    stats = stats,
    confusion = confusion
  )
}

#' Run posterior null-probability estimation and Bayesian FDR selection
#'
#' Computes posterior null probabilities, performs Bayesian FDR selection,
#' and evaluates the selected genes against known truth labels.
#'
#' @param fit_res An object returned by `fit_all_genes()`.
#' @param truth Binary vector of true SVG labels used for performance
#'   evaluation.
#' @param alpha Target Bayesian false discovery rate.
#' @param tol Nonnegative numerical tolerance used to define a zero
#'   coefficient.
#'
#' @return A list containing posterior null probabilities, Bayesian FDR
#'   selection results, and detection-performance summaries.
#'
#' @export
run_svg_bfdr_pipeline <- function(fit_res, truth, alpha = 0.05, tol = 0) {
  p0g <- compute_p0g(fit_res, tol = tol)
  bfdr_res <- select_svg_bfdr(p0g, alpha = alpha)
  eval_res <- evaluate_svg_detection(truth, bfdr_res$pred)
  
  list(
    p0g = p0g,
    bfdr = bfdr_res,
    evaluation = eval_res
  )
}

