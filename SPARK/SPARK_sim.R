library(SPARK)
args = commandArgs(trailingOnly = T) # contain inputs from command line
args = as.numeric(args)
seed = args[1]
# ============================================================
# 1. 3D version of spark.test
# ============================================================

spark.test.new <- function(object, 
                           kernel_mat = NULL, 
                           check_positive = TRUE, 
                           verbose = TRUE) {
  cat("\n### Running spark.test.new with 3D coordinates ###\n")
  
  ## the variables input
  if(!is.null(kernel_mat) & !is.list(kernel_mat)){ stop("SPARK.TEST::kernel_mat must be a list, please provide list type!")}
  ## testing one kernel at a time
  res_kernel <- list()
  res_pval <- NULL
  if(is.null(kernel_mat)){ #calculate kernel matrix using location information
    # Euclid distance, and compute the range of parameter l
    ED <- as.matrix(dist(object@location[ ,1:3]))
    lrang <- ComputeGaussianPL(ED, compute_distance=FALSE)[3:7]
    ## total 10 kernels, i.e., 5 gaussian and 5 periodic kernels
    for(ikernel in c(1:5) ){
      # Gaussian kernel
      cat(paste0("## testing Gaussian kernel: ",ikernel,"...\n"))
      kernel_mat <- exp(-ED^2/(2*lrang[ikernel]^2))
      object <- spark.test_each(object, kernel_mat=kernel_mat, check_positive=check_positive, verbose=verbose)
      res_pval <- cbind(res_pval, object@res_stest$sw)
      res_kernel <- c(res_kernel, object@res_stest)
      rm(kernel_mat)
      
      # Periodic kernel
      cat(paste0("## testing Periodic kernel: ",ikernel,"...\n"))
      kernel_mat <- cos(2*pi*ED/lrang[ikernel])
      object <- spark.test_each(object, kernel_mat=kernel_mat, check_positive=check_positive, verbose=verbose)
      res_pval <- cbind(res_pval, object@res_stest$sw)
      res_kernel <- c(res_kernel, object@res_stest)
      rm(kernel_mat)
    }# end for ikernel
    colnames(res_pval) <- paste0(c("GSP","COS"), rep(1:5,each=2))
    rownames(res_pval) <- rownames(object@counts)
    print(paste("rownames(object@counts):", rownames(object@counts)))
    print(paste("rownames(res_pval):", rownames(res_pval)))
  }else{# kernel_mat is a list provided by user
    for(ikernel in 1:length(kernel_mat) ){
      # pre-defined kernels
      cat(paste0("## testing pre-defined kernel: ",ikernel, "...\n"))
      object <- spark.test_each(object, kernel_mat=kernel_mat[[ikernel]], check_positive=check_positive, verbose=verbose)
      res_pval <- cbind(res_pval, object@res_stest$sw)
      res_kernel <- c(res_kernel, object@res_stest)
    }# end for ikernel
    colnames(res_pval) <- paste0("kernel", 1:length(kernel_mat) )
    print(paste("rownames(object@counts):", rownames(object@counts)))
    print(paste("rownames(res_pval):", rownames(res_pval)))
    rownames(res_pval) <- rownames(object@counts)
  }# end fi
  
  object@res_stest <- res_kernel
  ## integrate ten p-values into one
  num_pval <- ncol(res_pval)
  num_gene <- nrow(res_pval)
  ## compute the weight to p-values integrate
  if(is.null(object@weights)){
    weights <- matrix(rep(1.0/num_pval, num_pval*num_gene), ncol=num_pval )
  }else if(!is.matrix(object@weights)){
    weights <- as.matrix(object@weights)
  }else {
    weights <- object@weights
  }# end fi
  combined_pvalue <- CombinePValues(res_pval, weights=weights)
  object@res_mtest <- data.frame(res_pval, combined_pvalue = combined_pvalue,  adjusted_pvalue = p.adjust(combined_pvalue, method="BY") )
  # return the results
  return(object)
}# end function score test


# ============================================================
# 2. Settings
# ============================================================

# seed <- 1
p_cutoff <- 0.01

data_dir <- "/public3/home/scg5453/zhyuan/GeoSVG_sim/sim2"

out_dir <- "/public3/home/scg5453/zhyuan/GeoSVG_SPARK_output/sim2"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

sim_rds_file <- file.path(
  data_dir,
  paste0(
    "simulation2_data_seed",
    seed,
    ".rds"
  )
)

spark_rds_file <- file.path(
  out_dir,
  paste0(
    "spark_simulation2_3D_results_seed",
    seed,
    ".rds"
  )
)

spark_csv_file <- file.path(
  out_dir,
  paste0(
    "spark_simulation2_3D_table_seed",
    seed,
    ".csv"
  )
)

simulate_svg_continuous_data <- function(
    n = 10000,
    G = 1000,
    n_svg = 120,
    sigma_eps = 0.01,              # noise for SVG genes
    sigma_eps_non_svg = 0.01,      # stronger noise for non-SVG genes
    eps_cap = 5,                   # truncate log-scale noise to avoid extreme exp()
    zero_rate = 0,
    zero_rate_svg = zero_rate,
    zero_rate_non_svg = zero_rate,
    mu_mean = 0,
    mu_sd = 1,
    amp_gradient = 2,    # 1
    amp_periodic = 2,    # 1
    amp_gaussian = 2,    # 1
    amp_composite = 2,   # 1
    amp_granularity = 2, # 1
    seed = 123
) {
  seed = seed
  set.seed(seed)
  stopifnot(n_svg <= G)
  stopifnot(n_svg <= 120)
  stopifnot(sigma_eps >= 0)
  stopifnot(sigma_eps_non_svg >= 0)
  stopifnot(zero_rate_svg >= 0, zero_rate_svg <= 1)
  stopifnot(zero_rate_non_svg >= 0, zero_rate_non_svg <= 1)
  
  # ============================================================
  # 1. coordinates in [0, 1]^3
  # ============================================================
  # Intrinsic coordinates
  u <- runif(n, 0, 1)
  v <- runif(n, 0, 1)
  w <- runif(n, 0, 1)
  
  # Full angular range
  theta <- 2 * pi * v
  
  # Preserve the ring-like thickness of the original geometry
  rho <- 0.65 + 0.35 * u
  
  # Semi-major and semi-minor axes
  a <- 1.30
  b <- 0.80
  
  # Elliptical coordinates
  x_obs <- a * rho * cos(theta)
  y_obs <- b * rho * sin(theta)
  z_obs <- w
  
  # Rescale observed coordinates to [0, 1]
  rescale01 <- function(a) {
    (a - min(a)) / (max(a) - min(a))
  }
  
  S <- cbind(
    x = rescale01(x_obs),
    y = rescale01(y_obs),
    z = rescale01(z_obs)
  )
  
  # Intrinsic coordinates used to define spatial expression patterns
  x <- u
  y <- v
  z <- w
  
  xc <- x - 0.5
  yc <- y - 0.5
  zc <- z - 0.5
  
  # ------------------------------------------------------------
  # Helper functions for 3D pattern construction
  # ------------------------------------------------------------
  z_warp1 <- z + 0.10 * sin(2 * pi * x) + 0.08 * cos(2 * pi * y)
  z_warp2 <- z + 0.10 * sin(2 * pi * (x + y))
  z_warp3 <- z + 0.08 * cos(2 * pi * x) - 0.08 * sin(2 * pi * y)
  
  r_xy <- sqrt(xc^2 + yc^2)
  r_3d <- sqrt(xc^2 + yc^2 + zc^2)
  theta_xy <- atan2(yc, xc)
  
  soft_step <- function(v, width = 0.06) {
    tanh(v / width)
  }
  
  sheet <- function(surface, width = 0.07) {
    exp(-0.5 * (surface / width)^2)
  }
  
  gspot3d <- function(cx, cy, cz, sx = 0.12, sy = 0.12, sz = 0.12) {
    exp(-0.5 * (((x - cx) / sx)^2 +
                  ((y - cy) / sy)^2 +
                  ((z - cz) / sz)^2))
  }
  
  shell3d <- function(cx, cy, cz, r0 = 0.25, width = 0.05) {
    r <- sqrt((x - cx)^2 + (y - cy)^2 + (z - cz)^2)
    exp(-0.5 * ((r - r0) / width)^2)
  }
  
  center_f <- function(f) {
    as.numeric(f - mean(f))
  }
  
  # ============================================================
  # 2. construct 120 spatial patterns
  #    1-30     3D gradient
  #    31-60    3D periodic
  #    61-90    3D Gaussian localized
  #    91-120   complex
  # ============================================================
  f_list <- list()
  f_type <- character(0)
  
  add_f <- function(f, type) {
    f_list[[length(f_list) + 1]] <<- center_f(f)
    f_type[length(f_list)] <<- type
  }
  
  # ------------------------------------------------------------
  # A. 30 3D gradient patterns
  # ------------------------------------------------------------
  angles <- seq(0, 2 * pi, length.out = 15 + 1)[-16]
  
  for (a in angles) {
    # 3D oblique gradient with mild z-axis bending
    dx1 <- cos(a)
    dy1 <- sin(a)
    dz1 <- 0.45 * sin(2 * a)
    norm1 <- sqrt(dx1^2 + dy1^2 + dz1^2)
    
    dx1 <- dx1 / norm1
    dy1 <- dy1 / norm1
    dz1 <- dz1 / norm1
    
    # f1 <- dx1 * xc + dy1 * yc + dz1 * zc
    f1 <- dx1 * x + dy1 * y + dz1 * z
    f1 <- f1 + 0.15 * sin(2 * pi * z_warp1)
    add_f(amp_gradient * f1, "gradient")
    
    dx2 <- -cos(a)
    dy2 <- -sin(a)
    dz2 <- 0.45 * cos(2 * a)
    norm2 <- sqrt(dx2^2 + dy2^2 + dz2^2)
    
    dx2 <- dx2 / norm2
    dy2 <- dy2 / norm2
    dz2 <- dz2 / norm2
    
    # f2 <- dx2 * xc + dy2 * yc + dz2 * zc
    f2 <- dx2 * x + dy2 * y + dz2 * z
    f2 <- f2 + 0.15 * cos(2 * pi * z_warp2)
    add_f(amp_gradient * f2, "gradient")
  }
  
  # ------------------------------------------------------------
  # B. 30 3D periodic patterns
  # ------------------------------------------------------------
  periodic_dirs3 <- rbind(
    c(1, 0, 0.5),
    c(0, 1, 0.5),
    c(1, 1, 0.5),
    c(1, -1, 0.5),
    c(2, 1, 0.7),
    c(1, 2, -0.7),
    c(2, -1, 0.7),
    c(1, -2, -0.7),
    c(1, 1, 1),
    c(1, -1, 1)
  )
  
  freqs <- c(1, 2, 3)
  
  for (j in seq_len(nrow(periodic_dirs3))) {
    d <- periodic_dirs3[j, ]
    d <- d / sqrt(sum(d^2))
    
    u <- d[1] * x + d[2] * y + d[3] * z
    u <- (u - min(u)) / (max(u) - min(u) #+ 1e-8
    )
    
    for (fr in freqs) {
      phase <- ((j + fr) %% 6) * pi / 6
      
      f <- sin(2 * pi * fr * u + phase +
                 0.25 * sin(2 * pi * z_warp1))
      
      add_f(amp_periodic * f, "periodic")
    }
  }
  
  # ------------------------------------------------------------
  # C. 30 3D Gaussian localized patterns
  # ------------------------------------------------------------
  centers3 <- as.matrix(expand.grid(
    cx = c(0.25, 0.50, 0.75),
    cy = c(0.25, 0.50, 0.75),
    cz = c(0.25, 0.50, 0.75)
  ))
  
  extra_centers3 <- matrix(
    c(
      0.20, 0.50, 0.50,
      0.80, 0.50, 0.50,
      0.50, 0.20, 0.50
    ),
    ncol = 3,
    byrow = TRUE
  )
  
  centers3 <- rbind(centers3, extra_centers3)
  centers3 <- centers3[seq_len(30), , drop = FALSE]
  
  sx_vec <- c(0.08, 0.10, 0.12, 0.15, 0.18)
  sy_vec <- c(0.10, 0.14, 0.18, 0.08, 0.12)
  sz_vec <- c(0.10, 0.12, 0.16, 0.08, 0.14)
  
  for (j in seq_len(30)) {
    cx <- centers3[j, 1]
    cy <- centers3[j, 2]
    cz <- centers3[j, 3]
    
    sx <- sx_vec[(j - 1) %% length(sx_vec) + 1]
    sy <- sy_vec[(j - 1) %% length(sy_vec) + 1]
    sz <- sz_vec[(j - 1) %% length(sz_vec) + 1]
    
    # Slightly curved 3D localized region
    cx_z <- cx + 0.04 * sin(2 * pi * z + j * pi / 10)
    cy_z <- cy + 0.04 * cos(2 * pi * z + j * pi / 10)
    
    f <- exp(-0.5 * (((x - cx_z) / sx)^2 +
                       ((y - cy_z) / sy)^2 +
                       ((z - cz) / sz)^2))
    
    if (j %% 2 == 0) {
      f <- -f
    }
    
    add_f(amp_gaussian * f, "focal")
  }
  
  # ------------------------------------------------------------
  # D. 10 3D composite spatial patterns
  # ------------------------------------------------------------
  
  # 91. smooth 3D trend plus periodic modulation
  add_f(
    amp_composite * (
      0.6 * (xc + yc + 0.8 * zc) +
        0.7 * sin(4 * pi * x + 1.5 * pi * z)
    ),
    "complex"
  )
  
  # 92. 3D hotspot contrast
  add_f(
    amp_composite * (
      gspot3d(0.30, 0.30, 0.30, 0.13, 0.13, 0.15) -
        gspot3d(0.70, 0.70, 0.70, 0.13, 0.13, 0.15)
    ),
    "complex"
  )
  
  # 93. 3D shell-shaped enrichment
  add_f(
    amp_composite * shell3d(0.50, 0.50, 0.50, r0 = 0.30, width = 0.06),
    "complex"
  )
  
  # 94. 3D sinusoidal interaction
  add_f(
    amp_composite * (
      sin(4 * pi * x) * cos(4 * pi * y) * sin(2 * pi * z)
    ),
    "complex"
  )
  
  # 95. smooth boundary with z-dependent bending and local enrichment
  add_f(
    amp_composite * (
      soft_step(x + y + 0.5 * z - 1.25, width = 0.07) +
        0.6 * gspot3d(0.70, 0.30, 0.65, 0.12, 0.12, 0.15)
    ),
    "complex"
  )
  
  # 96. diagonal 3D wave plus smooth trend
  add_f(
    amp_composite * (
      sin(2 * pi * (x + y + 0.7 * z)) +
        0.5 * xc -
        0.3 * zc
    ),
    "complex"
  )
  
  # 97. multiple 3D hotspots with different signs
  add_f(
    amp_composite * (
      gspot3d(0.25, 0.75, 0.30, 0.12, 0.12, 0.14) +
        gspot3d(0.75, 0.25, 0.70, 0.12, 0.12, 0.14) -
        0.7 * gspot3d(0.50, 0.50, 0.50, 0.16, 0.16, 0.16)
    ),
    "complex"
  )
  
  # 98. warped ring in xy with z-dependent phase
  add_f(
    amp_composite * (
      sin(5 * pi * r_xy + 2 * pi * z_warp2)
    ),
    "complex"
  )
  
  # 99. quadrant-wise effect changing along z
  add_f(
    amp_composite * (
      as.numeric(x > 0.5 & y > 0.5) -
        as.numeric(x < 0.5 & y < 0.5) +
        0.5 * sin(2 * pi * z)
    ),
    "complex"
  )
  
  # 100. radial oscillation with vertical trend
  add_f(
    amp_composite * (
      sin(6 * pi * r_xy + 1.5 * pi * z) +
        0.4 * zc
    ),
    "complex"
  )
  
  # 101. warped laminar oscillation along z
  add_f(
    amp_granularity * sin(4 * pi * z_warp1),
    "complex"
  )
  
  # 102. oblique 3D bands
  add_f(
    amp_granularity * sin(4 * pi * (0.45 * x + 0.35 * y + 0.55 * z)),
    "complex"
  )
  
  # 103. balanced 3D sinusoidal interaction
  add_f(
    amp_granularity * (sin(2 * pi * x) * cos(2 * pi * y) * sin(2 * pi * z)),
    "complex"
  )
  
  # 104. two opposite warped slabs
  d1 <- z - 0.35 - 0.10 * sin(2 * pi * x)
  d2 <- z - 0.65 - 0.10 * cos(2 * pi * y)
  add_f(
    amp_granularity * (sheet(d1, width = 0.08) - sheet(d2, width = 0.08)),
    "complex"
  )
  
  # 105. smooth boundary with z-dependent bending
  add_f(
    amp_granularity * soft_step(
      x + y + 0.5 * z - 1.25 + 0.12 * sin(2 * pi * z),
      width = 0.08
    ),
    "complex"
  )
  
  # 106. radial oscillation with z modulation
  add_f(
    amp_granularity * sin(4 * pi * r_xy + 2 * pi * z),
    "complex"
  )
  
  # 107. crossed smooth layers
  add_f(
    amp_granularity * (
      sin(2 * pi * (x + z_warp2)) -
        0.8 * sin(2 * pi * (y - z_warp2))
    ),
    "complex"
  )
  
  # 108. quadrant-wise effect changing along z
  add_f(
    amp_granularity * (
      soft_step(x - 0.5, width = 0.08) * sin(2 * pi * z) +
        soft_step(y - 0.5, width = 0.08) * cos(2 * pi * z)
    ),
    "complex"
  )
  
  # 109. wavy sheet with opposite nearby layer
  d <- z - 0.5 - 0.12 * sin(2 * pi * x) - 0.08 * cos(2 * pi * y)
  add_f(
    amp_granularity * (
      sheet(d, width = 0.06) -
        0.6 * sheet(d - 0.12, width = 0.06)
    ),
    "complex"
  )
  
  # 110. low-frequency 3D composite oscillation
  add_f(
    amp_granularity * (
      0.6 * sin(2 * pi * x) +
        0.6 * cos(2 * pi * y) -
        0.8 * sin(2 * pi * z_warp3)
    ),
    "complex"
  )
  
  # 111. helicoidal spatial variation
  add_f(
    amp_granularity * sin(theta_xy + 2 * pi * z),
    "complex"
  )
  
  # 112. 3D shell-like oscillation
  add_f(
    amp_granularity * sin(5 * pi * r_3d),
    "complex"
  )
  
  # 113. alternating 3D hotspots
  add_f(
    amp_granularity * (
      gspot3d(0.30, 0.30, 0.30, 0.13, 0.13, 0.13) -
        gspot3d(0.70, 0.70, 0.70, 0.13, 0.13, 0.13) +
        0.6 * gspot3d(0.30, 0.70, 0.55, 0.12, 0.12, 0.16)
    ),
    "complex"
  )
  
  # 114. saddle-like spatial surface with z oscillation
  add_f(
    amp_granularity * (
      2.5 * (xc^2 - yc^2) +
        0.5 * sin(2 * pi * z)
    ),
    "complex"
  )
  
  # 115. bent diagonal layer
  add_f(
    amp_granularity * soft_step(
      x - y + 0.25 * sin(2 * pi * z) + 0.10 * cos(2 * pi * x),
      width = 0.07
    ),
    "complex"
  )
  
  # 116. moving hotspot across z
  cx_z <- 0.5 + 0.22 * sin(2 * pi * z)
  cy_z <- 0.5 + 0.22 * cos(2 * pi * z)
  add_f(
    amp_granularity * exp(-0.5 * (((x - cx_z) / 0.12)^2 +
                                    ((y - cy_z) / 0.12)^2)),
    "complex"
  )
  
  # 117. opposite moving hotspots across z
  cx1 <- 0.35 + 0.15 * sin(2 * pi * z)
  cy1 <- 0.35 + 0.15 * cos(2 * pi * z)
  cx2 <- 0.65 - 0.15 * sin(2 * pi * z)
  cy2 <- 0.65 - 0.15 * cos(2 * pi * z)
  add_f(
    amp_granularity * (
      exp(-0.5 * (((x - cx1) / 0.12)^2 + ((y - cy1) / 0.12)^2)) -
        exp(-0.5 * (((x - cx2) / 0.12)^2 + ((y - cy2) / 0.12)^2))
    ),
    "complex"
  )
  
  # 118. alternating z layers with xy-dependent phase
  add_f(
    amp_granularity * sin(4 * pi * z + 2 * pi * x * y),
    "complex"
  )
  
  # 119. smooth checkerboard attenuated by z
  add_f(
    amp_granularity * (
      sin(3 * pi * x) * sin(3 * pi * y) * (0.5 + z)
    ),
    "complex"
  )
  
  # 120. warped radial rings with vertical trend
  add_f(
    amp_granularity * (
      sin(5 * pi * r_xy + 1.5 * pi * z) +
        0.3 * zc
    ),
    "complex"
  )
  
  stopifnot(length(f_list) == 120)
  
  # ============================================================
  # 3. gene-specific baselines
  # ============================================================
  gene_names <- paste0("gene", seq_len(G))
  mu_g <- rnorm(G, mean = mu_mean, sd = mu_sd)
  
  # ============================================================
  # 4. true spatial effects
  # ============================================================
  F_true <- matrix(0, nrow = n, ncol = G)
  colnames(F_true) <- gene_names
  
  if (n_svg > 0) {
    for (g in seq_len(n_svg)) {
      F_true[, g] <- f_list[[g]]
    }
  }
  
  # ============================================================
  # 5. generate continuous positive expression
  #    eta_true is log-scale mean.
  #    SVG and non-SVG genes use different noise levels.
  # ============================================================
  eta_true <- matrix(rep(mu_g, each = n), nrow = n, ncol = G)
  eta_true <- eta_true + F_true
  colnames(eta_true) <- gene_names
  
  eps <- matrix(0, nrow = n, ncol = G)
  
  if (n_svg > 0) {
    eps[, seq_len(n_svg)] <- matrix(
      rnorm(n * n_svg, mean = 0, sd = sigma_eps),
      nrow = n,
      ncol = n_svg
    )
  }
  
  if (G > n_svg) {
    non_svg_idx <- (n_svg + 1):G
    
    eps[, non_svg_idx] <- matrix(
      rnorm(n * length(non_svg_idx), mean = 0, sd = sigma_eps_non_svg),
      nrow = n,
      ncol = length(non_svg_idx)
    )
  }
  
  # optional truncation to prevent exp() from producing extremely large outliers
  if (!is.null(eps_cap) && is.finite(eps_cap) && eps_cap > 0) {
    eps <- pmax(pmin(eps, eps_cap), -eps_cap)
  }
  
  Y <- exp(eta_true + eps)
  colnames(Y) <- gene_names
  
  # ============================================================
  # 6. impose zeros at the end
  # ============================================================
  zero_mask <- matrix(FALSE, nrow = n, ncol = G)
  colnames(zero_mask) <- gene_names
  
  if (zero_rate_svg > 0 && n_svg > 0) {
    for (g in seq_len(n_svg)) {
      idx_zero <- sample.int(n, size = floor(zero_rate_svg * n), replace = FALSE)
      Y[idx_zero, g] <- 0
      zero_mask[idx_zero, g] <- TRUE
    }
  }
  
  if (zero_rate_non_svg > 0 && G > n_svg) {
    for (g in (n_svg + 1):G) {
      idx_zero <- sample.int(n, size = floor(zero_rate_non_svg * n), replace = FALSE)
      Y[idx_zero, g] <- 0
      zero_mask[idx_zero, g] <- TRUE
    }
  }
  
  # ============================================================
  # 7. truth labels
  # ============================================================
  is_svg <- c(rep(1L, n_svg), rep(0L, G - n_svg))
  
  pattern_label <- rep("non_svg", G)
  pattern_id <- rep(NA_integer_, G)
  
  if (n_svg > 0) {
    pattern_label[seq_len(n_svg)] <- f_type[seq_len(n_svg)]
    pattern_id[seq_len(n_svg)] <- seq_len(n_svg)
  }
  
  truth <- data.frame(
    gene = gene_names,
    is_svg = is_svg,
    pattern_id = pattern_id,
    pattern_label = pattern_label,
    mu_g = mu_g,
    sigma_eps = ifelse(is_svg == 1L, sigma_eps, sigma_eps_non_svg),
    stringsAsFactors = FALSE
  )
  
  list(
    Y = Y,
    S = S,
    is_svg = is_svg,
    svg_index = if (n_svg > 0) seq_len(n_svg) else integer(0),
    truth = truth,
    pattern_label = pattern_label,
    mu_g = mu_g,
    eta_true = eta_true,
    F_true = F_true,
    eps = eps,
    zero_mask = zero_mask,
    f_list = f_list,
    f_type = f_type
  )
}

# ============================================================
# Example workflow
# ============================================================

# ------------------------------------------------------------
# 1 Generate and save simulation data
# ------------------------------------------------------------
# seed = 1
# alpha = 0.01
sim <- simulate_svg_continuous_data(
  n = 10000,
  G = 1000,
  n_svg = 120,
  sigma_eps = 0.01,  # 0.01
  sigma_eps_non_svg = 0.01, # 0.01
  zero_rate = 0,
  mu_mean = 0,
  mu_sd = 1,
  seed = seed
)

saveRDS(sim, file = file.path(data_dir, paste0("simulation2_data_seed", seed, ".rds")))

print(table(sim$truth$pattern_label[sim$truth$is_svg == 1]))

# ============================================================
# 3. Read simulation data
#
# sim$Y: spots x genes
# sim$S: spots x 3
# ============================================================

sim <- readRDS(
  sim_rds_file
)

# SPARK expression:
# genes x spots

Y <- t(
  as.matrix(sim$Y)
)

# 3D coordinates:
# spots x 3

location <- as.data.frame(
  sim$S[, 1:3, drop = FALSE]
)

colnames(location) <- c(
  "x",
  "y",
  "z"
)


# ============================================================
# 4. Set gene and spot IDs
# ============================================================

gene_ids <- colnames(
  sim$Y
)

if (is.null(gene_ids)) {
  
  gene_ids <- paste0(
    "gene_",
    seq_len(nrow(Y))
  )
}

spot_ids <- paste0(
  "spot_",
  seq_len(ncol(Y))
)

rownames(Y) <- gene_ids
colnames(Y) <- spot_ids

rownames(location) <- spot_ids


cat(
  "SPARK input:",
  nrow(Y),
  "genes x",
  ncol(Y),
  "spots\n"
)


# ============================================================
# 5. Create SPARK object
#
# Disable SPARK-specific filtering
# ============================================================

spark <- CreateSPARKObject(
  counts = Y,
  location = location,
  percentage = 0,
  min_total_counts = 0
)

# sim$Y is already continuous expression.
# Use it directly for the Gaussian SPARK model.

spark@scaled_counts <- Y


# ============================================================
# 6. Fit Gaussian SPARK model
# ============================================================

t1 <- Sys.time()

spark <- spark.vc(
  spark,
  covariates = NULL,
  fit.model = "gaussian",
  num_core = 10,
  verbose = FALSE
)


# ============================================================
# 7. Run 3D spatial test
# ============================================================

spark <- spark.test.new(
  spark,
  check_positive = TRUE,
  verbose = FALSE
)

t2 <- Sys.time()

runtime <- t2 - t1

print(runtime)


# ============================================================
# 8. Extract SPARK results
# ============================================================

spark_table <- data.frame(
  gene = rownames(
    spark@res_mtest
  ),
  spark@res_mtest,
  row.names = NULL,
  check.names = FALSE
)

spark_pred <- as.integer(
  spark_table$adjusted_pvalue <
    p_cutoff
)

spark_svg_idx <- which(
  spark_pred == 1
)

spark_svg_genes <- spark_table$gene[
  spark_svg_idx
]


# ============================================================
# 9. Evaluate SVG detection
# ============================================================

truth <- as.integer(
  sim$is_svg
)

names(truth) <- gene_ids

truth <- truth[
  spark_table$gene
]

TP <- sum(
  truth == 1 &
    spark_pred == 1
)

FP <- sum(
  truth == 0 &
    spark_pred == 1
)

FN <- sum(
  truth == 1 &
    spark_pred == 0
)

TN <- sum(
  truth == 0 &
    spark_pred == 0
)

evaluation <- data.frame(
  TP = TP,
  FP = FP,
  FN = FN,
  TN = TN,
  Power = TP / (TP + FN),
  Precision = TP / (TP + FP),
  empirical_FDR = FP / (TP + FP)
)

print(
  evaluation
)

cat(
  "Number of SVGs detected by SPARK:",
  length(spark_svg_genes),
  "\n"
)


# ============================================================
# 10. Save results
# ============================================================

spark_table$pred <- spark_pred
spark_table$truth <- truth

spark_table <- spark_table[
  order(
    spark_table$adjusted_pvalue
  ),
  ,
  drop = FALSE
]

write.csv(
  spark_table,
  file = spark_csv_file,
  row.names = FALSE
)

saveRDS(
  list(
    spark = spark,
    result = spark_table,
    pvalues = spark_table$combined_pvalue,
    adjusted_pvalues = spark_table$adjusted_pvalue,
    pred = spark_pred,
    svg_idx = spark_svg_idx,
    svg_genes = spark_svg_genes,
    evaluation = evaluation,
    cutoff = p_cutoff,
    runtime = runtime,
    seed = seed
  ),
  file = spark_rds_file
)

cat(
  "SPARK results saved to:\n",
  spark_rds_file,
  "\n"
)
