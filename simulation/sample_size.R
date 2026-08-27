library(GeoSVG3D)

# ============================================================
# Sample-size sensitivity analysis
# ============================================================

data_dir <- "D:/github/GeoSVG3D/simulation/data"
rds_dir <- "D:/github/GeoSVG3D/simulation/output/sample_size"
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

sim_full <- readRDS(
  file.path(data_dir, "simulation2_data_seed1.rds")
)

sample_sizes <- c(1000, 2500, 5000, 7500, 10000)
seed <- 11
alpha <- 0.01

# Use one fixed random ordering so that smaller datasets are nested
set.seed(seed)
sample_order <- sample(seq_len(nrow(sim_full$Y)))

for (n_sub in sample_sizes) {
  
  idx_sub <- if (n_sub == nrow(sim_full$Y)) {
    seq_len(nrow(sim_full$Y))
  } else {
    sample_order[seq_len(n_sub)]
  }
  
  sim <- sim_full
  sim$Y <- sim_full$Y[idx_sub, , drop = FALSE]
  sim$S <- sim_full$S[idx_sub, , drop = FALSE]
  
  set.seed(seed)
  t1 <- Sys.time()
  
  fit_res <- GeoSVG3D(
    Y = sim$Y, S = sim$S,
    n_basis = 5, knn = 100,
    n_iter = 2000, burn_in = 2000,
    kappa_g = 0, alpha_g = 10, xi_g0 = 0,
    prior = "SpSL-L", xi_prop_sd = 1,
    tau2_g_update = FALSE, tau2_g_init = 1,
    n_cores = 10, seed = seed
  )
  
  runtime <- Sys.time() - t1
  
  p0g <- compute_p0g(fit_res, tol = 0)
  bfdr_res <- select_svg_bfdr(p0g, alpha = alpha)
  evaluation <- evaluate_svg_detection(sim$is_svg, bfdr_res$pred)
  
  saveRDS(
    list(
      sim = sim,
      sample_idx = idx_sub,
      sample_size = n_sub,
      fit_res = fit_res,
      p0g = p0g,
      bfdr_res = bfdr_res,
      evaluation = evaluation,
      runtime = runtime,
      seed = seed
    ),
    file.path(
      rds_dir,
      paste0(
        "Geosvg3d_simulation2_n", n_sub,
        "_results_seed", seed, ".rds"
      )
    )
  )
  
  rm(sim, fit_res, p0g, bfdr_res, evaluation)
  gc()
}