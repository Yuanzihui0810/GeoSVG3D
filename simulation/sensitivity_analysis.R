# ============================================================
# GeoSVG-3D hyperparameter sensitivity analysis
# ============================================================

library(GeoSVG3D)
library(RANN)
library(dplyr)


# ============================================================
# 0. Paths and fixed settings
# ============================================================

data_dir <- "D:/github/GeoSVG3D/simulation/data"
sensitivity_dir <- "D:/github/GeoSVG3D/simulation/sensitivity_output"

dir.create(sensitivity_dir, recursive = TRUE, showWarnings = FALSE)

simulation_seed <- 1
mcmc_seed <- 11
bfdr_alpha <- 0.01

default_K <- 5
default_knn <- 100
default_kappa <- 0
default_alpha <- 10

sim_file <- file.path(
  data_dir,
  paste0("simulation2_data_seed", simulation_seed, ".rds")
)

sim <- readRDS(sim_file)
truth <- sim$is_svg


# ============================================================
# 1. Default bandwidth
# ============================================================

S_scaled <- scale_spatial_coords_01(sim$S)

nn0 <- RANN::nn2(S_scaled, k = default_knn + 1)
h0 <- stats::median(nn0$nn.dists[, -1, drop = FALSE])


# ============================================================
# 2. Helper functions
# ============================================================

num_tag <- function(x) {
  gsub("\\.", "p", format(x, scientific = FALSE, trim = TRUE))
}


get_result_file <- function(K, knn, h_factor, kappa_g, alpha_g) {
  file_tag <- paste0(
    "K", K,
    "_knn", knn,
    "_hfactor", num_tag(h_factor),
    "_kappa", num_tag(kappa_g),
    "_alpha", num_tag(alpha_g)
  )

  file.path(
    sensitivity_dir,
    paste0("GeoSVG3D_sensitivity_", file_tag, ".rds")
  )
}


compute_auc <- function(truth, score) {
  n1 <- sum(truth == 1)
  n0 <- sum(truth == 0)
  r <- rank(score, ties.method = "average")

  (
    sum(r[truth == 1]) -
      n1 * (n1 + 1) / 2
  ) / (n1 * n0)
}


calculate_metrics <- function(truth, pred, p0g) {
  TP <- sum(pred == 1 & truth == 1)
  FP <- sum(pred == 1 & truth == 0)
  TN <- sum(pred == 0 & truth == 0)
  FN <- sum(pred == 0 & truth == 1)

  precision <- if (TP + FP > 0) TP / (TP + FP) else NA_real_
  power <- TP / (TP + FN)
  empirical_FDR <- if (TP + FP > 0) FP / (TP + FP) else 0
  specificity <- TN / (TN + FP)
  FPR <- FP / (FP + TN)
  accuracy <- (TP + TN) / length(truth)

  data.frame(
    AUC = compute_auc(truth, 1 - p0g),
    TP = TP,
    FP = FP,
    TN = TN,
    FN = FN,
    Precision = precision,
    Power = power,
    Empirical_FDR = empirical_FDR,
    Specificity = specificity,
    FPR = FPR,
    Accuracy = accuracy,
    Detected_SVG = TP + FP
  )
}


fit_and_evaluate <- function(
    n_basis, knn, h, kappa_g, alpha_g, tau2_g = 1
) {
  set.seed(mcmc_seed)
  t1 <- Sys.time()

  fit_res <- GeoSVG3D(
    Y = sim$Y, S = sim$S,
    n_basis = n_basis, knn = knn, h = h,
    n_iter = 2000, burn_in = 2000,
    kappa_g = kappa_g, alpha_g = alpha_g,
    xi_g0 = 0,
    prior = "SpSL-L",
    xi_prop_sd = 1,
    tau2_g_update = FALSE,
    tau2_g_init = tau2_g,
    n_cores = 10,
    seed = mcmc_seed
  )

  runtime <- Sys.time() - t1

  p0g <- compute_p0g(fit_res, tol = 0)
  bfdr_res <- select_svg_bfdr(p0g, alpha = bfdr_alpha)
  evaluation <- evaluate_svg_detection(truth, bfdr_res$pred)
  metrics <- calculate_metrics(truth, bfdr_res$pred, p0g)

  list(
    fit_res = fit_res,
    p0g = p0g,
    bfdr_res = bfdr_res,
    evaluation = evaluation,
    metrics = metrics,
    AUC = metrics$AUC,
    runtime = runtime
  )
}


round_metrics <- function(x) {
  x %>%
    mutate(
      h = round(h, 6),
      AUC = round(AUC, 4),
      Precision = round(Precision, 4),
      Power = round(Power, 4),
      Empirical_FDR = round(Empirical_FDR, 4),
      Specificity = round(Specificity, 4),
      FPR = round(FPR, 4),
      Accuracy = round(Accuracy, 4)
    )
}


# ============================================================
# 3. Main hyperparameter sensitivity settings
# ============================================================

knn_grid <- c(50, 75, 100, 125, 150)
K_grid <- c(3, 4, 5, 6, 7)
h_factor_grid <- c(0.8, 0.9, 1, 1.1, 1.2)
kappa_grid <- c(0, 0.05, 0.10, 0.20, 0.50)
alpha_grid <- c(6, 8, 10, 12, 14)

settings_knn <- data.frame(
  sensitivity = "knn",
  K = default_K,
  knn = knn_grid,
  h = NA_real_,
  h_factor = NA_real_,
  kappa_g = default_kappa,
  alpha_g = default_alpha
)

settings_K <- data.frame(
  sensitivity = "K",
  K = K_grid,
  knn = default_knn,
  h = h0,
  h_factor = 1,
  kappa_g = default_kappa,
  alpha_g = default_alpha
)

settings_h <- data.frame(
  sensitivity = "bandwidth",
  K = default_K,
  knn = default_knn,
  h = h_factor_grid * h0,
  h_factor = h_factor_grid,
  kappa_g = default_kappa,
  alpha_g = default_alpha
)

settings_spectral <- expand.grid(
  kappa_g = kappa_grid,
  alpha_g = alpha_grid
)

settings_spectral$sensitivity <- "spectral"
settings_spectral$K <- default_K
settings_spectral$knn <- default_knn
settings_spectral$h <- h0
settings_spectral$h_factor <- 1

settings_spectral <- settings_spectral[
  ,
  c(
    "sensitivity", "K", "knn", "h", "h_factor",
    "kappa_g", "alpha_g"
  )
]

settings <- rbind(
  settings_knn,
  settings_K,
  settings_h,
  settings_spectral
)

settings$run_id <- seq_len(nrow(settings))


# ============================================================
# 4. Run main hyperparameter sensitivity analysis
# ============================================================

all_results <- vector("list", nrow(settings))
metrics_list <- vector("list", nrow(settings))

for (i in seq_len(nrow(settings))) {
  s <- settings[i, ]

  h_now <- if (s$sensitivity == "knn") NULL else s$h

  outfile <- get_result_file(
    K = s$K,
    knn = s$knn,
    h_factor = s$h_factor,
    kappa_g = s$kappa_g,
    alpha_g = s$alpha_g
  )

  if (file.exists(outfile)) {
    result <- readRDS(outfile)
  } else {
    fit <- fit_and_evaluate(
      n_basis = s$K,
      knn = s$knn,
      h = h_now,
      kappa_g = s$kappa_g,
      alpha_g = s$alpha_g
    )

    result <- list(
      parameters = list(
        K = s$K,
        knn = s$knn,
        h = h_now,
        h0 = h0,
        h_factor = s$h_factor,
        kappa_g = s$kappa_g,
        alpha_g = s$alpha_g
      ),
      simulation_file = sim_file,
      simulation_seed = simulation_seed,
      mcmc_seed = mcmc_seed,
      bfdr_alpha = bfdr_alpha,
      fit_res = fit$fit_res,
      p0g = fit$p0g,
      bfdr_res = fit$bfdr_res,
      evaluation = fit$evaluation,
      AUC = fit$AUC,
      runtime = fit$runtime
    )

    saveRDS(result, outfile)
  }

  all_results[[i]] <- list(
    run_id = s$run_id,
    sensitivity = as.character(s$sensitivity),
    K = s$K,
    knn = s$knn,
    h = h_now,
    h_factor = s$h_factor,
    kappa_g = s$kappa_g,
    alpha_g = s$alpha_g,
    AUC = result$AUC,
    evaluation = result$evaluation,
    runtime = result$runtime,
    result_file = outfile
  )

  metrics_list[[i]] <- data.frame(
    Sensitivity = s$sensitivity,
    K = s$K,
    knn = s$knn,
    h = s$h,
    h_factor = s$h_factor,
    kappa_g = s$kappa_g,
    alpha_g = s$alpha_g,
    calculate_metrics(
      truth = truth,
      pred = result$bfdr_res$pred,
      p0g = result$p0g
    ),
    stringsAsFactors = FALSE
  )

  saveRDS(
    list(
      h0 = h0,
      settings = settings,
      completed_results = all_results[seq_len(i)],
      simulation_file = sim_file,
      simulation_seed = simulation_seed,
      mcmc_seed = mcmc_seed,
      bfdr_alpha = bfdr_alpha
    ),
    file = file.path(
      sensitivity_dir,
      "GeoSVG3D_sensitivity_progress.rds"
    )
  )

  rm(result)
  gc()
}


# ============================================================
# 5. Save main sensitivity summaries
# ============================================================

saveRDS(
  list(
    simulation_file = sim_file,
    simulation_seed = simulation_seed,
    mcmc_seed = mcmc_seed,
    bfdr_alpha = bfdr_alpha,
    default_parameters = list(
      K = default_K,
      knn = default_knn,
      h = h0,
      kappa_g = default_kappa,
      alpha_g = default_alpha
    ),
    sensitivity_grids = list(
      knn = knn_grid,
      K = K_grid,
      h_factor = h_factor_grid,
      kappa_g = kappa_grid,
      alpha_g = alpha_grid
    ),
    settings = settings,
    results = all_results
  ),
  file = file.path(
    sensitivity_dir,
    "GeoSVG3D_sensitivity_summary.rds"
  )
)

sensitivity_table <- bind_rows(metrics_list)
sensitivity_table_display <- round_metrics(sensitivity_table)

saveRDS(
  sensitivity_table,
  file = file.path(
    sensitivity_dir,
    "GeoSVG3D_sensitivity_metrics.rds"
  )
)

write.csv(
  sensitivity_table_display,
  file = file.path(
    sensitivity_dir,
    "GeoSVG3D_sensitivity_metrics.csv"
  ),
  row.names = FALSE
)


# ============================================================
# 6. tau^2 sensitivity analysis
# ============================================================

tau2_grid <- c(0.25, 0.5, 1, 2, 4)

settings_tau <- data.frame(
  sensitivity = "tau2",
  K = default_K,
  knn = default_knn,
  h = h0,
  h_factor = 1,
  kappa_g = default_kappa,
  alpha_g = default_alpha,
  tau2_g = tau2_grid
)

settings_tau$run_id <- seq_len(nrow(settings_tau))

tau_results <- vector("list", nrow(settings_tau))

for (i in seq_len(nrow(settings_tau))) {
  s <- settings_tau[i, ]

  outfile <- file.path(
    sensitivity_dir,
    paste0(
      "GeoSVG3D_sensitivity_tau2_",
      num_tag(s$tau2_g),
      ".rds"
    )
  )

  if (file.exists(outfile)) {
    result <- readRDS(outfile)
  } else {
    fit <- fit_and_evaluate(
      n_basis = default_K,
      knn = default_knn,
      h = h0,
      kappa_g = default_kappa,
      alpha_g = default_alpha,
      tau2_g = s$tau2_g
    )

    result <- list(
      sensitivity = "tau2",
      parameters = list(
        K = default_K,
        knn = default_knn,
        h = h0,
        h_factor = 1,
        kappa_g = default_kappa,
        alpha_g = default_alpha,
        tau2_g = s$tau2_g
      ),
      simulation_file = sim_file,
      simulation_seed = simulation_seed,
      mcmc_seed = mcmc_seed,
      bfdr_alpha = bfdr_alpha,
      fit_res = fit$fit_res,
      p0g = fit$p0g,
      bfdr_res = fit$bfdr_res,
      evaluation = fit$evaluation,
      metrics = fit$metrics,
      AUC = fit$AUC,
      runtime = fit$runtime
    )

    saveRDS(result, outfile)
  }

  tau_results[[i]] <- result
}


# ============================================================
# 7. Save tau^2 sensitivity summaries
# ============================================================

tau_table <- bind_rows(
  lapply(
    tau_results,
    function(x) {
      data.frame(
        Sensitivity = "tau2",
        tau2_g = x$parameters$tau2_g,
        K = x$parameters$K,
        knn = x$parameters$knn,
        h = x$parameters$h,
        kappa_g = x$parameters$kappa_g,
        alpha_g = x$parameters$alpha_g,
        x$metrics,
        stringsAsFactors = FALSE
      )
    }
  )
)

tau_table_display <- round_metrics(tau_table)

saveRDS(
  tau_table,
  file = file.path(
    sensitivity_dir,
    "GeoSVG3D_tau2_sensitivity_metrics.rds"
  )
)

write.csv(
  tau_table_display,
  file = file.path(
    sensitivity_dir,
    "GeoSVG3D_tau2_sensitivity_metrics.csv"
  ),
  row.names = FALSE
)

saveRDS(
  list(
    simulation_file = sim_file,
    simulation_seed = simulation_seed,
    mcmc_seed = mcmc_seed,
    bfdr_alpha = bfdr_alpha,
    tau2_grid = tau2_grid,
    settings = settings_tau,
    results = tau_results
  ),
  file = file.path(
    sensitivity_dir,
    "GeoSVG3D_tau2_sensitivity_summary.rds"
  )
)
