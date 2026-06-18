library(plot3D)
library(pROC)

# ============================================================
# 1. Read saved results
# ============================================================

project_dir <- "D:/github/GL-neuronized-svg"
seed <- 1

geosvg_res <- readRDS(file.path(
  project_dir,
  "simulation/output",
  paste0("Geosvg3d_simulation_results_seed", seed, ".rds")
))

scbsp_res <- readRDS(file.path(
  project_dir,
  "scBSP",
  paste0("scbsp_simulation_results_seed", seed, ".rds")
))

sim <- geosvg_res$sim

truth <- sim$is_svg
pred_geosvg <- geosvg_res$bfdr_res$pred
pred_scbsp <- scbsp_res$pred

score_geosvg <- 1 - geosvg_res$p0g
score_scbsp <- 1 - scbsp_res$pvalues

fig_dir <- file.path(project_dir, "simulation/figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 2. Representative 3D expression patterns
# ============================================================


plot_one_simulated_Y_pattern_3d_plot3D <- function(
    sim,
    gene,
    output_file = NULL,
    point_cex = 0.5,
    max_points = 10000,
    theta = 40,
    phi = 20,
    seed = 1,
    width = 6.5,
    height = 5.5,
    dpi = 300,
    color_transform = c("log1p", "none")
) {
  color_transform <- match.arg(color_transform)
  
  if (!requireNamespace("plot3D", quietly = TRUE)) {
    stop("Please install plot3D first.")
  }
  
  S <- sim$S
  Y <- sim$Y
  
  stopifnot(ncol(S) >= 3)
  
  if (is.character(gene)) {
    gene_idx <- match(gene, colnames(Y))
    if (is.na(gene_idx)) stop("gene name not found.")
  } else {
    gene_idx <- gene
  }
  
  stopifnot(gene_idx >= 1, gene_idx <= ncol(Y))
  
  expr <- Y[, gene_idx]
  if (color_transform == "log1p") {
    expr <- log1p(expr)
  }
  
  set.seed(seed)
  idx <- seq_len(nrow(S))
  if (!is.null(max_points) && length(idx) > max_points) {
    idx <- sample(idx, max_points)
  }
  
  x <- S[idx, 1]
  y <- S[idx, 2]
  z <- S[idx, 3]
  val <- expr[idx]
  
  val_range <- range(val, na.rm = TRUE)
  if (!all(is.finite(val_range))) {
    stop("Expression values contain no finite values.")
  }
  
  if (diff(val_range) == 0) {
    clim_use <- val_range + c(-0.5, 0.5)
  } else {
    clim_use <- val_range
  }
  
  pal <- hcl.colors(100, "viridis")
  
  gene_name <- colnames(Y)[gene_idx]
  pattern_name <- sim$truth$pattern_label[gene_idx] #format_pattern_name(sim$truth$pattern_label[gene_idx])
  
  # New title format: gene1 (focal)
  title_text <- paste0(gene_name, " (", pattern_name, ")")
  
  if (!is.null(output_file)) {
    dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)
    ext <- tolower(tools::file_ext(output_file))
    
    if (ext == "pdf") {
      pdf(output_file, width = width, height = height)
    } else if (ext == "png") {
      png(output_file, width = width, height = height, units = "in", res = dpi, bg = "white")
    } else {
      stop("output_file must end with .pdf or .png")
    }
    
    on.exit(dev.off(), add = TRUE)
  }
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  
  layout(matrix(c(1, 2), nrow = 1), widths = c(4.5, 1.3))
  
  par(mar = c(2.5, 2.5, 2.8, 0.5),
      mgp = c(5, 0.8, 0)
  )
  plot3D::scatter3D(
    x = x,
    y = y,
    z = z,
    colvar = val,
    col = pal,
    clim = clim_use,
    pch = 16,
    cex = point_cex,
    theta = theta,
    phi = phi,
    bty = "g",
    ticktype = "detailed",
    xlab = "s1",
    ylab = "s2",
    zlab = "s3",
    main = title_text,
    colkey = FALSE,
    lighting = FALSE
  )
  
  # Manual color bar without the "Expression" title
  par(mar = c(3, 0.5, 2, 3))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  
  x_left <- 0.08
  x_right <- 0.22
  y0 <- 0.08
  y1 <- 0.95
  
  y_bar <- seq(y0, y1, length.out = length(pal) + 1)
  for (i in seq_along(pal)) {
    rect(x_left, y_bar[i], x_right, y_bar[i + 1], col = pal[i], border = NA)
  }
  
  ticks <- pretty(val, n = 5)
  ticks <- ticks[ticks >= val_range[1] & ticks <= val_range[2]]
  
  if (length(ticks) == 0) {
    ticks <- mean(val_range)
  }
  
  if (diff(val_range) == 0) {
    y_ticks <- rep((y0 + y1) / 2, length(ticks))
  } else {
    y_ticks <- y0 + (ticks - val_range[1]) / diff(val_range) * (y1 - y0)
  }
  
  segments(x_right, y_ticks, x_right + 0.03, y_ticks)
  text(
    x_right + 0.05,
    y_ticks,
    labels = formatC(ticks, digits = 3, format = "fg"),
    adj = c(0, 0.5),
    cex = 0.8
  )
  
  invisible(NULL)
}


plot_simulated_Y_patterns_3d_each_plot3D <- function(
    sim,
    genes = c(1, 31, 61, 91, 101, 111),
    output_dir,
    file_type = "pdf",
    point_cex = 0.5,
    theta = 305,
    phi = 40,
    color_transform = "log1p"
) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  for (gene in genes) {
    output_file <- file.path(
      output_dir,
      paste0(
        sprintf("gene%04d", gene), "_",
        sim$truth$pattern_label[gene], "_3d.",
        file_type
      )
    )
    
    plot_one_simulated_Y_pattern_3d_plot3D(
      sim = sim,
      gene = gene,
      output_file = output_file,
      point_cex = point_cex,
      theta = theta,
      phi = phi,
      color_transform = color_transform
    )
  }
}


plot_simulated_Y_patterns_3d_each_plot3D(
  sim,
  genes = c(1, 31, 61, 91, 101, 111),
  output_dir = file.path(fig_dir, "plot3D"),
  file_type = "pdf",
  point_cex = 0.5,
  theta = 305,
  phi = 40,
  color_transform = "log1p"
)

# ============================================================
# 3. Detection performance table
# ============================================================

get_metrics <- function(method, pred) {
  TP <- sum(truth == 1 & pred == 1)
  FP <- sum(truth == 0 & pred == 1)
  FN <- sum(truth == 1 & pred == 0)
  TN <- sum(truth == 0 & pred == 0)
  
  data.frame(
    method = method,
    TP = TP,
    FP = FP,
    FN = FN,
    TN = TN,
    power = TP / (TP + FN),
    precision = TP / (TP + FP)
  )
}

performance_table <- rbind(
  get_metrics("GeoSVG-3D", pred_geosvg),
  get_metrics("scBSP", pred_scbsp)
)

write.csv(
  performance_table,
  file.path(fig_dir, "sim_performance.csv"),
  row.names = FALSE
)

print(performance_table)


# ============================================================
# 4. ROC curves and AUC
# ============================================================

roc_geosvg <- roc(
  response = truth,
  predictor = score_geosvg,
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

roc_scbsp <- roc(
  response = truth,
  predictor = score_scbsp,
  levels = c(0, 1),
  direction = "<",
  quiet = TRUE
)

roc_geosvg_df <- data.frame(
  FPR = 1 - roc_geosvg$specificities,
  TPR = roc_geosvg$sensitivities
)
roc_geosvg_df <- roc_geosvg_df[
  order(roc_geosvg_df$FPR, roc_geosvg_df$TPR),
]

roc_scbsp_df <- data.frame(
  FPR = 1 - roc_scbsp$specificities,
  TPR = roc_scbsp$sensitivities
)
roc_scbsp_df <- roc_scbsp_df[
  order(roc_scbsp_df$FPR, roc_scbsp_df$TPR),
]

auc_geosvg <- as.numeric(auc(roc_geosvg))
auc_scbsp <- as.numeric(auc(roc_scbsp))

pdf(
  file.path(fig_dir, "sim_roc_auc.pdf"),
  width = 5.5,
  height = 5.5
)

par(mar = c(5, 5, 3, 2))

plot(
  roc_geosvg_df$FPR,
  roc_geosvg_df$TPR,
  type = "s",
  xlim = c(0, 0.5),
  ylim = c(0, 1),
  lwd = 2.5,
  col = "red",
  xlab = "False positive rate",
  ylab = "True positive rate",
  main = "ROC curve"
)

lines(
  roc_scbsp_df$FPR,
  roc_scbsp_df$TPR,
  type = "s",
  lwd = 2.5,
  col = "blue"
)

legend(
  "bottomright",
  legend = c(
    paste0("GeoSVG-3D, AUC = ", sprintf("%.4f", auc_geosvg)),
    paste0("scBSP, AUC = ", sprintf("%.4f", auc_scbsp))
  ),
  col = c("red", "blue"),
  lwd = 2.5,
  bty = "n"
)

dev.off()

auc_table <- data.frame(
  method = c("GeoSVG-3D", "scBSP"),
  AUC = c(auc_geosvg, auc_scbsp)
)

write.csv(
  auc_table,
  file.path(fig_dir, "sim_auc.csv"),
  row.names = FALSE
)

print(auc_table)
