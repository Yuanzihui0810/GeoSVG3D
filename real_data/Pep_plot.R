library(data.table)
library(ggplot2)
library(scatterplot3d)
library(plot3D)
# ============================================================
# 0. Paths
# ============================================================

rawdata_dir <- "D:/MERFISH/brain/Deep_STARMAP"

out_dir <- "D:/MERFISH/brain/Deep_STARMAP/STAGATE_imputation_output/Peptidergic_neurons"
prefix <- "Peptidergic_neurons"

project_dir <- "D:/github/GeoSVG3D"
real_data_dir <- file.path(project_dir, "real_data")

geosvg_rds <- file.path(real_data_dir, "output",
                        paste0(prefix, "_GeoSVG3D_results.rds"))

scbsp_rds <- file.path(project_dir, "scBSP/Peptidergic_neurons", 
                       paste0("scBSP_results_", prefix, ".rds"))

fig_dir <- file.path(real_data_dir, "figure")
tab_dir <- file.path(real_data_dir, "table")
expr_3d_dir <- file.path(fig_dir, "representative_svg_3d_expression")
expr_2d_dir <- file.path(fig_dir, "representative_svg_xy_expression")

dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(expr_3d_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(expr_2d_dir, recursive = TRUE, showWarnings = FALSE)

main_blue <- "#1F4E79"
main_red  <- "#B22222"
main_gray <- "#666666"

# ============================================================
# 1. Load STAGATE-imputed Peptidergic neuron data
# ============================================================

Y <- readMM(file.path(out_dir, paste0(prefix, "_Y_STAGATE_ReX.mtx")))
spot_ids <- readLines(file.path(out_dir, paste0(prefix, "_spot_ids.tsv")))
gene_ids <- readLines(file.path(out_dir, paste0(prefix, "_gene_ids.tsv")))

cell_ids_stagate <- paste0("cell_", as.integer(spot_ids) + 1L)
rownames(Y) <- cell_ids_stagate
colnames(Y) <- gene_ids

S <- read.csv(
  file.path(out_dir, paste0(prefix, "_S_3d_coordinates.csv")),
  row.names = 1,
  check.names = FALSE
)

rownames(S) <- cell_ids_stagate

stopifnot(nrow(Y) == nrow(S))

S <- as.matrix(S[, 1:3])
storage.mode(S) <- "double"

S_scaled <- scale_spatial_coords_01(S)

cat("Loaded Y:", dim(Y), "\n")
cat("Loaded S:", dim(S_scaled), "\n")

# ============================================================
# 1.2. Load raw count data
# ============================================================

expr_file    <- file.path(rawdata_dir, "Brain_Deep_STARmap_expression_matrix.csv")
spatial_file <- file.path(rawdata_dir, "Brain_Deep_STARmap_spatial.csv")

Raw_Exp <- fread(expr_file, data.table = FALSE)      # cells x genes
Spatial_raw <- fread(spatial_file, data.table = FALSE)

cell_id_all <- paste0("cell_", seq_len(nrow(Raw_Exp)))
rownames(Raw_Exp) <- cell_id_all
rownames(Spatial_raw) <- cell_id_all

# subset the same cell type
idx_pep <- which(Spatial_raw$FUSEmap_main_level == "Peptidergic neurons")
Raw_Exp_pep <- Raw_Exp[idx_pep, , drop = FALSE]

# align to the exact same cells and genes used in STAGATE/GeoSVG-3D
Y_raw <- Raw_Exp_pep[rownames(Y), colnames(Y), drop = FALSE]
Y_raw <- as.matrix(Y_raw)
storage.mode(Y_raw) <- "double"

stopifnot(all(rownames(Y_raw) == rownames(Y)))
stopifnot(all(colnames(Y_raw) == colnames(Y)))

cat("Raw count matrix dimension:", dim(Y_raw), "\n")

# ============================================================
# 2. Load GeoSVG-3D result
# ============================================================

if (!file.exists(geosvg_rds)) {
  stop("GeoSVG-3D RDS file not found: ", geosvg_rds)
}

geo_obj <- readRDS(geosvg_rds)

p0g <- geo_obj$p0g
bfdr_res <- geo_obj$bfdr_res

if (is.null(names(p0g))) {
  names(p0g) <- colnames(Y)
}

if (!is.null(geo_obj$svg_tab)) {
  geosvg_tab <- geo_obj$svg_tab
} else {
  geosvg_tab <- data.frame(
    gene = colnames(Y),
    p0g = as.numeric(p0g[colnames(Y)]),
    pred = as.integer(bfdr_res$pred),
    stringsAsFactors = FALSE
  )
}

colnames(geosvg_tab)[colnames(geosvg_tab) == "pred"] <- "geosvg_pred"

geosvg_tab <- geosvg_tab[, c("gene", "p0g", "geosvg_pred")]
geosvg_tab$gene <- as.character(geosvg_tab$gene)
geosvg_tab$p0g <- as.numeric(geosvg_tab$p0g)
geosvg_tab$geosvg_pred <- as.integer(geosvg_tab$geosvg_pred)

geosvg_tab <- geosvg_tab[order(geosvg_tab$p0g), ]

cat("Number of GeoSVG-3D SVGs:", sum(geosvg_tab$geosvg_pred == 1), "\n")

# ============================================================
# 3. Load scBSP result
# ============================================================

if (!file.exists(scbsp_rds)) {
  stop("scBSP RDS file not found: ", scbsp_rds)
}

scbsp_obj <- readRDS(scbsp_rds)
P_values <- as.data.frame(scbsp_obj$P_values)

colnames(P_values)[1:2] <- c("gene", "scbsp_pvalue")

scbsp_tab <- P_values[, c("gene", "scbsp_pvalue")]
scbsp_tab$gene <- as.character(scbsp_tab$gene)
scbsp_tab$scbsp_pvalue <- as.numeric(scbsp_tab$scbsp_pvalue)
scbsp_tab$scbsp_pred <- as.integer(scbsp_tab$scbsp_pvalue < 0.05)

scbsp_tab <- scbsp_tab[order(scbsp_tab$scbsp_pvalue), ]

cat("Number of scBSP SVGs:", sum(scbsp_tab$scbsp_pred == 1), "\n")

# ============================================================
# 4. Merge GeoSVG-3D and scBSP by gene name
# ============================================================

# Gene-level zero proportion in raw count data
gene_zero_prop_raw <- colMeans(Y_raw == 0, na.rm = TRUE)

zero_prop_df <- data.frame(
  gene = colnames(Y_raw),
  raw_zero_prop = as.numeric(gene_zero_prop_raw),
  stringsAsFactors = FALSE
)

head(zero_prop_df)
summary(zero_prop_df$raw_zero_prop)

merged_res <- merge(
  geosvg_tab,
  scbsp_tab,
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)

# Then merge raw zero proportion
merged_res <- merge(
  merged_res,
  zero_prop_df,
  by = "gene",
  all.x = TRUE,
  sort = FALSE
)

# Genes removed by SpFilter are treated as not selected by scBSP.
merged_res$scbsp_pred[is.na(merged_res$scbsp_pred)] <- 0L

merged_res$category <- with(
  merged_res,
  ifelse(geosvg_pred == 1 & scbsp_pred == 1, "Both",
         ifelse(geosvg_pred == 1 & scbsp_pred == 0, "GeoSVG-3D only",
                ifelse(geosvg_pred == 0 & scbsp_pred == 1, "scBSP only", "Neither")))
)

merged_res <- merged_res[order(merged_res$p0g), ]

write.csv(
  merged_res,
  file.path(tab_dir, paste0(prefix, "_GeoSVG3D_scBSP_merged_results.csv")),
  row.names = FALSE
)

# ============================================================
# 5. Summary tables
# ============================================================

method_summary <- data.frame(
  method = c("GeoSVG-3D", "scBSP"),
  n_detected_svg = c(
    sum(merged_res$geosvg_pred == 1),
    sum(merged_res$scbsp_pred == 1)
  )
)

overlap_summary <- data.frame(
  category = c("Both", "GeoSVG-3D only", "scBSP only", "Neither"),
  n_genes = c(
    sum(merged_res$category == "Both"),
    sum(merged_res$category == "GeoSVG-3D only"),
    sum(merged_res$category == "scBSP only"),
    sum(merged_res$category == "Neither")
  )
)

write.csv(
  method_summary,
  file.path(tab_dir, paste0(prefix, "_method_detected_svg_summary.csv")),
  row.names = FALSE
)

write.csv(
  overlap_summary,
  file.path(tab_dir, paste0(prefix, "_GeoSVG3D_scBSP_overlap_summary.csv")),
  row.names = FALSE
)

print(method_summary)
print(overlap_summary)

# ============================================================
# 6. Plot p0g distribution
# ============================================================

plot_p0g_df <- merged_res
plot_p0g_df$geosvg_status <- ifelse(
  plot_p0g_df$geosvg_pred == 1,
  "Selected",
  "Not selected"
)

p_p0g <- ggplot(plot_p0g_df, aes(x = p0g, fill = geosvg_status)) +
  geom_histogram(bins = 50, color = "white", linewidth = 0.2, alpha = 0.9) +
  scale_fill_manual(values = c("Selected" = main_blue, "Not selected" = "grey75")) +
  theme_bw(base_size = 13) +
  labs(
    x = expression(hat(p)[0*g]),
    y = "Number of genes",
    fill = "GeoSVG-3D"
  ) +
  theme(panel.grid.minor = element_blank())

ggsave(
  file.path(fig_dir, paste0(prefix, "_GeoSVG3D_p0g_distribution.pdf")),
  p_p0g,
  width = 5.5,
  height = 4
)

# ============================================================
# 7. BFDR selection curve
# ============================================================

p0g_vec <- merged_res$p0g
ord <- order(p0g_vec)
p_sorted <- p0g_vec[ord]
cum_bfdr <- cumsum(p_sorted) / seq_along(p_sorted)
k_selected <- sum(merged_res$geosvg_pred == 1)

bfdr_curve_df <- data.frame(
  rank = seq_along(p_sorted),
  p0g_sorted = p_sorted,
  cum_bfdr = cum_bfdr,
  selected = seq_along(p_sorted) <= k_selected
)

p_bfdr <- ggplot(bfdr_curve_df, aes(x = rank, y = cum_bfdr)) +
  geom_line(color = main_blue, linewidth = 0.8) +
  geom_hline(yintercept = 0.05, linetype = 2, color = main_red, linewidth = 0.7) +
  geom_vline(xintercept = k_selected, linetype = 2, color = main_gray, linewidth = 0.7) +
  theme_bw(base_size = 13) +
  labs(
    x = "Gene rank by posterior probability",
    y = "Estimated Bayesian FDR"
  ) +
  theme(panel.grid.minor = element_blank())

ggsave(
  file.path(fig_dir, paste0(prefix, "_GeoSVG3D_BFDR_curve.pdf")),
  p_bfdr,
  width = 5.5,
  height = 4
)

# ============================================================
# 8. Helper functions for 3D expression plots
# ============================================================
make_expr_color_info_imputed <- function(
    v,
    zero_col = "#6A51A3",
    positive_palette = c(
      "#2C7BB6",
      "#00A6CA",
      "#00CC66",
      "#FDE725",
      "#F46D43",
      "#D73027"
    ),
    n_col = 100,
    use_log1p = FALSE,
    upper_q = 0.99,
    force_zero_min = TRUE,
    point_alpha = 0.8
) {
  v_raw <- as.numeric(v)
  
  if (use_log1p && any(v_raw <= -1, na.rm = TRUE)) {
    stop("use_log1p = TRUE requires all expression values to be greater than -1.")
  }
  
  raw_min_data <- min(v_raw, na.rm = TRUE)
  
  # For non-negative imputed expression, use 0 as the lower bound.
  raw_min <- if (force_zero_min && raw_min_data >= 0) {
    0
  } else {
    raw_min_data
  }
  
  raw_max <- quantile(v_raw, upper_q, na.rm = TRUE)
  
  if (!is.finite(raw_max) || raw_max <= raw_min) {
    raw_max <- max(v_raw, na.rm = TRUE)
  }
  
  if (!is.finite(raw_max) || raw_max <= raw_min) {
    raw_max <- raw_min + 1e-8
  }
  
  v_clip <- pmin(pmax(v_raw, raw_min), raw_max)
  
  map_v <- if (use_log1p) log1p(v_clip) else v_clip
  map_min <- if (use_log1p) log1p(raw_min) else raw_min
  map_max <- if (use_log1p) log1p(raw_max) else raw_max
  
  pal <- colorRampPalette(c(zero_col, positive_palette))(n_col)
  pal_points <- grDevices::adjustcolor(pal, alpha.f = point_alpha)
  
  col_idx <- cut(
    map_v,
    breaks = seq(map_min, map_max, length.out = n_col + 1),
    include.lowest = TRUE,
    labels = FALSE
  )
  
  col_idx[is.na(col_idx)] <- 1L
  cols <- pal_points[col_idx]
  
  list(
    raw_min = raw_min,
    raw_max = raw_max,
    map_v = map_v,
    map_min = map_min,
    map_max = map_max,
    cols = cols,
    pal = pal,
    pal_points = pal_points
  )
}

draw_simple_colorbar <- function(
    pal,
    raw_min,
    raw_max,
    x_left = 0.04,
    x_right = 0.20,
    y0 = 0.10,
    y1 = 0.90,
    n_tick = 3,
    cex_tick = 0.8
) {
  n_col <- length(pal)
  
  y_bar <- seq(y0, y1, length.out = n_col + 1)
  
  for (i in seq_len(n_col)) {
    rect(
      x_left, y_bar[i],
      x_right, y_bar[i + 1],
      col = pal[i],
      border = NA
    )
  }
  
  tick_values <- seq(raw_min, raw_max, length.out = n_tick)
  y_ticks <- seq(y0, y1, length.out = n_tick)
  
  segments(
    x_right,
    y_ticks,
    x_right + 0.035,
    y_ticks,
    lwd = 0.8
  )
  
  text(
    x_right + 0.055,
    y_ticks,
    labels = formatC(tick_values, format = "fg", digits = 2),
    adj = c(0, 0.5),
    cex = cex_tick
  )
}

plot_one_gene_imputed_3d <- function(
    gene,
    Y_imp,
    S_scaled,
    out_file,
    point_cex = 0.55,
    theta = 40,
    phi = 25,
    use_log1p_for_color = FALSE,
    upper_q = 0.99,
    width = 6.2,
    height = 5.2,
    dpi = 600,
    order_by_expr = TRUE,
    point_alpha = 0.8,
    n_tick = 3
) {
  if (!requireNamespace("plot3D", quietly = TRUE)) {
    stop("Please install plot3D first.")
  }
  
  expr <- as.numeric(Y_imp[, gene])
  
  info <- make_expr_color_info_imputed(
    v = expr,
    use_log1p = use_log1p_for_color,
    upper_q = upper_q,
    point_alpha = point_alpha
  )
  
  idx <- seq_along(expr)
  if (order_by_expr) {
    idx <- order(info$map_v, decreasing = FALSE)
  }
  
  title_text <- paste0(
    toupper(substr(gene, 1, 1)),
    tolower(substr(gene, 2, nchar(gene)))
  )
  
  # ---------------- output device ----------------
  ext <- tolower(tools::file_ext(out_file))
  
  if (ext == "pdf") {
    pdf(out_file, width = width, height = height, useDingbats = FALSE)
  } else if (ext == "png") {
    png(
      filename = out_file,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      bg = "white"
    )
  } else {
    stop("out_file must end with .pdf or .png")
  }
  
  on.exit(dev.off())
  
  layout(matrix(c(1, 2), nrow = 1), widths = c(4.5, 1.8))
  
  # ---------------- main 3D plot ----------------
  par(mar = c(2.5, 2.5, 2.8, 0.5))
  
  plot3D::scatter3D(
    x = S_scaled[idx, 1],
    y = S_scaled[idx, 2],
    z = S_scaled[idx, 3],
    colvar = info$map_v[idx],
    col = info$pal_points,
    clim = c(info$map_min, info$map_max),
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
  
  # ---------------- legend panel ----------------
  par(mar = c(3, 0.2, 3, 4))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  
  draw_simple_colorbar(
    pal = info$pal,
    raw_min = info$raw_min,
    raw_max = info$raw_max,
    x_left = 0.02,
    x_right = 0.18,
    y0 = 0.08,
    y1 = 0.92,
    n_tick = n_tick,
    cex_tick = 0.85
  )
}


plot_one_gene_imputed_xy <- function(
    gene,
    Y_imp,
    S_scaled,
    out_file,
    point_cex = 0.45,
    use_log1p_for_color = FALSE,
    upper_q = 0.99,
    width = 5.5,
    height = 5,
    dpi = 600,
    order_by_expr = TRUE
) {
  expr <- as.numeric(Y_imp[, gene])
  
  # -------- palette --------
  zero_col <- "#6A51A3"
  positive_palette <- c(
    "#2C7BB6",  # blue
    "#00A6CA",  # cyan
    "#00CC66",  # green
    "#FDE725",  # yellow
    "#F46D43",  # orange
    "#D73027"   # red
  )
  
  n_col <- 100
  
  # continuous palette
  pal <- colorRampPalette(c(zero_col, positive_palette))(n_col)
  
  v <- expr
  if (use_log1p_for_color) {
    v <- log1p(v)
  }
  
  v_max <- as.numeric(quantile(v, upper_q, na.rm = TRUE))
  v_min <- min(v, na.rm = TRUE)
  
  if (!is.finite(v_max) || v_max <= v_min) {
    v_max <- max(v, na.rm = TRUE)
  }
  if (!is.finite(v_max) || v_max <= v_min) {
    v_max <- v_min + 1e-8
  }
  
  v_plot <- pmin(v, v_max)
  
  col_idx <- cut(
    v_plot,
    breaks = seq(v_min, v_max, length.out = n_col + 1),
    include.lowest = TRUE,
    labels = FALSE
  )
  
  cols <- pal[col_idx]
  
  # Put high-expression points on top
  idx <- seq_along(expr)
  if (order_by_expr) {
    idx <- order(v_plot, decreasing = FALSE)
  }
  
  x <- S_scaled[idx, 1]
  y <- S_scaled[idx, 2]
  cols <- cols[idx]
  
  title_text <- paste0(
    toupper(substr(gene, 1, 1)),
    tolower(substr(gene, 2, nchar(gene)))
  )
  
  # ---------------- output device ----------------
  ext <- tolower(tools::file_ext(out_file))
  
  if (ext == "pdf") {
    pdf(out_file, width = width, height = height, useDingbats = FALSE)
  } else if (ext == "png") {
    png(
      filename = out_file,
      width = width,
      height = height,
      units = "in",
      res = dpi,
      bg = "white"
    )
  } else {
    stop("out_file must end with .pdf or .png")
  }
  
  on.exit(dev.off())
  
  layout(matrix(c(1, 2), nrow = 1), widths = c(4.5, 1.4))
  
  # ---------------- xy plot ----------------
  par(mar = c(4, 4, 3, 0.5))
  
  plot(
    x, y,
    col = cols,
    pch = 16,
    cex = point_cex,
    asp = 1,
    xlab = "s1",
    ylab = "s2",
    main = title_text,
    axes = TRUE
  )
  
  grid(col = "grey85", lty = 2)
  box()
  
  # ---------------- legend ----------------
  par(mar = c(4, 0.2, 3, 4))
  plot.new()
  plot.window(xlim = c(0, 1), ylim = c(0, 1))
  
  x_left  <- 0.04
  x_right <- 0.20
  y0 <- 0.10
  y1 <- 0.90
  
  # continuous color bar, no title text
  y_bar <- seq(y0, y1, length.out = n_col + 1)
  
  for (i in seq_len(n_col)) {
    rect(
      x_left, y_bar[i],
      x_right, y_bar[i + 1],
      col = pal[i],
      border = NA
    )
  }
  
  # evenly spaced ticks, no ">="
  # you can change n_tick to 3 / 4 / 5 if you want
  n_tick <- 3
  tick_values <- seq(v_min, v_max, length.out = n_tick)
  y_ticks <- seq(y0, y1, length.out = n_tick)
  
  segments(x_right, y_ticks, x_right + 0.035, y_ticks)
  
  text(
    x_right + 0.055,
    y_ticks,
    labels = formatC(tick_values, format = "fg", digits = 2),
    adj = c(0, 0.5),
    cex = 0.8
  )
}

# ============================================================
# 9. Select representative genes
# ============================================================

rep_genes <- c(
  "ADCYAP1R1",
  "ADORA1",
  "ALCAM",
  "CALB2",
  "CBLN4",
  "CLU",
  "DKK3",
  "GPX3", 
  "HAP1"
)

for (g in rep_genes) {
  plot_one_gene_imputed_3d(
    gene = g,
    Y_imp = Y,
    S_scaled = S_scaled,
    out_file = file.path(expr_3d_dir, paste0(g, "_STAGATE_Y_3d.png")),
    point_cex = 0.5,
    theta = 300,
    phi = 35,
    use_log1p_for_color = FALSE,
    upper_q = 1,
    width = 6.5,
    height = 5,
    point_alpha = 0.8
  )
}

for (g in rep_genes) {
  plot_one_gene_imputed_xy(
    gene = g,
    Y_imp = Y,   
    S_scaled = S_scaled,
    out_file = file.path(expr_2d_dir, paste0(g, "_imputed_xy.png")),
    point_cex = 0.45,
    use_log1p_for_color = FALSE,
    upper_q = 1
  )
}