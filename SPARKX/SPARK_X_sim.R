library(SPARK)
library(Matrix)
library(GeoSVG3D)
# ============================================================
# Settings
# ============================================================

seed <- 1
p_cutoff <- 0.01

data_dir <- "D:/github/GeoSVG3D/simulation/data"
out_dir <- "D:/github/GeoSVG3D/SPARKX/output"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

sim_rds_file <- file.path(
  data_dir,
  paste0("simulation2_data_seed", seed, ".rds")
)

sparkx_rds_file <- file.path(
  out_dir,
  paste0("sparkx_simulation2_results.rds")
)

sparkx_csv_file <- file.path(
  out_dir,
  paste0("sparkx_simulation2_table.csv")
)


# ============================================================
# Read simulation data
# ============================================================

sim <- readRDS(sim_rds_file)

# coordinates: spots x 3
Coords <- as.matrix(
  sim$S[, 1:3, drop = FALSE]
)

# expression: genes x spots
Exp <- Matrix(
  t(sim$Y),
  sparse = TRUE
)

gene_names <- colnames(sim$Y)

if (is.null(gene_names)) {
  gene_names <- paste0("gene_", seq_len(ncol(sim$Y)))
}

rownames(Exp) <- gene_names

spot_names <- rownames(sim$Y)

if (is.null(spot_names)) {
  spot_names <- paste0("spot_", seq_len(nrow(sim$Y)))
}

colnames(Exp) <- spot_names
rownames(Coords) <- spot_names
colnames(Coords) <- c("x", "y", "z")


# ============================================================
# Run SPARK-X
# ============================================================

set.seed(500)
t1 <- Sys.time()

sparkx_fit <- sparkx(
  count_in = Exp,
  locus_in = Coords,
  numCores = 1,
  option = "mixture",
  verbose = FALSE
)

t2 <- Sys.time()

runtime <- t2 - t1

print(runtime)

print("Finish")
# ============================================================
# Extract results
# ============================================================

sparkx_tab <- data.frame(
  gene = rownames(sparkx_fit$res_mtest),
  sparkx_fit$res_mtest,
  row.names = NULL,
  check.names = FALSE
)

sparkx_tab$pred <- as.integer(
  sparkx_tab$adjustedPval < p_cutoff
)

truth <- sim$is_svg
names(truth) <- gene_names

truth_used <- truth[sparkx_tab$gene]

res_sparkx <- evaluate_svg_detection(
  truth = truth_used,
  pred = sparkx_tab$pred
)

sparkx_tab$truth <- truth_used

sparkx_svg_idx <- which(sparkx_tab$pred == 1)
sparkx_svg_genes <- sparkx_tab$gene[sparkx_svg_idx]

sparkx_tab <- sparkx_tab[
  order(sparkx_tab$adjustedPval),
  ,
  drop = FALSE
]


# ============================================================
# Display and save results
# ============================================================

print(res_sparkx$stats)
print(res_sparkx$confusion)
print(runtime)

write.csv(
  sparkx_tab,
  sparkx_csv_file,
  row.names = FALSE
)

saveRDS(
  list(
    sparkx_fit = sparkx_fit,
    result = sparkx_tab,
    combined_pvalues = sparkx_tab$combinedPval,
    adjusted_pvalues = sparkx_tab$adjustedPval,
    pred = sparkx_tab$pred,
    svg_idx = sparkx_svg_idx,
    svg_genes = sparkx_svg_genes,
    evaluation = res_sparkx,
    cutoff = p_cutoff,
    runtime = runtime,
    seed = seed
  ),
  file = sparkx_rds_file
)

cat("SPARK-X results saved to:\n", sparkx_rds_file, "\n")