library(Matrix)
library(data.table)
library(ggplot2)
library(scatterplot3d)

# ============================================================
# 0. Paths
# ============================================================

# Original Deep STARmap data
rawdata_dir <- "D:/MERFISH/brain/Deep_STARMAP"

# STAGATE-imputed data used as model input
imputation_dir <- file.path(
  rawdata_dir,
  "STAGATE_imputation_output",
  "Peptidergic_neurons"
)

prefix <- "Peptidergic_neurons"

# Output directories in the GitHub project
project_dir <- "D:/github/GL-neuronized-svg"
real_data_dir <- file.path(project_dir, "real_data")

tab_dir <- file.path(real_data_dir, "table")

dir.create(tab_dir, recursive = TRUE, showWarnings = FALSE)

main_blue <- "#1F4E79"

# ============================================================
# 1. Load STAGATE-imputed Peptidergic neurons data
#    Y: cells x genes
#    S: cells x 3 coordinates
# ============================================================

Y <- readMM(
  file.path(
    imputation_dir,
    paste0(prefix, "_Y_STAGATE_ReX.mtx")
  )
)

spot_ids <- readLines(
  file.path(
    imputation_dir,
    paste0(prefix, "_spot_ids.tsv")
  )
)

gene_ids <- readLines(
  file.path(
    imputation_dir,
    paste0(prefix, "_gene_ids.tsv")
  )
)

spot_ids_R <- as.integer(spot_ids) + 1L
cell_ids_stagate <- paste0("cell_", spot_ids_R)

rownames(Y) <- cell_ids_stagate
colnames(Y) <- gene_ids

S <- read.csv(
  file.path(
    imputation_dir,
    paste0(prefix, "_S_3d_coordinates.csv")
  ),
  row.names = 1,
  check.names = FALSE
)

# Make coordinate row names consistent with Y
rownames(S) <- cell_ids_stagate

stopifnot(nrow(Y) == nrow(S))
stopifnot(all(rownames(Y) == rownames(S)))

S <- as.matrix(S[, 1:3])
storage.mode(S) <- "double"

# Use the same scaled coordinates as GeoSVG-3D
S_scaled <- scale_spatial_coords_01(S)

cat("Imputed Y:", dim(Y), "\n")
cat("Coordinates:", dim(S_scaled), "\n")

# ============================================================
# 2. Load raw expression for before-imputation EDA
# ============================================================

expr_file <- file.path(
  rawdata_dir,
  "Brain_Deep_STARmap_expression_matrix.csv"
)

spatial_file <- file.path(
  rawdata_dir,
  "Brain_Deep_STARmap_spatial.csv"
)

Raw_Exp <- fread(
  expr_file,
  data.table = FALSE
)

Spatial_raw <- fread(
  spatial_file,
  data.table = FALSE
)

stopifnot(nrow(Raw_Exp) == nrow(Spatial_raw))

cell_id <- paste0("cell_", seq_len(nrow(Raw_Exp)))

rownames(Raw_Exp) <- cell_id
rownames(Spatial_raw) <- cell_id

idx_pep <- which(
  Spatial_raw$FUSEmap_main_level == "Peptidergic neurons"
)

Raw_Exp_pep <- Raw_Exp[
  idx_pep,
  ,
  drop = FALSE
]

# Align raw expression to the imputed cells and genes
common_cells <- intersect(
  rownames(Y),
  rownames(Raw_Exp_pep)
)

common_genes <- intersect(
  colnames(Y),
  colnames(Raw_Exp_pep)
)

cat("Matched raw cells:", length(common_cells), "\n")
cat("Matched raw genes:", length(common_genes), "\n")

Y_raw <- Raw_Exp_pep[
  rownames(Y),
  colnames(Y),
  drop = FALSE
]

Y_raw <- Matrix(
  as.matrix(Y_raw),
  sparse = TRUE
)

# ============================================================
# 3. EDA helper functions
# ============================================================

true_nonzero_count <- function(X, tol = 1e-12) {
  if (inherits(X, "Matrix")) {
    return(sum(abs(X@x) > tol, na.rm = TRUE))
  }
  
  sum(abs(X) > tol, na.rm = TRUE)
}

zero_prop <- function(X, tol = 1e-12) {
  1 - true_nonzero_count(X, tol = tol) / prod(dim(X))
}

row_nonzero <- function(X, tol = 1e-12) {
  if (inherits(X, "Matrix")) {
    return(Matrix::rowSums(abs(X) > tol))
  }
  
  rowSums(abs(X) > tol, na.rm = TRUE)
}

col_nonzero <- function(X, tol = 1e-12) {
  if (inherits(X, "Matrix")) {
    return(Matrix::colSums(abs(X) > tol))
  }
  
  colSums(abs(X) > tol, na.rm = TRUE)
}

make_cols <- function(
    v,
    n_col = 100,
    palette = "viridis"
) {
  pal <- hcl.colors(n_col, palette)
  
  if (all(is.na(v)) ||
      diff(range(v, na.rm = TRUE)) == 0) {
    return(rep("grey50", length(v)))
  }
  
  idx <- cut(
    v,
    breaks = seq(
      min(v, na.rm = TRUE),
      max(v, na.rm = TRUE),
      length.out = n_col + 1
    ),
    include.lowest = TRUE,
    labels = FALSE
  )
  
  pal[idx]
}

# ============================================================
# 4. EDA summary table
# ============================================================

eda_summary <- data.frame(
  metric = c(
    "n_cells",
    "n_genes",
    "raw_zero_proportion",
    "imputed_zero_proportion",
    "median_detected_genes_per_cell_raw",
    "median_detected_genes_per_cell_imputed",
    "median_gene_detection_fraction_raw",
    "median_gene_detection_fraction_imputed"
  ),
  value = c(
    nrow(Y),
    ncol(Y),
    zero_prop(Y_raw),
    zero_prop(Y),
    median(row_nonzero(Y_raw)),
    median(row_nonzero(Y)),
    median(col_nonzero(Y_raw) / nrow(Y_raw)),
    median(col_nonzero(Y) / nrow(Y))
  )
)

print(eda_summary)

write.csv(
  eda_summary,
  file.path(
    tab_dir,
    paste0(prefix, "_EDA_summary.csv")
  ),
  row.names = FALSE
)
