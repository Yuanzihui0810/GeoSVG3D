# install.packages("scBSP")
library(scBSP)
library(data.table)
library(Seurat)
#------------------------------simulated data --------------------------------------
# data_dir     <- "D:/github/GeoSVG3D/simulation/data"
# expr_file    <- file.path(data_dir, "log_counts.csv")
# spatial_file <- file.path(data_dir, "coords.csv")
# 
# Log_Exp <- fread(expr_file,    data.table = FALSE)
# Coords  <- fread(spatial_file, data.table = FALSE)
# 
# rownames(Log_Exp) <- Log_Exp[[1]];  Log_Exp <- Log_Exp[, -1]
# rownames(Coords)  <- Coords[[1]];   Coords  <- Coords[,  -1]
# 
# stopifnot(nrow(Log_Exp) == nrow(Coords))
# stopifnot(all(rownames(Log_Exp) == rownames(Coords))) 
# 
# Coords  <- as.matrix(Coords[, 1:3])       # spots × 3，x y z
# Log_Exp <- t(as.matrix(Log_Exp))          # genes × spots
# Log_Exp <- Matrix(Log_Exp, sparse = TRUE)

seed <- 1
p_cutoff <- 0.05

# RDS
sim_rds_file <- file.path("D:/github/GeoSVG3D/simulation/data", 
                          paste0("simulation_data_seed", seed, ".rds"))
# scBSP results are saved to:
scbsp_dir <- "D:/github/GeoSVG3D/scBSP"
dir.create(scbsp_dir, recursive = TRUE, showWarnings = FALSE)
scbsp_rds_file <- file.path(scbsp_dir, 
                            paste0("scbsp_simulation_results_seed", seed, ".rds"))

# Read data from rds
sim <- readRDS(sim_rds_file)

Coords <- as.matrix(sim$S[, 1:3, drop = FALSE])  # spots × 3
Log_Exp <- Matrix(t(sim$Y), sparse = TRUE)       # genes × spots

if (!is.null(colnames(sim$Y))) {
  rownames(Log_Exp) <- colnames(sim$Y)
}

t1 <- Sys.time()
P_values <- scBSP(Coords, Log_Exp, Exp_Norm = FALSE)
t2 <- Sys.time()
runtime <- t2 - t2

scbsp_pvalues <- P_values[, 2]
scbsp_pred <- as.integer(scbsp_pvalues < p_cutoff)

scbsp_svg_idx <- which(scbsp_pred == 1)
scbsp_svg_genes <- as.character(P_values[scbsp_svg_idx, 1])

res_scbsp <- evaluate_svg_detection(
  truth = sim$is_svg,
  pred = scbsp_pred
)

print(res_scbsp$stats)
print(res_scbsp$confusion)
print(runtime)

saveRDS(
  list(
    P_values = P_values,
    pvalues = scbsp_pvalues,
    pred = scbsp_pred,
    svg_idx = scbsp_svg_idx,
    svg_genes = scbsp_svg_genes,
    evaluation = res_scbsp,
    cutoff = p_cutoff,
    runtime = runtime,
    seed = seed
  ),
  file = scbsp_rds_file
)

cat("scBSP results saved to:\n", scbsp_rds_file, "\n")


#------------------------------real data --------------------------------------
# ============================================================
# 1. Load raw Deep-STARmap data
# ============================================================
#it is raw data
data_dir <- "D:/MERFISH/brain/Deep_STARMAP"

expr_file    <- file.path(data_dir, "Brain_Deep_STARmap_expression_matrix.csv")
spatial_file <- file.path(data_dir, "Brain_Deep_STARmap_spatial.csv")

expr <- fread(expr_file, data.table = FALSE)      # cells x genes
spatial <- fread(spatial_file, data.table = FALSE)

stopifnot(nrow(expr) == nrow(spatial))

# ============================================================
# 2. Create consistent cell IDs
#    This mimics Python row index behavior.
# ============================================================

cell_id <- paste0("cell_", seq_len(nrow(expr)))

rownames(expr) <- cell_id
rownames(spatial) <- cell_id

# ============================================================
# 3. Subset Peptidergic neurons
# ============================================================

cell_type <- "Peptidergic neurons"

idx_cell <- which(spatial$FUSEmap_main_level == cell_type)

expr_sub <- expr[idx_cell, , drop = FALSE]
spatial_sub <- spatial[idx_cell, , drop = FALSE]

stopifnot(nrow(expr_sub) == nrow(spatial_sub))
stopifnot(all(rownames(expr_sub) == rownames(spatial_sub)))

cat("Selected cell type:", cell_type, "\n")
cat("Number of cells:", nrow(expr_sub), "\n")
cat("Number of genes before scBSP filtering:", ncol(expr_sub), "\n")

# ============================================================
# 4. Prepare coordinates and expression matrix for scBSP
#    scBSP expects:
#    Coords: spots/cells x dimensions
#    Raw_Exp: genes x spots/cells
# ============================================================

Coords_sub <- as.matrix(spatial_sub[, 1:3])
storage.mode(Coords_sub) <- "double"

Raw_Exp_sub <- t(as.matrix(expr_sub))   # genes x cells
Raw_Exp_sub <- Matrix(Raw_Exp_sub, sparse = TRUE)

# Keep gene names
rownames(Raw_Exp_sub) <- colnames(expr_sub)
colnames(Raw_Exp_sub) <- rownames(expr_sub)

# ============================================================
# 5. Run scBSP filtering and p-value calculation
# ============================================================

Filtered_ExpMat <- SpFilter(Raw_Exp_sub)

cat("Number of genes after SpFilter:", nrow(Filtered_ExpMat), "\n")

t1 <- Sys.time()

# scbsp_bfdr_res <- which(P_values[,2]<0.05)
# scbsp_svg_genes <- P_values[scbsp_bfdr_res,1]
P_values <- scBSP(
  Coords_sub,
  Filtered_ExpMat
)

t2 <- Sys.time()
runtime_scbsp <- difftime(t2, t1, units = "mins")

cat("scBSP runtime:", runtime_scbsp, "minutes\n")

# ============================================================
# 6. Organize scBSP results
# ============================================================

P_values <- as.data.frame(P_values)

# Usually scBSP returns columns like gene and p-value.
# This standardizes the first two columns.
colnames(P_values)[1:2] <- c("gene", "scbsp_pvalue")

P_values$scbsp_pvalue <- as.numeric(P_values$scbsp_pvalue)
P_values$scbsp_pred <- as.integer(P_values$scbsp_pvalue < 0.05)

# P_values <- P_values[order(P_values$scbsp_pvalue), ]
scbsp_bfdr_res <- which(P_values[,2]<0.05)
scbsp_svg_genes <- P_values$gene[P_values$scbsp_pred == 1]

cat("Number of scBSP SVGs:", length(scbsp_svg_genes), "\n")

# ============================================================
# 7. Save scBSP results
# ============================================================

save_dir <- file.path("D:/github/GeoSVG3D/scBSP", gsub(" ", "_", cell_type))
dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

write.csv(
  P_values,
  file = file.path(save_dir, "scBSP_pvalues_Peptidergic_neurons.csv"),
  row.names = FALSE
)

saveRDS(
  list(
    cell_type = cell_type,
    cell_ids = rownames(expr_sub),
    genes_before_filter = colnames(expr_sub),
    genes_after_filter = rownames(Filtered_ExpMat),
    Coords_sub = Coords_sub,
    P_values = P_values,
    scbsp_svg_genes = scbsp_svg_genes,
    runtime_scbsp = runtime_scbsp
  ),
  file = file.path(save_dir, "scBSP_results_Peptidergic_neurons.rds")
)
