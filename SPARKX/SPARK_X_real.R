library(SPARK)
library(Matrix)

# ============================================================
# 1. Settings
# ============================================================

data_dir <- "D:/MERFISH/brain/Deep_STARMAP"

expr_file <- file.path(
  data_dir,
  "Brain_Deep_STARmap_expression_matrix.csv"
)

spatial_file <- file.path(
  data_dir,
  "Brain_Deep_STARmap_spatial.csv"
)

out_dir <- "D:/github/GeoSVG3D/SPARKX/output"

target_celltype <- "Peptidergic neurons"

dir.create(
  out_dir,
  recursive = TRUE,
  showWarnings = FALSE
)


# ============================================================
# 2. Read data
# ============================================================

expr <- read.csv(
  expr_file,
  check.names = FALSE
)

spatial <- read.csv(
  spatial_file,
  check.names = FALSE
)


# ============================================================
# 3. Select Peptidergic neurons
# ============================================================

keep <- spatial$FUSEmap_main_level == target_celltype

expr_sub <- expr[
  keep,
  ,
  drop = FALSE
]

spatial_sub <- spatial[
  keep,
  ,
  drop = FALSE
]


# ============================================================
# 4. Prepare SPARK-X input
#
# count_mat:
#   genes x cells
#
# location:
#   cells x 3
# ============================================================

count_mat <- Matrix(
  t(as.matrix(expr_sub)),
  sparse = TRUE
)

location <- as.matrix(
  spatial_sub[
    ,
    c("x", "y", "z")
  ]
)

storage.mode(location) <- "double"


# Add consistent cell IDs

cell_ids <- paste0(
  "cell_",
  seq_len(nrow(location))
)

rownames(count_mat) <- colnames(expr_sub)
colnames(count_mat) <- cell_ids

rownames(location) <- cell_ids
colnames(location) <- c(
  "x",
  "y",
  "z"
)


cat(
  "SPARK-X input:",
  nrow(count_mat), "genes x",
  ncol(count_mat), "cells\n"
)


# ============================================================
# 5. Run SPARK-X
# ============================================================

t1 <- Sys.time()

sparkx_res <- sparkx(
  count_in = count_mat,
  locus_in = location,
  numCores = 1,
  option = "mixture",
  verbose = FALSE
)

t2 <- Sys.time()

runtime_sparkx <- t2 - t1

print(runtime_sparkx)


# ============================================================
# 6. Extract results
# ============================================================

sparkx_tab <- data.frame(
  gene = rownames(sparkx_res$res_mtest),
  sparkx_res$res_mtest,
  row.names = NULL
)

sparkx_tab$pred <- as.integer(
  sparkx_tab$adjustedPval <= 0.05
)

sparkx_tab <- sparkx_tab[
  order(sparkx_tab$adjustedPval),
  ,
  drop = FALSE
]

svg_genes <- sparkx_tab$gene[
  sparkx_tab$pred == 1
]

cat(
  "The number of SVGs detected by SPARK-X:",
  length(svg_genes),
  "\n"
)


# ============================================================
# 7. Save results
# ============================================================

write.csv(
  sparkx_tab,
  file = file.path(
    out_dir,
    "Peptidergic_neurons_SPARKX_3D_svg_table.csv"
  ),
  row.names = FALSE
)

saveRDS(
  list(
    sparkx_res = sparkx_res,
    result = sparkx_tab,
    svg_genes = svg_genes,
    runtime_sparkx = runtime_sparkx
  ),
  file = file.path(
    out_dir,
    "Peptidergic_neurons_SPARKX_3D_results.rds"
  )
)