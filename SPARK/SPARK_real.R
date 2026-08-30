library(SPARK)

#extend the function "spark.test" to 3 dimensions.
spark.test.new <- function(object, 
                       kernel_mat = NULL, 
                       check_positive = TRUE, 
                       verbose = TRUE) {
  
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
  }else{# kernel_mat is a list provided by user
    for(ikernel in 1:length(kernel_mat) ){
      # pre-defined kernels
      cat(paste0("## testing pre-defined kernel: ",ikernel,"...\n"))
      object <- spark.test_each(object, kernel_mat=kernel_mat[[ikernel]], check_positive=check_positive, verbose=verbose)
      res_pval <- cbind(res_pval, object@res_stest$sw)
      res_kernel <- c(res_kernel, object@res_stest)
    }# end for ikernel
    colnames(res_pval) <- paste0("kernel", 1:length(kernel_mat) )
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

out_dir <- "D:/github/GeoSVG3D/SPARK/output"

target_celltype <- "Peptidergic neurons"
celltype_col <- "FUSEmap_main_level"

x_col <- "x"
y_col <- "y"
z_col <- "z"

prefix <- "Peptidergic_neurons"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# ============================================================
# 2. Read data and select Peptidergic neurons
# ============================================================

expr <- read.csv(
  expr_file,
  # row.names = 1,
  check.names = FALSE
)

spatial <- read.csv(
  spatial_file,
  # row.names = 1,
  check.names = FALSE
)

# Select Peptidergic neurons by the same row indices
keep <- spatial[[celltype_col]] == target_celltype

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

# Give the selected cells common IDs
cell_ids <- paste0(
  "cell_",
  seq_len(nrow(expr_sub))
)

# ============================================================
# 3. Prepare SPARK input
#
# counts:   genes x cells
# location: cells x 3
# ============================================================

# SPARK count matrix: genes x cells
counts <- t(
  as.matrix(expr_sub)
)

colnames(counts) <- cell_ids

# SPARK coordinates: cells x 3
location <- spatial_sub[
  ,
  c(x_col, y_col, z_col),
  drop = FALSE
]

rownames(location) <- cell_ids
colnames(location) <- c("x", "y", "z")

stopifnot(
  identical(
    colnames(counts),
    rownames(location)
  )
)


cat(
  "SPARK input:",
  nrow(counts), "genes x",
  ncol(counts), "cells\n"
)


# ============================================================
# 4. Create SPARK object
#
# Disable additional SPARK filtering
# ============================================================

spark <- CreateSPARKObject(
  counts = counts,
  location = location,
  percentage = 0,
  min_total_counts = 0
)

spark@lib_size <- apply(spark@counts, 2, sum)

# Release unused objects before constructing large kernels

rm(
  expr,
  spatial,
  expr_sub,
  spatial_sub,
  counts,
  location
)

invisible(gc())


# ============================================================
# 5. Fit count-based model under H0
# ============================================================

t1 <- Sys.time()

spark <- spark.vc(
  spark,
  covariates = NULL,
  lib_size = spark@lib_size,
  num_core = 10,
  verbose = FALSE
)

# ============================================================
# 7. Test spatial patterns
# ============================================================

spark <- spark.test.new(
  spark,
  check_positive = TRUE,
  verbose = FALSE
)

t2 <- Sys.time()

runtime_spark <- t2 - t1

print(runtime_spark)

# ============================================================
# 7. Extract results
# ============================================================

spark_res <- data.frame(
  gene = rownames(spark@res_mtest),
  spark@res_mtest,
  row.names = NULL,
  check.names = FALSE
)

spark_res$pred <- as.integer(
  spark_res$adjusted_pvalue <= 0.05
)

spark_res <- spark_res[
  order(spark_res$adjusted_pvalue),
  ,
  drop = FALSE
]

svg_genes <- spark_res$gene[
  spark_res$pred == 1
]

cat(
  "The number of SVGs detected by SPARK:",
  length(svg_genes),
  "\n"
)


# ============================================================
# 8. Save results
# ============================================================

write.csv(
  spark_res,
  file = file.path(
    out_dir,
    paste0(
      prefix,
      "_SPARK_3D_svg_table.csv"
    )
  ),
  row.names = FALSE
)

saveRDS(
  list(
    spark = spark,
    result = spark_res,
    svg_genes = svg_genes,
    runtime_spark = runtime_spark,
    kernel_length_scales = lrang
  ),
  file = file.path(
    out_dir,
    paste0(
      prefix,
      "_SPARK_3D_results.rds"
    )
  )
)