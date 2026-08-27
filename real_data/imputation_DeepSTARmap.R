# In R:
library(Matrix)
alpha = 0.01
input_dir <- "D:/MERFISH/brain/Deep_STARMAP/STAGATE_imputation_output/Peptidergic_neurons"
out_dir <- "D:/github/GeoSVG3D/real_data/output"
prefix <- "Peptidergic_neurons"
Y <- readMM(file.path(input_dir, paste0(prefix, "_Y_STAGATE_ReX.mtx")))
spot_ids <- readLines(file.path(input_dir, paste0(prefix, "_spot_ids.tsv")))
gene_ids <- readLines(file.path(input_dir, paste0(prefix, "_gene_ids.tsv")))
rownames(Y) <- spot_ids
colnames(Y) <- gene_ids
S <- read.csv(file.path(input_dir, paste0(prefix, "_S_3d_coordinates.csv")),
              row.names = 1, check.names = FALSE)

cat("Files read sucessfully.\n")

t1 = Sys.time()

fit_res <- GeoSVG3D(
  Y = Y, S = S, n_basis = 5, knn = 100,
  n_iter = 2000, burn_in = 2000,
  kappa_g = 0,
  alpha_g = 10,
  xi_g0 = 0,
  prior = "SpSL-L",
  xi_prop_sd = 1,
  tau2_g_update = FALSE,
  tau2_g_init = 1,
  n_cores = 10,
  seed = 1
)

t2 = Sys.time()
print(t2 - t1)

p0g <- compute_p0g(fit_res, tol = 0)
bfdr_res <- select_svg_bfdr(p0g, alpha = alpha)

svg_tab <- data.frame(
  gene = colnames(Y),
  p0g = p0g,
  pred = bfdr_res$pred,
  stringsAsFactors = FALSE
)

svg_tab <- svg_tab[order(svg_tab$p0g), ]

svg_idx <- which(bfdr_res$pred == 1)
svg_genes <- colnames(Y)[svg_idx]

print(paste("The number of SVGs detected: ", length(svg_idx)))

saveRDS(
  list(
    fit_res = fit_res,
    p0g = p0g,
    bfdr_res = bfdr_res,
    svg_tab = svg_tab,
    svg_genes = svg_genes,
    runtime_geosvg3d = t2 - t1
  ),
  file = file.path(out_dir, paste0(prefix, "_GeoSVG3D_results.rds"))
)

write.csv(
  svg_tab,
  file = file.path(out_dir, paste0(prefix, "_GeoSVG3D_svg_table.csv")),
  row.names = FALSE
)