# GeoSVG3D

The method constructs a spatial graph from observed 3D coordinates, derives graph Laplacian eigenvectors as data-adaptive spatial basis functions, and represents each gene-specific spatial effect as a sparse combination of these bases. Neuronized priors are used for Bayesian variable selection, and Bayesian false discovery rate control is applied to posterior gene-level evidence to determine the final SVG set.

## Main features

- Designed for three-dimensional spatial transcriptomics data.
- Learns spatial basis functions from the observed tissue geometry.
- Uses a sparse weighted k-nearest-neighbor graph and graph Laplacian eigenvectors.
- Applies SPDE-inspired spectral scaling to distinguish low- and high-frequency spatial components.
- Supports several neuronized prior specifications, including:
  - Bayesian Lasso-type shrinkage (`"BL"`);
  - horseshoe-type shrinkage (`"HS"`);
  - spike-and-slab-type specifications (`"SpSL-L"`, `"SpSL-C"`, and `"SpSL-G"`).
- Uses MCMC sampling for posterior inference.
- Computes posterior probabilities of no spatial effect for each gene.
- Selects SVGs using Bayesian FDR control.
- Supports gene-wise parallel computation.

## Input data

The main fitting function expects:

- `Y`: an \(n \times G\) numeric expression matrix, with cells or spots in rows and genes in columns;
- `S`: an \(n \times 3\) numeric matrix of three-dimensional spatial coordinates.

The expression values should be continuous and appropriately preprocessed. With the raw count data and spatial coordinates, a typical workflow may include quality control, library-size normalization, and a `log1p` transformation. STAGATE is subsequently used to impute zero expression values before the resulting continuous expression matrix was supplied to GeoSVG3D. The complete preprocessing and imputation workflow is provided in real_data/STAGATE_zero_imputation_DeepSTARmap.ipynb. 

The rows of `Y` and `S` must correspond to the same cells or spots and must appear in the same order.

## Quick start

The following example uses a small simulated dataset to illustrate the package interface.

```r
library(GeoSVG3D)

set.seed(123)

n <- 100
G <- 20

S <- matrix(
  runif(n * 3),
  nrow = n,
  ncol = 3,
  dimnames = list(NULL, c("x", "y", "z"))
)

Y <- matrix(
  rnorm(n * G),
  nrow = n,
  ncol = G
)

fit <- fit_all_genes(
  Y = Y,
  S = S,
  K = 5,
  knn = 10,
  N = 2000,
  BURN = 2000,
  thin = 1,
  prior = "SpSL-L",
  ncores = 1,
  seed = 123
)

p0g <- compute_p0g(fit, tol = 0)

bfdr_result <- select_svg_bfdr(
  p0g = p0g,
  alpha = 0.05
)

which(bfdr_result$pred == 1L)
```
## Parallel computation

Genes are conditionally independent in the fitting stage and can be processed in parallel:

```r
fit <- fit_all_genes(
  Y = Y,
  S = S,
  K = 5,
  knn = 10,
  N = 2000,
  BURN = 2000,
  prior = "BL",
  ncores = 4,
  seed = 123
)
```
## Datasets

The Deep-STARmap dataset used in the manuscript is available from [Zenodo](https://zenodo.org/records/16783355).