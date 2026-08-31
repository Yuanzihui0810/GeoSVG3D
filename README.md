# GeoSVG-3D: Bayesian detection of spatially variable genes in three-dimensional spatial transcriptomics

The R package **GeoSVG3D** is developed to implement GeoSVG-3D, a fully Bayesian framework for detecting spatially variable genes (SVGs) in three-dimensional spatial transcriptomics data. GeoSVG-3D constructs a neighborhood graph of cells based on their three-dimensional spatial coordinates and derives geometry-aware basis functions to characterize local connectivity and overall tissue organization. Gene expression is then modeled as a sparse expansion over these basis functions using neuronized priors, allowing complex spatial dependence to be learned adaptively from the data rather than being restricted to prespecified kernels, covariance structures, or spatial pattern classes. Bayesian false discovery rate control is further applied to gene-wise posterior evidence to determine the final SVG set.

More technical details can be found in:

Zihui Yuan., Yinqiao Yan, and Xiangyu Luo (2026) Geometry-aware Bayesian detection of spatially variable genes in three dimensions. *Statistical Learning and Data Science*. Accepted.

## Installation

``` r
if(!require(devtools))
  install.packages(devtools)
devtools::install_github("Yuanzihui0810/GeoSVG3D")
```

## Main functions

| Function            | Description                                                                                                                                                                                                                                                                                                                                                |
|------------------------------------|------------------------------------|
| `GeoSVG3D()`        | Main fitting function for GeoSVG-3D. It internally constructs a sparse weighted k-nearest-neighbor graph from the three-dimensional spatial coordinates using `build_W_knn()`, derives geometry-aware basis functions from the graph Laplacian, applies SPDE-inspired spectral scaling, and fits the gene-specific Bayesian model using neuronized priors. |
| `compute_p0g()`     | Computes the posterior probability of no spatial effect for each gene from a fitted GeoSVG-3D model.                                                                                                                                                                                                                                                       |
| `select_svg_bfdr()` | Selects spatially variable genes based on the gene-specific posterior null probabilities using Bayesian false discovery rate control.                                                                                                                                                                                                                      |

## Main arguments

The main fitting function expects:

-   `Y`: an $n \times G$ numeric expression matrix, with cells or spots in rows and genes in columns;
-   `S`: an $n \times 3$ numeric matrix of three-dimensional spatial coordinates.

The rows of `Y` and `S` must correspond to the same cells or spots and appear in the same order. Spatial coordinates are internally scaled to `[0, 1]` along each dimension.

The expression values should be continuous and appropriately preprocessed. For the real-data analysis presented in the paper, raw count data are processed by quality control, library-size normalization, and a `log1p` transformation. STAGATE is subsequently used to impute zero expression values before the resulting continuous expression matrix is supplied to GeoSVG3D. The complete preprocessing and imputation workflow is provided in real_data/STAGATE_zero_imputation_DeepSTARmap.ipynb.


The other main arguments of the fitting function `GeoSVG3D()` are described below.

-   `n_basis`: the number of graph Laplacian basis functions used to represent spatial variation.
-   `knn`: the number of nearest neighbors used to construct the spatial graph.
-   `h`: the bandwidth of the Gaussian kernel used to define graph weights. If `NULL`, the median of the corresponding nearest-neighbor distances is used.
-   `n_iter`: the number of retained MCMC iterations after burn-in.
-   `burn_in`: the number of burn-in MCMC iterations.
-   `thin`: the thinning value for retained MCMC samples. "thin = n" means the interval between the two retained posterior samples is n-1.
-   `kappa_g`: the spectral range parameter used in SPDE-inspired scaling.
-   `alpha_g`: the spectral smoothness parameter used in SPDE-inspired scaling.
-   `prior`: the neuronized prior specification, including Bayesian Lasso-type (`"BL"`), horseshoe-type (`"HS"`), and spike-and-slab priors,  where the spike is a Dirac point mass at zero ($\delta_0$) and the slab is Laplace-like (`"SpSL-L"`), Cauchy (`"SpSL-C"`), or Gaussian (`"SpSL-G"`).
-   `n_cores`: the number of CPU cores used for gene-wise parallel computation.
-   `seed`: the random seed used for reproducible computation.

## Quick start

The following example uses a small simulated dataset to illustrate the package interface.

``` r
library(GeoSVG3D)

set.seed(123)

n <- 200
G <- 20

Y <- matrix(rnorm(n * G), nrow = n, ncol = G)
S <- matrix(runif(n * 3), nrow = n, ncol = 3, dimnames = list(NULL, c("x", "y", "z")))

fit <- GeoSVG3D(
  Y = Y, S = S, n_basis = 5, knn = 100, n_iter = 2000, burn_in = 2000, thin = 1,
  prior = "SpSL-L", n_cores = 1, seed = 123
)

p0g <- compute_p0g(fit, tol = 0)

bfdr_result <- select_svg_bfdr(p0g = p0g, alpha = 0.05)

which(bfdr_result$pred == 1L)
```

## Parallel computation

GeoSVG-3D fits the model separately for each gene, allowing gene-wise computations to be performed in parallel.

``` r
fit <- GeoSVG3D(
  Y = Y, S = S, n_basis = 5, knn = 100, n_iter = 2000, burn_in = 2000, thin = 1,
  prior = "SpSL-L", n_cores = 4, seed = 123
)
```

## Datasets

The Deep-STARMAP dataset used in the manuscript is available at <https://zenodo.org/records/16783355>.

## Remarks

If you have any questions regarding this package, please contact Zihui Yuan at [yuanzihui0810@gmail.com](mailto:yuanzihui0810@gmail.com).