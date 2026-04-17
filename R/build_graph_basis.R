library(Matrix)     
library(RANN)       
library(RSpectra)
library(Rcpp)
library(parallel)

compute_laplacian <- function(W, normalized = FALSE) {
  d <- rowSums(W)
  
  if (!normalized) {
    D <- Diagonal(x = d)
    L <- D - W
  } else {
    # normalized Laplacian
    d_inv_sqrt <- 1 / sqrt(d + 1e-8)
    D_inv_sqrt <- Diagonal(x = d_inv_sqrt)
    L <- Diagonal(n = length(d)) - D_inv_sqrt %*% W %*% D_inv_sqrt
  }
  
  return(L)
}

compute_eigen <- function(L, K) {
  # find the smallest K+1 eigenvalues (remove λ0=0)
  eig <- RSpectra::eigs_sym(L, k = K + 1, which = "SM")
  
  values <- eig$values
  vectors <- eig$vectors
  
  # remove the first trivial eigenvector
  list(
    lambda = values[2:(K + 1)],
    U = vectors[, 2:(K + 1), drop = FALSE]
  )
}

build_W_knn <- function(S, k = 10, h = NULL, symmetric = TRUE) {
  # S: n x 3 matrix (coordinates)
  n <- nrow(S)
  
  # kNN search (fast C implementation)
  nn <- RANN::nn2(S, k = k + 1)  # include itself
  
  idx <- nn$nn.idx[, -1, drop = FALSE] # remove itself
  dist <- nn$nn.dists[, -1, drop = FALSE]
  
  if (is.null(h)) {
    h <- median(dist)
  }
  
  # Gaussian kernel
  weights <- exp(-(dist^2) / (h^2))
  
  i <- rep(1:n, each = k)
  j <- as.vector(t(idx))
  x <- as.vector(t(weights))
  
  W <- sparseMatrix(i = i, j = j, x = x, dims = c(n, n))
  
  if (symmetric) {
    W <- (W + t(W)) / 2
  }
  
  return(W)
}