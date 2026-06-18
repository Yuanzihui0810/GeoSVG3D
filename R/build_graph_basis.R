#' Compute a graph Laplacian
#'
#' @param W A square adjacency matrix.
#' @param normalized Logical; whether to construct the normalized Laplacian.
#'
#' @return A graph Laplacian matrix.
#' @noRd
compute_laplacian <- function(W, normalized = FALSE) {
  d <- Matrix::rowSums(W)
  
  if (!normalized) {
    D <- Matrix::Diagonal(x = d)
    L <- D - W
  } else {
    # normalized Laplacian
    # d_inv_sqrt <- 1 / sqrt(d + 1e-8)
    d_inv_sqrt <- rep(0, length(d))
    idx <- d > 0
    d_inv_sqrt[idx] <- 1 / sqrt(d[idx])
    
    D_inv_sqrt <- Matrix::Diagonal(x = d_inv_sqrt)
    L <- Matrix::Diagonal(n = length(d)) - D_inv_sqrt %*% W %*% D_inv_sqrt
  }
  
  return(L)
}

#' Compute graph Laplacian eigenvectors
#'
#' @param L A graph Laplacian matrix.
#' @param K Number of eigenvectors.
#' @param tol Numerical tolerance.
#'
#' @return A list containing eigenvalues and eigenvectors.
#' @noRd
compute_eigen <- function(L, K, tol = 1e-10) {
  # find the smallest K non-zero eigenvalues (remove λ_*=0)
  eig <- RSpectra::eigs_sym(L, k = K + 1, which = "LM", sigma = 0)
  # using "LR" can get the same result as matrix L is a real symmetric matrice 
  # eig <- RSpectra::eigs_sym(L, k = K + 1, which = "SM")
  
  values <- eig$values
  vectors <- eig$vectors
  values <- values[1: K]
  vectors <- vectors[, 1: K, drop = FALSE]
  
  flag = (abs(values) <= tol)
  values[flag] <- tol
  
  return(list(
            lambda = values,
            U = vectors
          ))
}

#' Construct a spatial k-nearest-neighbor graph
#'
#' Constructs a symmetric spatial graph using binary weights or
#' Gaussian-kernel weights.
#'
#' @param S Numeric matrix with observations in rows and spatial
#'   coordinates in columns.
#' @param k Number of nearest neighbors.
#' @param h Gaussian-kernel bandwidth. When `NULL`, the median
#'   nearest-neighbor distance is used.
#' @param binary Logical; whether to use binary edge weights.
#'
#' @return A symmetric sparse adjacency matrix.
#'
#' @export
build_W_knn <- function(S, k = 100, h = NULL, binary = FALSE) {
  # S: n x 3 matrix (coordinates)
  n <- nrow(S)
  
  # kNN search (fast C implementation)
  nn <- RANN::nn2(S, k = k + 1)  # include itself
  
  idx <- nn$nn.idx[, -1, drop = FALSE] # remove itself
  dist <- nn$nn.dists[, -1, drop = FALSE]
  
  i <- rep(1:n, each = k)
  j <- as.vector(t(idx))
  
  if(!binary){
    # Gaussian kernel
    if (is.null(h)) {
      h = stats::median(dist)
    } else if (h == 0){
      stop("The bandwidth h should not be zero.")
    }
    
    weights <- exp(-(dist^2) / (h^2))
    x <- as.vector(t(weights))
  }else{
    x <- rep(1, length(i))
  }
  
  W <- Matrix::sparseMatrix(i = i, j = j, x = x, dims = c(n, n))
 
  W <- (W + Matrix::t(W)) / 2
  
  return(W)
}