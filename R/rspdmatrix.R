#' Generate a Random Sparse Positive Definite Matrix
#'
#' Generates a random symmetric positive definite matrix with approximately a specified
#' level of sparsity. Useful for creating synthetic covariance matrices for simulation studies.
#'
#' @param D Integer. Dimension of the matrix (number of rows and columns).
#' @param lambda Numeric scalar between 0 and 1. Approximate sparsity level
#'   (fraction of zero entries).
#'   A value of 0 produces a dense matrix; a value of 1 produces mostly zeros.
#' @param epsilon Numeric. Small positive value added to diagonal to ensure positive definiteness.
#'   Default: \code{1e-4}.
#'
#' @return A \eqn{D \times D}{D x D} symmetric positive definite matrix with approximate
#'   sparsity \eqn{\lambda}{lambda}.
#'
#' @details
#' The algorithm:
#' 1. Fill the lower triangle with random uniform values in \eqn{(0,2)}.
#' 2. Randomly set entries to zero with probability \eqn{\lambda}{lambda}.
#' 3. Mirror to upper triangle to create symmetry.
#' 4. Add \eqn{\epsilon + r_i}{epsilon + row_sum_i} to diagonal entry \eqn{(i,i)},
#'    where \eqn{r_i}{row_sum_i} is the sum of row \eqn{i}, ensuring positive definiteness.
#'
#' This produces sparse matrices suitable for testing graphical LASSO and other
#' sparsity-inducing methods, though the sparsity pattern differs from the near-diagonal
#' structure typically learned by graphical LASSO from real data.
#'
#' @examples
#' \dontrun{
#' # Generate a 5x5 sparse positive definite matrix with ~30% sparsity
#' M <- rspdmatrix(D = 5, lambda = 0.3)
#' all(eigen(M)$values > 0)  # Verify positive definiteness
#' }
#'
#' @importFrom stats rbinom runif
#' @export

rspdmatrix <- function(D, lambda, epsilon = 1e-4) {
  if (!is.numeric(D) || length(D) != 1 || !is.finite(D) || D < 1) {
    stop("D must be a positive scalar.")
  }
  if (!is.numeric(lambda) || length(lambda) != 1 || !is.finite(lambda) || lambda < 0 || lambda > 1) {
    stop("lambda must be a finite scalar in [0, 1].")
  }
  if (!is.numeric(epsilon) || length(epsilon) != 1 || !is.finite(epsilon) || epsilon <= 0) {
    stop("epsilon must be a positive finite scalar.")
  }

  D <- as.integer(D)
  rMat <- matrix(rep(NA, D * D), nrow = D)

  pSum <- (D * (D + 1)) / 2

  rMat[lower.tri(rMat, diag = TRUE)] <- stats::runif(pSum, 0, 2)

  ZeroIndex <- stats::rbinom(D * D, 1, lambda)

  rMat[ZeroIndex == 1] <- 0

  rMat[upper.tri(rMat)] <- t(rMat)[upper.tri(rMat)]

  Eigens <- rowSums(rMat) + epsilon

  diag(rMat) <- Eigens

  return(rMat)
}
