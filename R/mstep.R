#' M-Step of the EM Algorithm
#'
#' Maximizes the expected complete log-likelihood to update parameter estimates.
#' This is the maximization step of the EM algorithm.
#'
#' @param x An \eqn{n \times d}{n x d} matrix of observations.
#' @param Tnk An \eqn{n \times k}{n x k} matrix of posterior responsibilities from the E-step,
#'   where \eqn{T_{i,j}} is the probability that observation \eqn{i} belongs to component \eqn{j}.
#'
#' @return A list containing updated parameter estimates:
#'   \item{Tau}{Updated mixture proportions, a numeric vector of length \eqn{k}}
#'   \item{Mu}{List of updated mean vectors, each of length \eqn{d}}
#'   \item{Sigma}{List of updated covariance matrices, each \eqn{d \times d}{d x d}}
#'
#' @details
#' The updates are computed as:
#' \deqn{\tau_j^{(t+1)} = \frac{1}{n} \sum_{i=1}^{n} T_{i,j}}
#' \deqn{\mu_j^{(t+1)} = \frac{\sum_{i=1}^{n} T_{i,j} x_i}{\sum_{i=1}^{n} T_{i,j}}}
#' \deqn{\Sigma_j^{(t+1)} = \frac{\sum_{i=1}^{n} T_{i,j} (x_i - \mu_j^{(t+1)})(x_i - \mu_j^{(t+1)})'}{\sum_{i=1}^{n} T_{i,j}}}
#'
#' @importFrom stats cov.wt
#' @keywords internal

mstep <- function(x, Tnk) {
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (!is.numeric(x) || nrow(x) < 2 || ncol(x) < 1) {
    stop("x must be a numeric matrix with at least 2 rows.")
  }
  if (!is.matrix(Tnk) || !is.numeric(Tnk) || nrow(Tnk) != nrow(x) || ncol(Tnk) < 1) {
    stop("Tnk must be a numeric matrix with nrow(Tnk) == nrow(x) and at least one column.")
  }
  if (!all(is.finite(Tnk)) || any(Tnk < 0)) {
    stop("Tnk must contain finite non-negative values.")
  }
  if (any(colSums(Tnk) <= 0)) {
    stop("Each column of Tnk must have positive sum.")
  }

  X <- lapply(
    seq_len(ncol(Tnk)),
    function(k) cov.wt(x, Tnk[, k], method = "ML")
  )

  Tau <- colMeans(Tnk)
  Mu <- lapply(X, "[[", "center")
  Sigma <- lapply(X, "[[", "cov")
  Theta <- list(Mu = Mu, Sigma = Sigma, Tau = Tau)

  return(Theta)
}
