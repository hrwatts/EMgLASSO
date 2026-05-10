#' E-Step of the EM Algorithm
#'
#' Computes posterior responsibilities (cluster membership probabilities) for each observation
#' given the current parameter estimates. This is the expectation step of the EM algorithm.
#'
#' @param x An \eqn{n \times d}{n x d} matrix of observations.
#' @param Theta A list containing current parameter estimates with elements `Tau`
#'   (mixture proportions), `Mu` (component mean vectors), and `Sigma`
#'   (component covariance matrices).
#' @param W A list of covariance matrices to use for density evaluation (typically equals Theta$Sigma).
#'   Used to allow flexible specification of the covariance structure.
#'
#' @return A list with `LogLike`, an \eqn{n \times k}{n x k} matrix of weighted densities,
#'   and `Tnk`, an \eqn{n \times k}{n x k} matrix of posterior responsibilities.
#'
#' @details
#' The posterior responsibility is computed via Bayes' rule:
#' \deqn{T_{i,j} = \frac{\tau_j f(x_i; \mu_j, \Sigma_j)}{\sum_{\ell=1}^{k} \tau_\ell f(x_i; \mu_\ell, \Sigma_\ell)}}
#'
#' @importFrom mvtnorm dmvnorm
#' @keywords internal

estep <- function(x, Theta, W) {
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (!is.numeric(x) || nrow(x) < 1 || ncol(x) < 1) {
    stop("x must be a non-empty numeric matrix.")
  }
  if (!is.list(Theta) || is.null(Theta$Tau) || is.null(Theta$Mu) || is.null(Theta$Sigma)) {
    stop("Theta must contain Tau, Mu, and Sigma.")
  }

  K <- length(Theta$Tau)
  if (K < 1 || length(Theta$Mu) != K || length(W) != K) {
    stop("Tau, Mu, and W must have the same positive length.")
  }

  LogLike <- with(
    Theta,
    do.call(cbind, lapply(seq_len(K), function(k) Tau[k] * dmvnorm(x, Mu[[k]], W[[k]])))
  )

  denom <- rowSums(LogLike)
  bad_rows <- !is.finite(denom) | denom <= 0

  Tnk <- LogLike / denom
  if (any(bad_rows)) {
    Tnk[bad_rows, ] <- 1 / K
  }
  Tnk[!is.finite(Tnk)] <- 1 / K

  return(list(LogLike = LogLike, Tnk = Tnk))
}
