#' EM Algorithm for Gaussian Mixture Models with Sparse Covariance Estimation
#'
#' Fits a Gaussian mixture model using the EM algorithm and estimates sparse
#' covariance matrices for each component via graphical LASSO with BIC-selected penalties.
#'
#' @param x An \eqn{n \times d}{n x d} matrix of observations from a \eqn{d}-dimensional
#'   Gaussian mixture model with \eqn{n} observations.
#' @param iTau A numeric vector of length \eqn{k} specifying initial mixture proportions.
#'   Should sum to 1. If not, it will be normalized.
#' @param Tol Numeric scalar. Convergence tolerance for the EM algorithm. The algorithm
#'   stops when the Frobenius norm of changes in all parameters falls below \code{Tol}.
#'   Default: \code{1e-4}.
#' @param MaxIt Integer. Maximum number of EM iterations to perform. Default: \code{250}.
#'
#' @return A list with two elements: `Theta`, a list containing estimated mixture
#'   proportions (`Tau`), mean vectors (`Mu`), and sparse covariance matrices (`Sigma`);
#'   and `BIC`, a \eqn{k \times r}{k x r} matrix of BIC values for each component and
#'   graphical LASSO penalty parameter, where \eqn{k} is the number of components and
#'   \eqn{r} is the number of penalty values tested.
#'
#' @details
#' The algorithm proceeds in two stages:
#'
#' \strong{Stage 1: EM Algorithm}
#' The EM algorithm estimates the mixture proportions, means, and covariance matrices
#' by iteratively performing:
#' - \emph{E-step}: Compute posterior responsibilities (cluster membership probabilities)
#' - \emph{M-step}: Update parameter estimates using weighted samples
#'
#' \strong{Stage 2: Graphical LASSO}
#' After EM convergence, graphical LASSO is applied to each component's covariance matrix
#' to impose sparsity. The penalty parameter is selected via BIC for each component independently.
#'
#' @references
#' Dempster, A.P., Laird, N.M., and Rubin, D.B. (1977). Maximum likelihood from incomplete
#'   data via the EM algorithm. \emph{Journal of the Royal Statistical Society}, Series B, 39, 1-38.
#'
#' Friedman, J., Hastie, T., and Tibshirani, R. (2008). Sparse inverse covariance estimation
#'   with the graphical lasso. \emph{Biostatistics}, 9(3), 432-441.
#'
#' @examples
#' \dontrun{
#' # Generate synthetic data
#' set.seed(123)
#' Tau <- c(0.5, 0.5)
#' Mu <- list(c(0, 0), c(5, 5))
#' Sigma <- list(diag(2), diag(2))
#' X <- rmmnorm(N = 100, D = 2, Tau = Tau, Mu = Mu, Sigma = Sigma)
#'
#' # Fit emglasso
#' result <- emglasso(X, iTau = Tau, Tol = 1e-4, MaxIt = 100)
#' }
#'
#' @importFrom mvtnorm dmvnorm
#' @importFrom stats cov.wt
#' @importFrom glasso glasso
#' @export

emglasso <- function(x, iTau, Tol = 1e-4, MaxIt = 250) {
  if (!is.matrix(x)) {
    x <- as.matrix(x)
  }
  if (!is.numeric(x) || nrow(x) < 2 || ncol(x) < 1) {
    stop("x must be a numeric matrix with at least 2 rows and 1 column.")
  }
  if (!all(is.finite(x))) {
    stop("x contains non-finite values.")
  }
  if (!is.numeric(iTau) || length(iTau) < 1 || any(!is.finite(iTau)) || any(iTau < 0)) {
    stop("iTau must be a non-negative numeric vector with finite values.")
  }
  if (sum(iTau) <= 0) {
    stop("iTau must have positive sum.")
  }
  if (!is.numeric(Tol) || length(Tol) != 1 || !is.finite(Tol) || Tol <= 0) {
    stop("Tol must be a positive finite scalar.")
  }
  if (!is.numeric(MaxIt) || length(MaxIt) != 1 || !is.finite(MaxIt) || MaxIt < 1) {
    stop("MaxIt must be a positive scalar.")
  }
  MaxIt <- as.integer(MaxIt)

  Tau <- iTau / sum(iTau)
  K <- length(iTau)
  init_size <- min(nrow(x), 36)

  X <- lapply(
    1:K,
    function(k) cov.wt(x[sample(nrow(x), init_size, FALSE), ], method = "ML")
  )

  Mu <- lapply(X, "[[", "center")
  Sigma <- lapply(X, "[[", "cov")
  Theta <- list(Mu = Mu, Sigma = Sigma, Tau = Tau)
  Converged <- FALSE
  stopped_early <- FALSE
  estep_fn <- utils::getFromNamespace("estep", "emglasso")
  mstep_fn <- utils::getFromNamespace("mstep", "emglasso")
  lstep_fn <- utils::getFromNamespace("lstep", "emglasso")

  for (iter in 1:MaxIt) {
    W <- Theta$Sigma
    Tnk <- estep_fn(x, Theta, W)$Tnk

    if (!all(is.finite(Tnk))) {
      warning("Non-finite responsibilities produced in E-step; returning latest stable estimate.")
      stopped_early <- TRUE
      break
    } else if (!all(colSums(Tnk) != 0)) {
      warning("At least one component received zero total responsibility; returning latest stable estimate.")
      stopped_early <- TRUE
      break
    }

    OldTheta <- Theta
    Theta <- mstep_fn(x, Tnk)

    TauError <- norm(as.matrix(unlist(Theta$Tau) - unlist(OldTheta$Tau)))
    MuError <- norm(as.matrix(unlist(Theta$Mu) - unlist(OldTheta$Mu)))
    SigmaError <- norm(as.matrix(unlist(Theta$Sigma) - unlist(OldTheta$Sigma)))

    Converged <- TauError <= Tol && MuError <= Tol && SigmaError <= Tol

    if (Converged) {
      break
    }
  }

  if (!Converged && !stopped_early) {
    warning("EM did not converge within MaxIt iterations; returning last iterate.")
  }

  Sigma <- lstep_fn(x, Theta)
  Index <- order(Theta$Tau)

  return(list(
    Theta = list(
      Tau = Theta$Tau[Index],
      Mu = Theta$Mu[Index],
      Sigma = Sigma$W[Index]
    ),
    BIC = Sigma$BIC[Index]
  ))
}
