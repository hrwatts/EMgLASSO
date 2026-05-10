#' L-Step: Sparse Covariance Estimation via Graphical LASSO
#'
#' Applies graphical LASSO regularization to estimate sparse precision (inverse covariance)
#' matrices. The regularization parameter is selected separately for each component using BIC.
#'
#' @param x An \eqn{n \times d}{n x d} matrix of observations.
#' @param Theta A list of current parameter estimates with elements `Tau`
#'   (mixture proportions), `Mu` (component mean vectors), and `Sigma`
#'   (component covariance matrices).
#' @param lTol Numeric. Step size for the penalty grid. Default: \code{0.01}.
#' @param lMax Numeric. Maximum penalty value. Default: \code{3}.
#'
#' @return A list with `W`, the estimated sparse covariance matrices for each component,
#'   and `BIC`, a \eqn{k \times r}{k x r} matrix of BIC values over the penalty grid.
#'
#' @details
#' For each component and penalty value \eqn{\rho}, the Bayesian Information Criterion is:
#' \deqn{BIC_\rho = -2 \log L(\hat{W}_\rho) + d_\rho \log(n)}
#' where \eqn{d_\rho} is the number of non-zero off-diagonal elements in the estimated
#' precision matrix and \eqn{n} is the sample size.
#'
#' The penalty that minimizes BIC is selected for each component independently.
#'
#' @references
#' Friedman, J., Hastie, T., and Tibshirani, R. (2008). Sparse inverse covariance estimation
#'   with the graphical lasso. \emph{Biostatistics}, 9(3), 432-441.
#'
#' @importFrom glasso glasso
#' @keywords internal

lstep <- function(x, Theta, lTol = 0.01, lMax = 3) {
  N <- nrow(x)
  K <- length(Theta$Tau)
  Rho <- seq(0.01, lMax, by = lTol)
  R <- length(Rho)
  BIC <- matrix(NA_real_, nrow = K, ncol = R)

  for (r in seq_len(R)) {
    W <- lapply(
      seq_len(K),
      function(k) glasso(Theta$Sigma[[k]], Rho[r], penalize.diagonal = FALSE)
    )

    D2 <- lapply(
      seq_len(K),
      function(k) sum((W[[k]]$wi != 0) & (col(Theta$Sigma[[k]]) < row(Theta$Sigma[[k]])))
    )

    BIC[, r] <- vapply(
      seq_len(K),
      function(k) -2 * (W[[k]]$loglik) + D2[[k]] * log(N),
      numeric(1)
    )
  }

  ORI <- vapply(
    seq_len(K),
    function(k) {
      finite_idx <- which(is.finite(BIC[k, ]))
      if (length(finite_idx) == 0) {
        return(1L)
      }
      finite_idx[which.min(BIC[k, finite_idx])]
    },
    integer(1)
  )
  W <- lapply(
    seq_len(K),
    function(k) glasso(Theta$Sigma[[k]], Rho[ORI[k]], penalize.diagonal = FALSE)$w
  )

  return(list(W = W, BIC = BIC))
}
