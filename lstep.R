lstep = function( x, Theta, lTol=0.01, lMax = 3 ){
#' L-Step
#'
#' Graphical LASSO for sparse matrix estimation
#'
#' @param x N obervations from D-dimensional mixed normal distribution
#' @param Theta initial parameter estimates
#' @param lTol tolerance for matrix sparsity
#' @param lMax maximum penalty for gLASSO
#' @return W - covariance matrix pivot estimate
#' @return LHat - optimized log-likelihood values for Theta
#' @return BIC - Bayesian Info Criterion for W
#' @import glasso
#' @export
#'
#'
  N <- length( x[,1] )

  K <- length( Theta$Tau )

  D <- length( Theta$Mu[[1]] )

  Rho <- seq( 0.01, lMax, by = lTol )

  R <- length( Rho )

  BIC <- rep(list(rep(list(0),K)),R)

  for(r in seq(1,R)){

    W <- lapply( 1:K,	function(k) glasso( Theta$Sigma[[k]], Rho[r], penalize.diagonal=F ))  
    # run glasso for each rho

    D2 <-lapply( 1:K, function(k) +
    sum((W[[k]]$wi!=0) & (col(Theta$Sigma[[k]])<row(Theta$Sigma[[k]])) ) ) 
    # nonzero parameters

    BIC[[r]] <- lapply( 1:K, function(k) -2*(W[[k]]$loglik) + D2[[k]]*log(N))

  }

  BICmat <- matrix(unlist(BIC), ncol=R)

  ORI <- lapply(1:K, function(k) which.min(BICmat[k,]) )

  W <- lapply( 1:K, function( k ) glasso(Theta$Sigma[[k]],Rho[(ORI[[k]])])$w )
  # minimize BIC

  return(list(W=W, BIC=BICmat) ) # return W, BIC

}
