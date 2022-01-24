emglasso = function( x, iTau, Tol = 1e-4, MaxIt = 250 ){
#' EM with GLASSO
#'
#' @param x N obervations from D-dimensional mixed normal distribution
#' @param iTau initial estimate for Tau
#' @param Tol tolerance for convergence
#' @param MaxIt maximum iterations allowed for EM
#' @return Theta - parameter estimates for x
#' @return Rho - optimal penalty for gLASSO
#' @importFrom mvtnorm dmvnorm
#' @importFrom stats cov.wt
#' @import glasso
#' @export
#'

  #### initialization ####

  N <- length( x[,1] )	# sample Size

  Tau <- iTau # initialize mixture probabilities

  K <- length(iTau)

  X <- lapply( 1:K, function(k) cov.wt( x[sample( nrow( x ), 36, F ),], method="ML" ) )	#randomly sample K clusters

  Mu <- lapply( X, "[[", "center" )	# initial mean estimates

  Sigma <- lapply( X, "[[", "cov" )	# initial covariance estimates

  Theta <- list( Mu = Mu, Sigma = Sigma, Tau = Tau )	# initial parameter vector

  Converged <- F	# initial convergence parameter



  #### convergence loop ####

  for( iter in 1:MaxIt ){

    #### E-Step ####

    W <- Theta$Sigma

    Tnk <- estep( x, Theta, W )$Tnk	# weights

    if( !all( !is.na(Tnk) ) ){

      print( "NaN returned" )

      break

    }

    else if( !all( colSums( Tnk )!=0 ) ){

      print( "all zero column" )

      break

    }

    #### M-step ####

    OldTheta <- Theta	# previous iteration

    Theta <- mstep( x, Tnk ) # maximize likelihood function

    TauError <- norm( as.matrix( unlist( Theta$Tau )-unlist( OldTheta$Tau ) ) )	# Tau difference vector

    MuError <- norm( as.matrix( unlist( Theta$Mu ) - unlist( OldTheta$Mu ) ) )	# Mu difference vector

    SigmaError <- norm( as.matrix( unlist( Theta$Sigma ) - unlist( OldTheta$Sigma ) ) )	# Sigma difference vector

    Converged <- ( ( TauError <= Tol ) && ( MuError <= Tol ) && ( SigmaError <= Tol ) )	# TRUE if Tau, Mu, Sigma all converge

    if( Converged==TRUE ) break	# break if iterations are sufficiently close

  }



  #### L-step ####

  Sigma <- lstep( x, Theta )

  #### return list ####

  Index=order(Theta$Tau)

  return( list( Theta=list(Tau = Theta$Tau[Index], Mu = Theta$Mu[Index], +
  Sigma = Sigma$W[Index]), BIC=Sigma$BIC[Index] ) )

}
