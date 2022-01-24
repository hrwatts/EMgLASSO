mstep = function( x,Tnk ){
#' M-Step
#'
#' Maximization Step of the EM Algorithm
#'
#' @param x N obervations from D-dimensional mixed normal distribution
#' @param Tnk - normalized loglikelihood values
#' @return Theta - new parameter estimates
#' @importFrom stats cov.wt
#'

  X <- lapply( 1:ncol( Tnk ), function( k ) cov.wt( x, Tnk[,k], method="ML" ) )	# maximum likelihood estimates weighted with E-step

  Tau <-colMeans( Tnk )	# new tau estimate

  Mu <- lapply( X,"[[","center" )	# new mu estimate

  Sigma <- lapply( X,"[[","cov" )	# new sigma estimate

  Theta <- list( Mu = Mu, Sigma = Sigma, Tau=Tau )	# new parameter estimate

  return( Theta )

}
