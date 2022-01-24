estep = function( x, Theta, W ){
#' E-Step
#'
#' Expectation Step of the EM Algorithm
#'
#' @param x N obervations from D-dimensional mixed normal distribution
#' @param Theta initial parameter estimates
#' @param W covariance matrix pivot
#' @return LogLike - log-likelihood values for Theta
#' @return Tnk - normalized loglikelihood values
#' @importFrom mvtnorm dmvnorm
#'
#'
  Theta$Tau <- Theta$Tau

  K<-length( Theta$Tau )

  LogLike <- with( Theta, do.call( cbind, lapply( 1:K, function(k) +
  Tau[k]*dmvnorm( x, Mu[[k]], W[[k]] ) ) ) )	# posteriors

  Tnk <- LogLike/rowSums( LogLike )	# normalize

  return( list(LogLike=LogLike,Tnk=Tnk) )

}
