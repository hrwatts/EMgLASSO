rmmnorm=function(N,D,Tau,Mu,Sigma){
#' Random MMVN Sample Generator
#'
#' Generates observations from a random mixed multivarate normal (MMVN) distribution
#'
#' @param N Number of observations
#' @param D Dimension of MMVN distribution
#' @param Tau Mixture ratios
#' @param Mu List of mean vectors for each component in Tau
#' @param Sigma List of covariance matrices for each component in Tau
#' @return X - N observations from Tau-mixed D-dimensional MMVN distribution
#' @importFrom MASS mvrnorm
#' @export
#'
  K <- length(Tau) # number of mixture components

  Normal <- rep(list(0),K) # initialize list of observations

  for(ii in seq(1,K)){

    Normal[[ii]]<-mvrnorm(n = N%*%Tau[ii], Mu[[ii]], Sigma[[ii]]); 
    # use mvrnorm for each mixture component

  }

  X <- do.call("rbind",Normal)

}
