rspdmatrix=function(D,lambda,epsilon=1e-4){
#' Random Sparse Positive Definite Matrix
#' @param D number of rows
#' @param lambda sparsity of matrix
#' @param epsilon radius of eigenvalues
#' @return rMat random positive definite DxD matrix of sparsity lambda
#' @export
  rMat <- matrix(rep(NA,D*D), nrow=D) # initialize NA matrix

  pSum <- (D*(D+1))/2 # partial sums give number of elements in matrix triangle

  rMat[lower.tri(rMat,diag=T)] <- runif(pSum,0,2) 
  # fill matrix with standard uniform numbers

  ZeroIndex <- rbinom(D*D,1,lambda) # index for zeroes

  rMat[ZeroIndex==1] <- 0 # eliminate lambda% of elements

  rMat[upper.tri(rMat)] <- t(rMat)[upper.tri(rMat)] # ensure symmetric

  Eigens <- rowSums(rMat)+epsilon # diagonal elements > rowSum

  diag(rMat) <- Eigens # ensure positive definite

  return(rMat)
  }
