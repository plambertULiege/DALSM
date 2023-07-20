## Function to generate Penalty matrix for additive terms
#' Compute the penalty matrix associated to a vector containing fixed (non-penalized) parameters and equal-size sub-vectors of penalized spline parameters
#'
#' @param nfixed the number of fixed (i.e. non-penalized) parameters
#' @param lambda a vector of \code{p} penalty parameters where each component is associated to a sub-vector of spline parameters of length \code{J}
#' @param Pd.x a penalty matrix of size \code{J} associated to a given sub-vector of spline parameters
#'
#' @return A block diagonal penalty matrix of size \code{(nfixed+pJ)} given by Blockdiag(diag(0,\code{nfixed}), diag(\code{lambda}).kron.\code{Pd.x})
#' @export
#'
#' @examples
Pcal.fun = function(nfixed,lambda,Pd.x){
  ntot = nfixed + ncol(Pd.x)*length(lambda)
  Pcal = matrix(0,ntot,ntot)
  Pcal[(nfixed+1):ntot,(nfixed+1):ntot] = kronecker(diag(lambda,length(lambda)),Pd.x)
  return(Pcal)
}
