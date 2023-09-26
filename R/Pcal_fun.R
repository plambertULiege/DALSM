#' Generation of the penalty matrix for an additive model based on P-splines
#' @description Compute the penalty matrix associated to a vector containing fixed (non-penalized) parameters and equal-size sub-vectors of penalized B-spline parameters.
#'
#' @param nfixed the number of fixed (i.e. non-penalized) parameters.
#' @param lambda a vector of \code{p} penalty parameters where each component is associated to a sub-vector of spline parameters of length \code{J}.
#' @param Pd.x a penalty matrix of size \code{J} associated to a given sub-vector of spline parameters.
#'
#' @return A block diagonal penalty matrix of size \code{(nfixed+pJ)} given by Blockdiag(diag(0,\code{nfixed}), diag(\code{lambda}).kron.\code{Pd.x}).
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references Lambert, P. (2021). Fast Bayesian inference using Laplace approximations
#' in nonparametric double additive location-scale models with right- and
#' interval-censored data.
#' \emph{Computational Statistics and Data Analysis}, 161: 107250.
#' <doi:10.1016/j.csda.2021.107250>
#'
#' @export
#'
#' @examples
#' D = diff(diag(5),diff=2) ## Difference penalty matrix of order 2 for 5 P-spline parameters
#' P = t(D) %*% D ## Penalty matrix of order 2 for 5 B-spline parameters
#' lambda = c(100,10) ## Penalty parameters for 2 additive terms with 5 B-spline parameters each
#' Pcal.fun(3,lambda,P) ## Global penalty matrix for 3 additional non-penalized parameters
Pcal.fun = function(nfixed,lambda,Pd.x){
  ntot = nfixed + ncol(Pd.x)*length(lambda)
  Pcal = matrix(0,ntot,ntot)
  Pcal[(nfixed+1):ntot,(nfixed+1):ntot] = kronecker(diag(lambda,length(lambda)),Pd.x)
  return(Pcal)
}
