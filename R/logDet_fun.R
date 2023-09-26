#' Log-determinant of a positive-definite matrix
#'
#' @param x positive definite matrix.
#'
#' @return log of det(x).
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
#' A = matrix(1:4,ncol=2)
#' logDet.fun(A)
#'
logDet.fun = function(x){
  if(any(is.nan(x))) return(NA)
  if(any(!is.finite(x))) return(NA)
  eival = svd(x)$d
  idx = which(eival > 0)
  return(sum(log(eival[idx])))
}
