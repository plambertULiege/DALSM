#' Specification of knots in a cubic P-spline model
#' @description
#' Specification of knots for a cubic B-spline basis with \code{K} elements in
#' a P-spline model. The knots should support the data contained in vector \code{x}
#' and are by default assumed equidistant. Alternatively, they can be based
#' on the data quantiles. The penalty matrix of the selected penalty order
#' (3 by default) is also returned.
#'
#' @param x data that the knots should upport.
#' @param xmin desired minimum value for the knots.
#' @param xmax desired maximum value for the knots.
#' @param equid.knots logical indicating if equidistant knots are desired (Default: TRUE). If FALSE, the quantile of \code{x} are used to select the knots.
#' @param pen.order penalty order if \code{equid.knots} is \code{TRUE}.
#' @param K number of B-splines in the basis.
#'
#' @return A list containing:
#' \itemize{
#' \item{\code{xmin} : \verb{ }}{specified minimum value for the knots, except if \eqn{\min(x) < xmin}, in which case the default value \eqn{\min(x)-sd(x)} is returned.}
#' \item{\code{xmin} : \verb{ }}{specified maximum value for the knots, except if \eqn{xmax < \max(x)}, in which case the default value \eqn{\max(x)+sd(x)} is returned.}
#' \item{\code{K} : \verb{ }}{number of B-splines.}
#' \item{\code{knots} : \verb{ }}{equidistant knots on (xmin,xmax) if \code{equidistant.knots} is TRUE, based on quantiles of \code{x} otherwise.}
#' \item{\code{Pd} : \verb{ }}{\eqn{K\times K} penalty matrix of order \code{pen.order}.}
#' \item{\code{pen.order} : \verb{ }}{a reminder of the requested penalty order (Default: 3).}
#' }
#' @export
#'
#' @examples
#' x = runif(500)
#' obj = qknots(x,xmin=0,xmax=1,K=13)
#' print(obj)
qknots = function(x, xmin=NULL,xmax=NULL, equid.knots = TRUE, pen.order=3, K=25){
  if (is.null(xmin)) xmin = min(x) - sd(x)
  if (is.null(xmax)) xmax = max(x) + sd(x)
  if (xmin > min(x)) xmin = min(x) - sd(x)
  if (xmax < max(x)) xmax = max(x) + sd(x)
  cnt = 0
  if (xmin < min(x)) cnt = cnt + 1
  if (xmax > max(x)) cnt = cnt + 1
  ##
  if (equid.knots){
    knots = seq(xmin,xmax,length=K-2)
    Dd = diff(diag(K),dif=pen.order) ## Penalty order
    Pd = t(Dd)%*%Dd #+ 1e-6*diag(KK)
  } else {
    knots = unique(c(xmin,quantile(x,p=seq(0,1,length=K-2-cnt)),xmax))
    ngrid = 501
    xgrid = seq(xmin,xmax,length=ngrid)
    d2B = D2Bsplines(xgrid,knots)
    Pd = 1/ngrid * (t(d2B)%*%d2B)
    pen.order = 2 ## Penalty based on 2nd derivatives of B-splines
  }
  ##
  ans = list(xmin=xmin,xmax=xmax,knots=knots,Pd=Pd,pen.order=pen.order)
  return(ans)
}
