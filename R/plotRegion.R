## -----------------------------------
## Philippe LAMBERT (ULiege, Oct 2018
## Email:  p.lambert@uliege.be
## Web: http://www.statsoc.ulg.ac.be
## -----------------------------------
#' Plot a credible region for a curve together with its point estimates
#' @description Plot a credible region (in grey) for a curve together with its point estimates.
#'
#' @usage plotRegion(x,mat,
#'        add=FALSE, xlim=range(x), ylim=range(mat),
#'        lwd=2, xlab="", ylab="", main="", ...)
#' @param x vector of values where the curve is evaluated.
#' @param mat cbind(f.hat,f.low,f.up) is a matrix containing the point estimates <f.hat> and the liming values <f.low> and <f.up> for the credible region.
#' @param add Logical indicating if the plot must be added to the active plot (Default: FALSE).
#' @param xlim range of <x> values for which the plot should be provided.
#' @param ylim range of curve values that should be considered for the plot.
#' @param lwd line width for the plot (Default: 2).
#' @param xlab x-label.
#' @param ylab y-label.
#' @param main main title.
#' @param ... optional plotting parameters.
#'
#' @return No returned value (just a plot)
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references Lambert, P. (2021). Fast Bayesian inference using Laplace approximations
#' in nonparametric double additive location-scale models with right- and
#' interval-censored data.
#' \emph{Computational Statistics and Data Analysis}, 161: 107250.
#' <doi:10.1016/j.csda.2021.107250>
#'
#' @seealso \code{\link{DALSM}}, \code{\link{DALSM.object}}.
#'
#' @export
#'
#' @examples
#' require(DALSM)
#' data(DALSM_IncomeData)
#' resp = DALSM_IncomeData[,1:2]
#' fit = DALSM(y=resp,
#'             formula1 = ~twoincomes+s(age)+s(eduyrs),
#'             formula2 = ~twoincomes+s(age)+s(eduyrs),
#'             data = DALSM_IncomeData)
#' obj = DALSM_additive(fit)
#' ## par(mfrow=c(1,2),mar=c(4,5,1,1))
#' with(obj$f.loc.grid$age, plotRegion(x, y.mat,
#'      xlab="age", ylab=expression('f'[1]^{~mu}*(age))))
#' with(obj$f.disp.grid$age, plotRegion(x, y.mat,
#'      xlab="age", ylab=expression('f'[1]^{~sigma}*(age))))
plotRegion = function(x,mat,add=FALSE,xlim=range(x),ylim=range(mat),
                      lwd=2,xlab="",ylab="",main="",...){
  f = mat[,1] ; f.low = mat[,2] ; f.up = mat[,3]
  if (add==FALSE) plot(x,f,type="n",ylim=ylim,xlim=xlim,
                       lwd=lwd,xlab=xlab,ylab=ylab,main=main,...)
  polygon(c(x,rev(x)),c(f.low,rev(f.up)),col="grey",border=F)
  lines(x,f,lwd=lwd)
}
