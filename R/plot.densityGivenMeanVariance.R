#' Plot the error density of a fitted DALSM regression model
#' @description Plot the error density corresponding to the fit of a DALSM regression model. The histogram of the observed residuals can be superposed to it.
#' @param x a <densityGivenMeanVariance.object> (typically the $derr element of a DALSM.object list)
#' @param xlim interval of values where the density should be plotted
#' @param breaks (Optional) breaks for the histogram of the observed residuals
#' @param histRC Logical (Default: FALSE) indicating whether the histogram of the right-censored residuals should be highlighted
#' @param xlab Optional label for the x-axis
#' @param ylab Optional label for the y-axis
#' @param main Plot main title
#' @param ... Optional additional plot parameters
#'
#' @return No returned value (just plots)
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references Lambert, P. (2021). Fast Bayesian inference using Laplace approximations
#' in nonparametric double additive location-scale models with right- and
#' interval-censored data.
#' \emph{Computational Statistics and Data Analysis}, 161: 107250.
#' \url{https://doi.org/10.1016/j.csda.2021.107250}
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
#' plot(fit$derr)  ## Plot the estimated error density
#' print(fit$derr) ## ... and provide some descriptive elements on it

plot.densityGivenMeanVariance = function(x,xlim=range(fit$bins),breaks=NULL,histRC=FALSE,xlab="",ylab="Density",main="",...){
  fit = x
  dens.grid = fit$ddist(fit$ugrid) ## Fitted density at bin midpoints
  ##
  if (is.null(breaks)) breaks="Sturges"
  temp = with(fit, hist(rep(ugrid,apply((C/apply(C,1,sum)),2,sum)),breaks=breaks,plot=FALSE))
  temp = with(fit, hist(rep(ugrid,apply((C/apply(C,1,sum)),2,sum)),breaks=breaks,freq=FALSE,
                        xlim=xlim,ylim=1.05*c(0,max(dens.grid,temp$density)),
                        xlab=xlab,main=main,col="grey"))
  n.RC = sum(fit$event==0)
  if (histRC & n.RC>0) temp = with(fit, hist(rep(ugrid,apply((C[event==0,]/apply(C[event==0,],1,sum)),2,sum)),breaks=breaks,freq=FALSE,col="gray87",add=T))
  with(fit, curve(ddist,lwd=2,add=T,n=501))
}
