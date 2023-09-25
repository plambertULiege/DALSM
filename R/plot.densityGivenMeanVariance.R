#' Plot the density estimate in a \code{densityGivenMeanVariance.object}
#' @description Plot the density estimate obtained by \code{densityGivenMeanVariance} from censored data with given mean and variance.
#'
#' @usage \method{plot}{densityGivenMeanVariance}(x,
#'        xlim=range(fit$bins),breaks=NULL,histRC=FALSE,
#'        xlab="",ylab="Density",main="",...)
#'
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
#' @seealso \code{\link{densityGivenMeanVariance.object}}, \code{\link{print.densityGivenMeanVariance}}, \code{\link{plot.densityGivenMeanVariance}}, \code{\link{Dens1d.object}}, \code{\link{Dens1d}}.
#'
#' @export
#'
#' @examples
#' require(DALSM)
#' ## Example 1: density estimation from IC data
#' n = 500 ## Sample size
#' x = 3 + rgamma(n,10,2) ## Exact generated data
#' width = runif(n,1,3) ## Width of the IC data (mean width = 2)
#' w = runif(n) ## Positioning of the exact data within the interval
#' xmat = cbind(x-w*width,x+(1-w)*width) ## Generated IC data
#' head(xmat)
#' obj.data = Dens1d(xmat,ymin=0) ## Prepare the data for estimation
#' ## Density estimation with fixed mean and variance
#' obj = densityGivenMeanVariance(obj.data,Mean0=3+10/2,Var0=10/4)
#' plot(obj) ## Histogram of thepseudo-data with the density estimate
#' curve(dgamma(x-3,10,2), ## ... compared to the true density (in red)
#'       add=TRUE,col="red",lwd=2,lty=2)
#' legend("topright",col=c("black","red","grey"),lwd=c(2,2,20),lty=c(1,2,1),
#'        legend=c("Fitted density","True density","Pseudo-data"),bty="n")
#' print(obj) ## ... with summary statistics
#'
#' ## Example 2: estimation of the error density in a DALSM model
#' data(DALSM_IncomeData)
#' resp = DALSM_IncomeData[,1:2]
#' fit = DALSM(y=resp,
#'             formula1 = ~twoincomes+s(age)+s(eduyrs),
#'             formula2 = ~twoincomes+s(age)+s(eduyrs),
#'             data = DALSM_IncomeData)
#' plot(fit$derr)  ## Plot the estimated error density
#' legend("topright",col=c("black","grey"),lwd=c(2,20),lty=c(1,1),
#'        legend=c("Estimated error density","Pseudo-residuals"),bty="n")
#' print(fit$derr) ## ... and provide summary statistics for it

plot.densityGivenMeanVariance = function(x,xlim=range(fit$bins),breaks=NULL,histRC=FALSE,xlab="",ylab="Density",main="",...){
  fit = x
  dens.grid = fit$ddist(fit$ugrid) ## Fitted density at bin midpoints
  ##
  if (is.null(breaks)) breaks="Sturges"
  temp = with(fit, hist(rep(ugrid,ncol(C)*apply((C/apply(C,1,sum)),2,sum)),breaks=breaks,plot=FALSE))
  temp = with(fit, hist(rep(ugrid,ncol(C)*apply((C/apply(C,1,sum)),2,sum)),breaks=breaks,freq=FALSE,
                        xlim=xlim,ylim=1.05*c(0,max(dens.grid,temp$density)),
                        xlab=xlab,main=main,col="grey"))
  n.RC = sum(fit$event==0)
  if (histRC & n.RC>0) temp = with(fit, hist(rep(ugrid,apply((C[event==0,]/apply(C[event==0,],1,sum)),2,sum)),breaks=breaks,freq=FALSE,col="gray87",add=T))
  with(fit, curve(ddist,lwd=2,add=T,n=501))
}
