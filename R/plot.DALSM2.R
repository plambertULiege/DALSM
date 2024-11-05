## -------------------------------------------------------
## Philippe LAMBERT (ULiege, Oct 2018, Revised in Oct 2024
## Email:  p.lambert@uliege.be
## Web: http://www.statsoc.ulg.ac.be
## -------------------------------------------------------
#' Plot visual information on a \code{DALSM.object}
#' @description Visualize the estimated additive terms and error density corresponding to a Double Additive Location-Scale Model (DALSM) object.
#'
#' @usage \method{plot}{DALSM}(x, ngrid=300, ci.level=.95, pages=0, select=NULL,
#'        fill=TRUE, pointwise=TRUE, mar=c(4,5,1,1),
#'        xlim0=NULL, ylim0=NULL, xlim1=NULL, ylim1=NULL, xlim2=NULL, ylim2=NULL,
#'        equal.ylims=TRUE, ...)
#'
#' @param x a \code{\link{DALSM.object}}.
#' @param ngrid (optional) number of points to make the plots for the fitted additive terms (default: 300).
#' @param ci.level (optional) nominal level for the plotted pointwise credible intervals (default: .95).
#' @param pages The number of pages over which to spread the output. For example, if pages=1 then all terms will be plotted on one page with the layout performed automatically. Set to 0 to have the routine leave all graphics settings as they are (default: 0).
#' @param select Allows the plot for a single model term to be selected for printing. e.g. if you just want the plot for the second smooth term set select=2. The plot of the standardized error density \eqn{f_0(z)} is provided when select=0 (default: NULL).
#' @param fill Logical indicating whether credible regions should be greyed out. (Default: TRUE).
#' @param pointwise Logical indicating whether only pointwise credible intervals should be plotted for additive terms. When FALSE, simultaneous credible regions are also provided. (Default: TRUE).
#' @param mar A numerical vector of the form c(bottom, left, top, right) which gives the number of lines of margin to be specified on the four sides of the plot. (Default: c(4,5,1,1)).
#' @param xlim0 Vector of length 2 specifying x-axis limits when plotting the estimated error density function \eqn{f_0(z)}.
#' @param ylim0 Vector of length 2 specifying y-axis limits when plotting the estimated error density function \eqn{f_0(z)}.
#' @param xlim1 Vector of length 2 specifying (common) x-axis limits when plotting the fitted additive term(s) entering the conditional mean submodel. (Default: NULL).
#' @param ylim1 Vector of length 2 specifying (common) y-axis limits when plotting the fitted additive term(s) entering the conditional mean submodel. (Default: NULL).
#' @param xlim2 Vector of length 2 specifying (common) x-axis limits when plotting the fitted additive term(s) in the log of the conditional standard deviation submodel. (Default: NULL).
#' @param ylim2 Vector of length 2 specifying (common) y-axis limits when plotting the fitted additive term(s) entering the conditional standard deviation submodel. (Default: NULL).
#' @param equal.ylims logical indicating if the same y-limits must be used when plotting the fitted additive terms from the same submodel. It can be overriden by non-NULL values for \code{ylim1} or \code{ylim2}. (Default: TRUE).
#' @param ... additional generic plotting arguments.
#'
#' @details Plot the fitted additive terms and the estimated standardized error density contained in the \code{\link{DALSM.object}} \code{x}.
#'
#' @return In addition to the plots, an invisible list containing the following is returned:
#' \itemize{
#' \item \code{J1} : number of additive terms in the location sub-model.
#' \item \code{x.loc} : a \code{ngrid} by \code{J1} matrix containing a regular grid of \code{ngrid} covariate values where the corresponding additive term in location is evaluated.
#' \item \code{f.loc} : a \code{ngrid} by \code{J1} matrix containing the \code{J1} fitted location additive terms evaluated at \code{x.loc}.
#' \item \code{se.loc} : a \code{ngrid} by \code{J1} matrix containing the the pointwise standard errors of the fitted location additive terms evaluated at \code{x.loc}.
#' \item \code{J2} : number of additive terms in the dispersion sub-model.
#' \item \code{x.disp} : a \code{ngrid} by \code{J2} matrix containing a regular grid of \code{ngrid} covariate values where the corresponding additive term in dispersion is evaluated.
#' \item \code{f.disp} : a \code{ngrid} by \code{J2} matrix containing the \code{J2} fitted dispersion additive terms evaluated at \code{x.disp}.
#' \item \code{se.disp} : a \code{ngrid} by \code{J2} matrix containing the pointwise standard errors of the fitted dispersion additive terms evaluated at \code{x.disp}.
#' }
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references Lambert, P. (2021). Fast Bayesian inference using Laplace approximations
#' in nonparametric double additive location-scale models with right- and
#' interval-censored data.
#' \emph{Computational Statistics and Data Analysis}, 161: 107250.
#' <doi:10.1016/j.csda.2021.107250>
#'
#' @examples
#' require(DALSM)
#' data(DALSM_IncomeData)
#' resp = DALSM_IncomeData[,1:2]
#' fit = DALSM(y=resp,
#'             formula1 = ~twoincomes+s(age)+s(eduyrs),
#'             formula2 = ~twoincomes+s(age)+s(eduyrs),
#'             data = DALSM_IncomeData)
#' plot(fit, select=0) ## Estimated error density
#' plot(fit, select=1:4, pages=1) ## Estimated additive terms
#'
#' @seealso \code{\link{DALSM}}, \code{\link{DALSM.object}}, \code{\link{print.DALSM}}
#'
#' @export
#'
plot.DALSM = function(x, ngrid=300, ci.level=.95, pages=0, select=NULL,
                       fill=TRUE, pointwise=TRUE, mar=c(4,5,1,1),
                       xlim0=NULL, ylim0=NULL, xlim1=NULL, ylim1=NULL, xlim2=NULL, ylim2=NULL,
                       equal.ylims=TRUE,...){
                      ## mfrow.loc=NULL, mfrow.disp=NULL,
                      ## equal.ylims=TRUE, true.loc=NULL,true.disp=NULL,
                      ## error.lim = NULL, add.residuals=FALSE, true.derr=NULL,
                      ## new.dev=TRUE, ...){
  obj.fit = x
  alpha = 1-ci.level
  z.alpha = qnorm(1-.5*alpha)
  ## ## Make sure to restore graphical user's options after function call
  ## oldpar <- par(no.readonly = TRUE) ; on.exit(par(oldpar))
  ## ## Plot fitted density
  ## if (new.dev) dev.new()
  ## par(mfrow=c(1,1))
  ##
  npar.tot = with(obj.fit,length(psi1)+length(psi2))
  E.error = with(obj.fit$derr, mean.dist) ## Mean of the current density estimate
  V.error = with(obj.fit$derr, var.dist) ## Variance of the current density estimate
  ##
  perc.obs = obj.fit$perc.obs ## Percentage exactly observed
  perc.IC = obj.fit$perc.IC ## Percentage Interval Censored
  perc.RC = obj.fit$perc.RC ## Percentage Right Censored
  ##
  ## if (is.null(error.lim)) error.lim = c(obj.fit$rmin,obj.fit$rmax)
  ## if (add.residuals){
  ##   plot.densLPS(obj.fit$derr,histRC=TRUE,breaks=25,xlim=error.lim,
  ##                                 xlab="Std.residuals",
  ##                                 main=paste("n = ",obj.fit$n," (Unc: ",perc.obs,"% ; IC: ",perc.IC,"% ; RC: ",perc.RC,"%)",sep=""))
  ##   if (!is.null(true.derr)) curve(true.derr,lwd=2,col="red",add=T)
  ## } else {
  ##   rmin = error.lim[1]
  ##   rmax = error.lim[2]
  ##   rseq = seq(.99*rmin,.99*rmax,length=501) ;
  ##   ylim = 1.05*range(c(obj.fit$derr$ddist(rseq),dnorm(rseq))) ## Make sure that dnorm() will be fully pictured
  ##   with(obj.fit$derr, curve(ddist,rmin,rmax,lwd=2,type="l",
  ##                            ylim=ylim,
  ##                            main="",
  ##                            xlab="Error",ylab="Error density",col="blue")
  ##   )
  ##   curve(dnorm,add=T,lwd=2,lty="dotted")
  ##   legend('topright',col=c("blue","black"),lty=c("solid","dotted"),lwd=2,legend=c("Estimated density","N(0,1)"),bty='n')
  ##   if (!is.null(true.derr)) with(obj.fit, curve(true.derr,rmin,rmax,lwd=2,col="red",add=T))
  ## }
  ##
  ## Compute the additive terms
  ## --------------------------
  fhat = DALSM_additive(obj.fit, ci.level=ci.level)
  ##
  ## Plot organization (extracted from mgcv::plot.gam)
  ## -------------------------------------------------
  if (is.null(select)){
      n.plots = 1 + fhat$J1 + fhat$J2 ## Number of plots: f0, F0 + additive terms
  } else {
      n.plots = length(select)
      ## if (all(select == 0)) n.plots = 1
  }
  if (pages > n.plots) pages = n.plots
  if (pages < 0) pages = 0
  if (pages != 0) {
      ppp = n.plots%/%pages
      if (n.plots%%pages != 0) {
          ppp = ppp + 1
          while (ppp * (pages - 1) >= n.plots) pages = pages - 1
      }
      c = r = trunc(sqrt(ppp))
      if (c < 1) r = c = 1
      if (c * r < ppp) c = c + 1
      if (c * r < ppp) r = r + 1
      oldpar = par(mfrow = c(r, c))
  }
  else {
      ppp = 1
      oldpar = par()
  }
  if ((pages == 0 && prod(par("mfcol")) < n.plots && dev.interactive()) ||
      pages > 1 && dev.interactive())
      ask = TRUE
  else ask = FALSE
  if (!is.null(select) && (all(select!=0))) {
      ask = FALSE
  }
  if (ask) {
      oask = devAskNewPage(TRUE)
      on.exit(devAskNewPage(oask))
  }
  ## Plot baseline f0
  ## ----------------
  if (is.null(xlim0)) xlim0 = with(obj.fit$derr, c(ymin,ymax))
  ## Plot only f0 if select=0
    if (is.null(select) | (any(select==0))){
      par(mar=mar)
      curve(obj.fit$derr$ddist(x),
            xlim=xlim0,ylim=ylim0,
            xlab=bquote(z), ylab="Error density", type="n",...)
            ## xlab="Standardized error", ylab="Error density", type="n",...)
      grid(lwd=.5,lty=1)
      curve(obj.fit$derr$ddist(x),lwd=1.5,col="blue",
            xlim=xlim0,ylim=ylim0, add=TRUE)
      curve(dnorm,add=T,lwd=2,lty="dotted")
      legend('topright',col=c("blue","black"),lty=c("solid","dotted"),lwd=2,legend=c("Estimated","N(0,1)"),bty='n')

  }
  ## Plot Additive terms
  ## -------------------
  plotAdd = function(x,y,y2=NULL,col=1,colfill=c("#CCCCCC80","#E5E5E580"),las=1,...) {
      matplot(x, y,type="n",col=col,las=las,
              xlim=xlims,ylim=ylims,xlab=xlab,ylab=ylab,
              lwd=c(1.5,1,1),lty=c(1,2,2),...)
      grid(lwd=.5,lty=1)
      if (!fill){
          if (!is.null(y2)){
              matplot(x,y2,type="l",add=TRUE,col="grey",las=las,
                      xlim=xlims,ylim=ylims,xlab=xlab,ylab=ylab,
                      lwd=c(1.5,1,1),lty=c(1,3,3),...)
          }
          matplot(x,y,type="l",add=TRUE,col=col,las=las,
                  xlim=xlims,ylim=ylims,xlab=xlab,ylab=ylab,
                  lwd=c(1.5,1,1),lty=c(1,2,2),...)
      } else {
          if (!is.null(y2)){
              plotRegion(x,y2,add=TRUE,col=col,colfill=colfill[2],las=las,
                         xlim=xlims,ylim=ylims,
                         xlab=xlab,ylab=ylab,lwd=2,...)
          }
          plotRegion(x,y,add=TRUE,col=col,colfill=colfill[1],las=las,
                     xlim=xlims,ylim=ylims,
                     xlab=xlab,ylab=ylab,lwd=2,...)
      }
  }
  ##
  if (fhat$J1 > 0){
      xlims = ylims = NULL
      if (equal.ylims){
          if (pointwise){
              ylims =  range(lapply(fhat$f.loc.grid, function(x) range(x$y.mat)))
          } else {
              ylims =  range(lapply(fhat$f.loc.grid, function(x) range(x$y.mat2)))
          }
      }
      if (!is.null(ylim1)) ylims = ylim1
      if (!is.null(xlim1)) xlims = xlim1
      for (j in 1:fhat$J1){
          if ((is.null(select)) || (any(select == j))){
              par(mar=mar)
              xlab = names(fhat$f.loc.grid)[j]
              ylab = bquote('f'[.(j)]^{~mu}*(.(xlab)))
              if (!pointwise) {
                  with(fhat$f.loc.grid[[j]], plotAdd(x,y.mat,y.mat2,...))
              } else {
                  with(fhat$f.loc.grid[[j]], plotAdd(x,y.mat,...))
              }
              ## if (obj.fit$regr1$has.ref[j]){
              ##     xt = obj.fit$regr1$ref.values[j]
              ##     rug(xt,ticksize=-.015,lwd=1)
              ##     rug(xt,ticksize=.015,lwd=1)
              ## }
          }
      }
  }
  ##
  if (fhat$J2 > 0){
      xlims = ylims = NULL
      if (equal.ylims){
          if (pointwise){
              ylims =  range(lapply(fhat$f.disp.grid, function(x) range(x$y.mat)))
          } else {
              ylims =  range(lapply(fhat$f.disp.grid, function(x) range(x$y.mat2)))
          }
      }
      if (!is.null(ylim2)) ylims = ylim2
      if (!is.null(xlim2)) xlims = xlim2
      for (j in 1:fhat$J2){
          if ((is.null(select)) || (any(select == fhat$J1+j))){
              par(mar=mar)
              xlab = names(fhat$f.disp.grid)[j]
              ylab = bquote('f'[.(j)]^{~sigma}*(.(xlab)))
              if (!pointwise) {
                  with(fhat$f.disp.grid[[j]], plotAdd(x,y.mat,y.mat2,...))
              } else {
                  with(fhat$f.disp.grid[[j]], plotAdd(x,y.mat,...))
              }
              ## if (obj$regr2$has.ref[j]){
              ##     xt = obj$regr2$ref.values[j]
              ##     rug(xt,ticksize=-.015,lwd=1)
              ##     rug(xt,ticksize=.015,lwd=1)
              ## }
          }
      }
  }
  ##
  if (pages > 0) par(oldpar)
  return(invisible(fhat))
}
