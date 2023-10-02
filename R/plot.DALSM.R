## -----------------------------------
## Philippe LAMBERT (ULiege, Oct 2018
## Email:  p.lambert@uliege.be
## Web: http://www.statsoc.ulg.ac.be
## -----------------------------------
#' Plot visual information on a \code{DALSM.object}
#' @description Visualize the estimated additive terms and error density corresponding to a Double Additive Location-Scale Model (DALSM) object.
#'
#' @usage \method{plot}{DALSM}(x,
#'        mfrow.loc=NULL, mfrow.disp=NULL,
#'        nx=101, equal.ylims=TRUE, true.loc=NULL,true.disp=NULL, ci.level=NULL,
#'        error.lim = NULL, add.residuals=FALSE, true.derr=NULL, new.dev=TRUE, ...)
#'
#' @param x a \code{\link{DALSM.object}}.
#' @param mfrow.loc (optional) window layout to plot the additive terms for location.
#' @param mfrow.disp (optional) window layout to plot the additive terms for dispersion.
#' @param nx (optional) number of points to make the plots for the fitted additive terms (default: 101).
#' @param equal.ylims logical indicating if the same y-limits must be used when plotting the fitted additive terms (default: TRUE).
#' @param true.loc (optional) list of functions to be superposed to the corresponding estimated additive terms in the location submodel (default: NULL).
#' @param true.disp (optional) list of functions to be superposed to the corresponding estimated additive terms in the dispersion submodel (default: NULL).
#' @param ci.level (optional) nominal level for the plotted pointwise credible intervals (default: x$ci.level).
#' @param error.lim (optional) plotting interval for the estimated standardized error density in the DALSM model (default: support of the fitted standardized error density).
#' @param add.residuals logical requesting to add the (possibly censored) standardized residuals to the plot of the fitted standardized error density (default: FALSE).
#' @param true.derr (optional) density function to superpose to the estimated standardized error density when plotting (default: NULL).
#' @param new.dev (optional) logical indicating whether a new plotting device must be opened for each graph (default: TRUE).
#' @param ... additional generic plotting arguments.
#'
#' @details Plot the fitted additive terms and the estimated standardized error density contained in the \code{\link{DALSM.object}} \code{x}.
#'
#' @return In addition to the plots, an invisible list containing the following is returned:
#' \itemize{
#' \item{\code{J1} : \verb{ }}{number of additive terms in the location sub-model.}
#' \item{\code{x.loc} : \verb{ }}{a \code{nx} by \code{J1} matrix containing a regular grid of \code{nx} covariate values where the corresponding additive term in location is evaluated.}
#' \item{\code{f.loc} : \verb{ }}{a \code{nx} by \code{J1} matrix containing the \code{J1} fitted location additive terms evaluated at \code{x.loc}.}
#' \item{\code{se.loc} : \verb{ }}{a \code{nx} by \code{J1} matrix containing the the pointwise standard errors of the fitted location additive terms evaluated at \code{x.loc}.}
#' \item{\code{J2} : \verb{ }}{number of additive terms in the dispersion sub-model.}
#' \item{\code{x.disp} : \verb{ }}{a \code{nx} by \code{J2} matrix containing a regular grid of \code{nx} covariate values where the corresponding additive term in dispersion is evaluated.}
#' \item{\code{f.disp} : \verb{ }}{a \code{nx} by \code{J2} matrix containing the \code{J2} fitted dispersion additive terms evaluated at \code{x.disp}.}
#' \item{\code{se.disp} : \verb{ }}{a \code{nx} by \code{J2} matrix containing the pointwise standard errors of the fitted dispersion additive terms evaluated at \code{x.disp}.}
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
#' plot(fit)
#'
#' @seealso \code{\link{DALSM}}, \code{\link{DALSM.object}}, \code{\link{print.DALSM}}
#'
#' @export
#'
plot.DALSM = function(x,
                      mfrow.loc=NULL, mfrow.disp=NULL,
                      nx = 101, equal.ylims=TRUE, true.loc=NULL,true.disp=NULL, ci.level=NULL,
                      error.lim = NULL, add.residuals=FALSE, true.derr=NULL,
                      new.dev=TRUE, ...){
  obj.fit = x
  if (is.null(ci.level)) ci.level = x$ci.level
  alpha = 1-ci.level
  z.alpha = qnorm(1-.5*alpha)
  ## Make sure to restore graphical user's options after function call
  oldpar <- par(no.readonly = TRUE) ; on.exit(par(oldpar))
  ## Plot fitted density
  if (new.dev) dev.new()
  par(mfrow=c(1,1))
  ##
  npar.tot = with(obj.fit,length(psi1)+length(psi2))
  E.error = with(obj.fit$derr, mean.dist) ## Mean of the current density estimate
  V.error = with(obj.fit$derr, var.dist) ## Variance of the current density estimate
  ##
  perc.obs = round(100*obj.fit$derr$n.uncensored/obj.fit$n,1) ## Percentage exactly observed
  perc.IC = round(100*sum(obj.fit$derr$n.IC)/obj.fit$n,1) ## Percentage Interval Censored
  perc.RC = round(100*sum(1-obj.fit$derr$event)/obj.fit$n,1) ## Percentage Right Censored
  ##
  if (is.null(error.lim)) error.lim = c(obj.fit$rmin,obj.fit$rmax)
  if (add.residuals){
    plot.densLPS(obj.fit$derr,histRC=TRUE,breaks=25,xlim=error.lim,
                                  xlab="Std.residuals",
                                  main=paste("n = ",obj.fit$n," (Unc: ",perc.obs,"% ; IC: ",perc.IC,"% ; RC: ",perc.RC,"%)",sep=""))
    if (!is.null(true.derr)) curve(true.derr,lwd=2,col="red",add=T)
  } else {
    rmin = error.lim[1]
    rmax = error.lim[2]
    rseq = seq(.99*rmin,.99*rmax,length=501) ;
    ylim = 1.05*range(c(obj.fit$derr$ddist(rseq),dnorm(rseq))) ## Make sure that dnorm() will be fully pictured
    with(obj.fit$derr, curve(ddist,rmin,rmax,lwd=2,type="l",
                             ylim=ylim,
                             main="",
                             ## main=paste("n = ",obj.fit$n," (Unc: ",perc.obs,"% ; IC: ",perc.IC,"% ; RC: ",perc.RC,"%)",sep=""),
                             xlab="Error",ylab="Error density",col="blue")
    )
    curve(dnorm,add=T,lwd=2,lty="dotted")
    legend('topright',col=c("blue","black"),lty=c("solid","dotted"),lwd=2,legend=c("Estimated density","N(0,1)"),bty='n')
    if (!is.null(true.derr)) with(obj.fit, curve(true.derr,rmin,rmax,lwd=2,col="red",add=T))
  }
  ## Plot fitted additive terms
  obj.add = DALSM_additive(obj.fit)
  Bx1 = with(obj.fit,regr1$Bx) ; Bx2 = with(obj.fit,regr2$Bx)
  pen.order.1 = with(obj.fit,regr1$pen.order) ; pen.order.2 = with(obj.fit,regr2$pen.order)
  knots.1 = with(obj.fit,regr1$knots.x) ; knots.2 = with(obj.fit,regr2$knots.x)
  nfixed1 = with(obj.fit,regr1$nfixed) ; nfixed2 = with(obj.fit,regr2$nfixed)
  K1 = with(obj.fit,regr1$K) ; K2 = with(obj.fit,regr2$K)
  J1 = with(obj.fit,regr1$J) ; J2 = with(obj.fit,regr2$J)
  psi1.cur = obj.fit$psi1 ; psi2.cur = obj.fit$psi2
  ##
  ## ... for location
  ## x.loc = f.loc = se.loc = NULL
  if (J1 > 0){
    if (is.null(mfrow.loc)){
      if (J1==1) mfrow.loc = c(1,1)
      else mfrow.loc = c(ceiling(J1/2),2)
    }
    if (new.dev) dev.new()
    par(mfrow=mfrow.loc,mar=c(4,5,1,1))
    Cst1 = coverage.loc = rep(0,J1)
    ## temp = matrix(nrow=nx,ncol=J1) ; rownames(temp) = 1:nx
    ## x.loc = f.loc = se.loc = temp
    ## addloc.lab = x$regr1$additive.lab
    ## colnames(x.loc) = addloc.lab ## paste("x.loc.",1:J1,sep="")
    ## colnames(f.loc) = paste("f(",addloc.lab,")",sep="") ## paste("f.loc.",1:J1,sep="")
    ## colnames(se.loc) = paste("se.",addloc.lab,sep="") ## paste("se.loc.",1:J1,sep="")
    for (j in 1:J1){
      ## xlow = min(knots.1[[j]]) ; xup = max(knots.1[[j]])
      ## x.loc[,j] =
      xx = obj.add$f.loc.grid[[j]]$x ## seq(xlow,xup,length=nx)
      # BB.1 = centeredBasis.gen(xx,knots.1[[j]],cm=Bx1[[j]]$cm,pen.order=pen.order.1)$B
      # idx = nfixed1 + (j-1)*K1 + (1:K1)
      # S.loc = obj.fit$Cov.psi1
      # S11 = S.loc[idx,idx]
      ## f.loc[,j] = fj.hat = BB.1%*%psi1.cur[idx]
      ## se.loc[,j] = fj.se = sqrt(pmax(0,rowSums((BB.1%*%S11)*BB.1)))
      ## fj.min = fj.hat-z.alpha*fj.se ; fj.max = fj.hat+z.alpha*fj.se
      f.mat = obj.add$f.loc.grid[[j]]$y.mat
      fj.min = f.mat[,2] ; fj.max = f.mat[,3]
      ##
      # ## Global credible region (essai)
      ## omega = sqrt(pmax(0,rowSums(BB.1*(BB.1%*%S11)))) / sqrt(qchisq(1-alpha,K1))
      ## Theta.up = psi1.cur[idx] + S11 %*% t((1/omega)*(BB.1))  ## K1 x nx matrix
      ## Theta.low = psi1.cur[idx] - S11 %*% t((1/omega)*(BB.1)) ## K1 x nx matrix
      ## cred.low = colSums(Theta.low*t(BB.1))
      ## cred.up = colSums(Theta.up*t(BB.1))
      xlab = names(obj.add$f.loc)[j]
      ## xlab = colnames(x.loc)[j] ## colnames(obj.fit$regr1$X)[j]
      ## if (equal.ylims) ylims = c(min(cred.low),max(cred.up))
      ## else  ylims = c(min(f.mat),max(f.mat))
      ylims = c(min(f.mat),max(f.mat))
      plotRegion(xx,f.mat, ##cbind(fj.hat,fj.min,fj.max),
                 xlab=xlab,ylab=bquote('f'[.(j)]^{~mu}*(.(xlab))), ##colnames(f.loc)[j], ## paste("f1.",j,"(x)",sep=""),
                 ylim=ylims,
                 main="") ## bquote('f'[.(j)]^{~mu}*(.(xlab[j]))))
                 ## paste("Loc:  ",colnames(f.loc)[j],sep="")) ## main=paste("Loc:  f1.",j,"(x)",sep=""))
      with(obj.fit,rug(jitter(regr1$X[,j])))
      # if (equal.ylims){
      #   lines(xx,cred.low,col="blue",lwd=2)
      #   lines(xx,cred.up,col="blue",lwd=2)
      # }
      if (!is.null(true.loc)){
        Cst1[j] = integrate(true.loc[[j]],0,1)$val
        true.loc.j = true.loc[[j]](xx)-Cst1[j]
        lines(xx,true.loc.j,lwd=2,col="red")
        coverage.loc[j] = mean((true.loc.j > fj.min) * (true.loc.j < fj.max))
      }
      abline(a=0,b=0,lty="dotted")
    }
  }
  ## ... for dispersion
  x.disp = f.disp = se.disp = NULL
  if (J2 > 0){
    if (is.null(mfrow.disp)){
      if (J2==1) mfrow.disp = c(1,1)
      else mfrow.disp = c(ceiling(J2/2),2)
    }
    if (new.dev) dev.new()
    par(mfrow=mfrow.disp,mar=c(4,5,1,1))
    Cst2 = coverage.disp = rep(0,J2)
    ## temp = matrix(nrow=nx,ncol=J2) ; rownames(temp) = 1:nx
    ## x.disp = f.disp = se.disp = temp
    ## adddisp.lab = x$regr2$additive.lab
    ## colnames(x.disp) = adddisp.lab ## paste("x.disp.",1:J2,sep="")
    ## colnames(f.disp) = paste("f(",adddisp.lab,")",sep="") ## paste("f.disp.",1:J2,sep="")
    ## colnames(se.disp) = paste("se.",adddisp.lab,sep="") ## paste("se.disp.",1:J2,sep="")
    for (j in 1:J2){
      ## xlow = min(knots.2[[j]]) ; xup = max(knots.2[[j]])
      xx = obj.add$f.disp.grid[[j]]$x ## seq(xlow,xup,length=nx)
      ## x.disp[,j] = xx
      # BB.2 = centeredBasis.gen(xx,knots.2[[j]],cm=Bx2[[j]]$cm,pen.order=pen.order.2)$B
      # idx = nfixed2 + (j-1)*K2 + (1:K2)
      # S.disp = obj.fit$Cov.psi2
      # S11 = S.disp[idx,idx]
      ## f.disp[,j] = fj.hat = BB.2%*%psi2.cur[idx]
      ## se.disp[,j] = fj.se = sqrt(pmax(0,rowSums((BB.2%*%S11)*BB.2)))
      ## fj.min = fj.hat-z.alpha*fj.se ; fj.max = fj.hat+z.alpha*fj.se
      f.mat = obj.add$f.disp.grid[[j]]$y.mat
      fj.min = f.mat[,2] ; fj.max = f.mat[,3]
      ## Global credible region (essai)
      # omega = sqrt(pmax(0,rowSums(BB.2*(BB.2%*%S11)))) / sqrt(qchisq(1-alpha,K1))
      # Theta.up = psi2.cur[idx] + S11 %*% t((1/omega)*(BB.2))  ## K2 x length(xx) matrix
      # Theta.low = psi2.cur[idx] - S11 %*% t((1/omega)*(BB.2)) ## K2 x length(xx) matrix
      # cred.low = colSums(Theta.low*t(BB.2))
      # cred.up = colSums(Theta.up*t(BB.2))
      xlab = names(obj.add$f.disp)[j]
      ## xlab = colnames(x.disp)[j] ## colnames(obj.fit$regr2$X)[j]
      # if (equal.ylims) ylims = c(min(cred.low),max(cred.up))
      # else  ylims = c(min(f.mat),max(f.mat))
      ylims = c(min(f.mat),max(f.mat))
      plotRegion(xx,f.mat, ## cbind(fj.hat,fj.min,fj.max),
                 xlab=xlab, ylab=bquote('f'[.(j)]^{~sigma}*(.(xlab))),
                 ##ylab=colnames(f.disp)[j], ##paste("f2.",j,"(x)",sep=""),
                 ylim=ylims,
                 main="") ##paste("Disp:  ",colnames(f.disp)[j],sep="")) ## f2.",j,"(x)",sep=""))
      with(obj.fit,rug(jitter(regr2$X[,j])))
      # if (equal.ylims){
      #   lines(xx,cred.low,col="blue",lwd=2)
      #   lines(xx,cred.up,col="blue",lwd=2)
      # }
      if (!is.null(true.disp)){
        Cst2[j] = integrate(true.disp[[j]],0,1)$val
        true.disp.j = true.disp[[j]](xx)-Cst2[j]
        lines(xx,true.disp.j,lwd=2,col="red")
        coverage.disp[j] = mean((true.disp.j > fj.min) * (true.disp.j < fj.max))
      }
      abline(a=0,b=0,lty="dotted")
    }
  }
  ## Return observed coverages if the true additive terms are known
  ans = obj.add
  return(invisible(ans))
}
