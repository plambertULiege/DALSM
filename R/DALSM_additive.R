## ------------------------------------
## Philippe LAMBERT (ULiege, June 2019)
## Email:  p.lambert@uliege.be
## Web: http://www.statsoc.ulg.ac.be
## ------------------------------------
#' Extraction of the estimated additive terms in a \code{\link{DALSM.object}}
#'
#' @description Extract the estimated additive terms with their standard errors from a \code{\link{DALSM.object}} resulting from the fit of a \code{\link{DALSM}} model.
#' In addition to the estimated functions in the location and dispersion submodels, their values on a regular grid covering the observed covariate values
#' are reported together with credible intervals. The mean effective coverage of these pointwise credible intervals
#' for the additive terms with respect to given (optional) reference functions (such as the ones for the 'true' additive terms used
#' to generate data in a simulation study) can also be computed.
#' @usage DALSM_additive(obj.DALSM, ngrid=101,
#'        true.loc=NULL, true.disp=NULL, ci.level=NULL,
#'        verbose=FALSE)
#'
#' @param obj.DALSM a \code{\link{DALSM.object}}
#' @param ngrid (optional) grid size of covariate values where the additive terms are calculated (default: 101).
#' @param true.loc (optional) list of functions containing the 'true' additive terms in the location sub-model.
#' @param true.disp (optional) list of functions containing the 'true' additive terms in the dispersion sub-model.
#' @param ci.level (optional) level of credible intervals.
#' @param verbose logical indicating whether the computed corverages should be printed out (default: TRUE).
#'
#' @return It returns an invisible list containing:
#' \itemize{
#' \item{\code{J1} : \verb{ }}{number of additive terms in the location sub-model.}
#' \item{\code{labels.loc} : \verb{ }}{labels of the additive terms in the location sub-model.}
#' \item{\code{f.loc.grid} : \verb{ }}{list of length \code{J1} with, for each additive term, a list of length 3 with elements 'x': a vector of \code{ngrid} values for the covariate ; 'y.mat': a matrix with 3 columns (est,low,up) giving the additive term and its pointwise credible region ; se: the standard error of the additive term on the x-grid.}
#' \item{\code{f.loc} : \verb{ }}{a list of length \code{J1} with, for each additive term <x>, a list with f.loc$x: a function computing the additive term f.loc(x) for a given covariate value 'x' ; attributes(f.loc$x): support, label, range.}
#' \item{\code{se.loc} : \verb{ }}{a list of length \code{J1} with, for each additive term <x>, a list with se.loc$x: a function computing the s.e. of f(x) for a given covariate value 'x' ; attributes(se.loc$x): support, label, range.}
#' \item{\code{coverage.loc} : \verb{ }}{if \code{true.loc} is provided: a vector of length \code{J1} giving the average effective coverage of pointwise credible intervals for each of the additive terms in the location sub-model.}
#' \item{\code{J2} : \verb{ }}{number of additive terms in the dispersion sub-model.}
#' \item{\code{labels.disp} : \verb{ }}{labels of the additive terms in the dispersion sub-model.}
#' \item{\code{f.disp.grid} : \verb{ }}{list of length \code{J2} with, for each additive term, a list of length 3 with elements 'x': a vector of \code{ngrid} values for the covariate ; 'y.mat': a matrix with 3 columns (est,low,up) giving the additive term and its pointwise credible region ; se: the standard error of the additive term on the x-grid.}
#' \item{\code{f.disp} : \verb{ }}{a list of length \code{J2} with, for each additive term <x>, a list with f.disp$x: a function computing the additive term f.disp(x) for a given covariate value 'x' ; attributes(f.disp$x): support, label, range.}
#' \item{\code{se.disp} : \verb{ }}{a list of length \code{J2} with, for each additive term <x>, a list with se.disp$x: a function computing the s.e. of f(x) for a given covariate value 'x' ; attributes(se.disp$x): support, label, range.}
#' \item{\code{coverage.disp} : \verb{ }}{if <true.disp> is provided: a vector of length \code{J2} giving the average effective coverage of pointwise credible intervals for each of the additive terms in the dispersion sub-model.}
#' }
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references Lambert, P. (2021). Fast Bayesian inference using Laplace approximations
#' in nonparametric double additive location-scale models with right- and
#' interval-censored data.
#' \emph{Computational Statistics and Data Analysis}, 161: 107250.
#' <doi:10.1016/j.csda.2021.107250>
#'
#' @seealso \code{\link{DALSM.object}}, \code{\link{DALSM}}, \code{\link{print.DALSM}}, \code{\link{plot.DALSM}}.
#'
#' @export
#' @examples
#' require(DALSM)
#' data(DALSM_IncomeData)
#' resp = DALSM_IncomeData[,1:2]
#' fit = DALSM(y=resp,
#'             formula1 = ~twoincomes+s(age)+s(eduyrs),
#'             formula2 = ~twoincomes+s(age)+s(eduyrs),
#'             data = DALSM_IncomeData)
#' obj = DALSM_additive(fit)
#' str(obj)
#'
#' ## Visualize the estimated additive terms for Age
#' ## ... in the location submodel
#' with(obj$f.loc.grid$age, matplot(x,y.mat,
#'                                  xlab="Age",ylab="f.loc(Age)",
#'                                  type="l",col=1,lty=c(1,2,2)))
#' ## ... and in the dispersion submodel
#' with(obj$f.disp.grid$age, matplot(x,y.mat,
#'                                  xlab="Age",ylab="f.disp(Age)",
#'                                  type="l",col=1,lty=c(1,2,2)))
#' ## Also report their values for selected age values
#' obj$f.loc$age(c(30,40,50)) ; obj$f.disp$age(c(30,40,50))
#' ## ... together with their standard errors
#' obj$se.loc$age(c(30,40,50)) ; obj$se.disp$age(c(30,40,50))
DALSM_additive <- function(obj.DALSM,
                          ngrid=101,
                          true.loc=NULL,true.disp=NULL,
                          ci.level=NULL,
                          verbose=FALSE){
  if (is.null(ci.level)) ci.level = obj.DALSM$ci.level
  alpha = 1-ci.level
  z.alpha = qnorm(1-.5*alpha)
  ##
  npar.tot = with(obj.DALSM,length(psi1)+length(psi2))
  E.error = with(obj.DALSM$derr, mean.dist) ## Mean of the current density estimate
  V.error = with(obj.DALSM$derr, var.dist) ## Variance of the current density estimate
  ## Some stats on censoring (exactly observed, RC, IC)
  perc.obs = round(100*obj.DALSM$derr$n.uncensored/obj.DALSM$n,1) ## Percentage exactly observed
  perc.IC = round(100*sum(obj.DALSM$derr$n.IC)/obj.DALSM$n,1) ## Percentage Interval Censored
  perc.RC = round(100*sum(1-obj.DALSM$derr$event)/obj.DALSM$n,1) ## Percentage Right Censored
  ## Stats on fitted additive terms
  Bx1 = with(obj.DALSM,regr1$Bx) ; Bx2 = with(obj.DALSM,regr2$Bx)
  pen.order.1 = with(obj.DALSM,regr1$pen.order) ; pen.order.2 = with(obj.DALSM,regr2$pen.order)
  knots.1 = with(obj.DALSM,regr1$knots.x) ; knots.2 = with(obj.DALSM,regr2$knots.x)
  nfixed1 = with(obj.DALSM,regr1$nfixed) ; nfixed2 = with(obj.DALSM,regr2$nfixed)
  K1 = with(obj.DALSM,regr1$K) ; K2 = with(obj.DALSM,regr2$K)
  J1 = with(obj.DALSM,regr1$J) ; J2 = with(obj.DALSM,regr2$J)
  psi1.cur = obj.DALSM$psi1 ; psi2.cur = obj.DALSM$psi2
  ## ... for location
  f.loc.grid = f.loc = se.loc = list()
  ## f.loc = se.loc = NULL
  if (J1 > 0){
    Cst1 = coverage.loc = rep(0,J1)
    add.lab = labels.loc = obj.DALSM$regr1$additive.lab
    ## temp = matrix(nrow=ngrid,ncol=J1) ; rownames(temp) = 1:ngrid
    ## x.loc = f.loc = se.loc = temp
    ## addloc.lab = obj.DALSM$regr1$additive.lab
    ## colnames(x.loc) = addloc.lab ## paste("x.loc.",1:J1,sep="")
    ## colnames(f.loc) = paste("f(",addloc.lab,")",sep="") ## paste("f.loc.",1:J1,sep="")
    ## colnames(se.loc) = paste("se.",addloc.lab,sep="") ## paste("se.loc.",1:J1,sep="")
    names(coverage.loc) = paste("f(",add.lab,")",sep="") ## colnames(f.loc)
    for (j in 1:J1){
      xlow = min(knots.1[[j]]) ; xup = max(knots.1[[j]])
      x.grid = seq(xlow,xup,length=ngrid)
      ## x.loc[,j] = x.grid
      BB.1 = centeredBasis.gen(x.grid,knots.1[[j]],cm=Bx1[[j]]$cm,pen.order=pen.order.1)$B
      idx = nfixed1 + (j-1)*K1 + (1:K1)
      S.loc = obj.DALSM$Cov.psi1
      S11 = S.loc[idx,idx]
      ##
      fj.grid = BB.1%*%psi1.cur[idx]
      fj.se = sqrt(pmax(0,rowSums((BB.1%*%S11)*BB.1)))
      fj.min = fj.grid-z.alpha*fj.se ; fj.max = fj.grid+z.alpha*fj.se
      ## f.loc[,j] = fj.grid = BB.1%*%psi1.cur[idx]
      ## se.loc[,j] = fj.se = sqrt(pmax(0,rowSums((BB.1%*%S11)*BB.1)))
      ## fj.min = fj.grid-z.alpha*fj.se ; fj.max = fj.grid+z.alpha*fj.se
      ##
      f.loc.grid[[add.lab[j]]]$x = x.grid
      mat = cbind(est=fj.grid, low=fj.min, up=fj.max) ; colnames(mat) = c("est","low","up")
      f.loc.grid[[add.lab[j]]]$y.mat = mat ## cbind(est=fj.grid, low=fj.min, up=fj.max)
      f.loc.grid[[add.lab[j]]]$se = fj.se
      ##
      f.loc[[add.lab[j]]]  = splinefun(x.grid, fj.grid)
      se.loc[[add.lab[j]]] = splinefun(x.grid, fj.se)
      attr(f.loc[[add.lab[j]]],"support") = attr(se.loc[[add.lab[j]]],"support") = c(xlow,xup)
      attr(f.loc[[add.lab[j]]],"label") = attr(se.loc[[add.lab[j]]],"label") = add.lab[j]
      attr(f.loc[[add.lab[j]]],"range") = range(fj.grid)
      attr(se.loc[[add.lab[j]]],"range") = range(fj.se)
      ##
      ## Global credible region (essai)
      omega = sqrt(pmax(0,rowSums(BB.1*(BB.1%*%S11)))) / sqrt(qchisq(1-alpha,K1))
      Theta.up = psi1.cur[idx] + S11 %*% t((1/omega)*(BB.1)) ## K1 x length(x.grid) matrix
      Theta.low = psi1.cur[idx] - S11 %*% t((1/omega)*(BB.1)) ## K1 x length(x.grid) matrix
      cred.low = colSums(Theta.low*t(BB.1))
      cred.up = colSums(Theta.up*t(BB.1))
      if (!is.null(true.loc)){
        Cst1[j] = integrate(true.loc[[j]],xlow,xup)$val
        true.loc.j = true.loc[[j]](x.grid)-Cst1[j]
        coverage.loc[j] = mean((true.loc.j > fj.min) * (true.loc.j < fj.max))
      }
    }
  }
  ## ... for dispersion
  f.disp.grid = f.disp = se.disp = list()
  ## f.disp = se.disp = NULL
  if (J2 > 0){
    Cst2 = coverage.disp = rep(0,J2)
    add.lab = labels.disp = obj.DALSM$regr2$additive.lab
    ## temp = matrix(nrow=ngrid,ncol=J2) ; rownames(temp) = 1:ngrid
    ## x.disp = f.disp = se.disp = temp
    ## adddisp.lab = obj.DALSM$regr2$additive.lab
    ## colnames(x.disp) = adddisp.lab ## paste("x.disp.",1:J2,sep="")
    ## colnames(f.disp) = paste("f(",adddisp.lab,")",sep="") ## paste("f.disp.",1:J2,sep="")
    ## colnames(se.disp) = paste("se.",adddisp.lab,sep="") ## paste("se.disp.",1:J2,sep="")
    names(coverage.disp) = paste("f(",add.lab,")",sep="") ## colnames(f.disp)
    for (j in 1:J2){
      xlow = min(knots.2[[j]]) ; xup = max(knots.2[[j]])
      x.grid = seq(xlow,xup,length=ngrid)
      BB.2 = centeredBasis.gen(x.grid,knots.2[[j]],cm=Bx2[[j]]$cm,pen.order=pen.order.2)$B
      idx = nfixed2 + (j-1)*K2 + (1:K2)
      S.disp = obj.DALSM$Cov.psi2
      S11 = S.disp[idx,idx]
      ##
      fj.grid = BB.2%*%psi2.cur[idx]
      fj.se = sqrt(pmax(0,rowSums((BB.2%*%S11)*BB.2)))
      fj.min = fj.grid-z.alpha*fj.se ; fj.max = fj.grid+z.alpha*fj.se
      ## f.disp[,j] = fj.grid = BB.2%*%psi2.cur[idx]
      ## se.disp[,j] = fj.se = sqrt(pmax(0,rowSums((BB.2%*%S11)*BB.2)))
      ## fj.min = fj.grid-z.alpha*fj.se ; fj.max = fj.grid+z.alpha*fj.se
      ##
      f.disp.grid[[add.lab[j]]]$x = x.grid
      mat = cbind(est=fj.grid, low=fj.min, up=fj.max) ; colnames(mat) = c("est","low","up")
      f.disp.grid[[add.lab[j]]]$y.mat = mat ## cbind(est=fj.grid, low=fj.min, up=fj.max)
      f.disp.grid[[add.lab[j]]]$se = fj.se
      ##
      f.disp[[add.lab[j]]]  = splinefun(x.grid, fj.grid)
      se.disp[[add.lab[j]]] = splinefun(x.grid, fj.se)
      attr(f.disp[[add.lab[j]]],"support") = attr(se.disp[[add.lab[j]]],"support") = c(xlow,xup)
      attr(f.disp[[add.lab[j]]],"label") = attr(se.disp[[add.lab[j]]],"label") = add.lab[j]
      attr(f.disp[[add.lab[j]]],"range") = range(fj.grid)
      attr(se.disp[[add.lab[j]]],"range") = range(fj.se)
      ##
      ## Global credible region (essai)
      omega = sqrt(pmax(0,rowSums(BB.2*(BB.2%*%S11)))) / sqrt(qchisq(1-alpha,K1))
      Theta.up = psi2.cur[idx] + S11 %*% t((1/omega)*(BB.2)) ## K1 x length(x.grid) matrix
      Theta.low = psi2.cur[idx] - S11 %*% t((1/omega)*(BB.2)) ## K1 x length(x.grid) matrix
      cred.low = colSums(Theta.low*t(BB.2))
      cred.up = colSums(Theta.up*t(BB.2))
      if (!is.null(true.disp)){
        Cst2[j] = integrate(true.disp[[j]],xlow,xup)$val
        true.disp.j = true.disp[[j]](x.grid)-Cst2[j]
        coverage.disp[j] = mean((true.disp.j > fj.min) * (true.disp.j < fj.max))
      }
    }
  }
  ## Return observed coverages if the true additive terms are known
  ans = list(J1=J1) ## Number of additive terms for location

  if (J1 > 0){
    ans = append(ans,
            list(
             labels.loc = labels.loc, ## Labels of the additive terms for location
             f.loc.grid=f.loc.grid, ## Fitted location additive terms on a grid with pointwise standard errors and credible intervals
             f.loc=f.loc, ## Estimated functions for the location additive terms
             se.loc=se.loc)) ## Estimated standard error functions for the location additive terms
  }
  ans = append(ans,
               list(J2=J2)) ## Number of additive terms for dispersion
  if (J2 > 0){
    ans = append(ans,
            list(
              labels.disp = labels.disp, ## Labels of the additive terms for dispersion
              f.disp.grid=f.disp.grid, ## Fitted dispersion additive terms on a grid with pointwise standard errors and credible intervals
              f.disp=f.disp, ## Estimated functions for the dispersion additive terms
              se.disp=se.disp)) ## Estimated standard error functions for the dispersionadditive terms
  }
  ans$coverage.loc = ans$coverage.disp = NULL
  if (!is.null(true.loc)&(J1>0)){
        ans$coverage.loc = coverage.loc
  }
  if (!is.null(true.disp)&(J2>0)){
      ans$labels.disp = labels.disp
      ans$coverage.disp = coverage.disp
  }
  if (!is.null(ans$coverage.loc) | !is.null(ans$coverage.disp)){
      if (verbose){
          cat("Observed coverages for additive terms (with nominal coverage ",1-alpha,")\n",sep="")
          print(lapply(ans[c("coverage.loc","coverage.disp")], function(x) round(x,3)))
      }
  }
  return(invisible(ans))
}
