#' Fit a double additive location-scale model (DALSM) with a flexible error distribution
#' @description Fit a location-scale regression model with a flexible error distribution and
#' additive terms in location (=mean) and dispersion (= log(sd)) using Laplace P-splines
#' from potentially right- and interval-censored response data.
#' @usage DALSM(y, formula1,formula2, w, data,
#'        K1=10, K2=10, pen.order1=2, pen.order2=2,
#'        a.pen=1, b.pen=1e-4, lambda1.min=1, lambda2.min=1,
#'        phi.0=NULL,psi1.0=NULL,psi2.0=NULL,
#'        psi.method=c("LM","NR","BFGS","NM"),
#'        lambda0.0=NULL,lambda1.0=NULL,lambda2.0=NULL,
#'        lambda.method=c("LPS","LPS2","none"), logscale=FALSE,
#'        density.method=c("LPS","Schall"),
#'        REML=FALSE, diag.only=TRUE, Normality=FALSE, sandwich=TRUE,
#'        joint.computation=FALSE,
#'        K.error=20, rmin=NULL,rmax=NULL,
#'        ci.level=.95,
#'        grad.tol=1e-2, RDM.tol=1e-4, fun.tol=1e-3,
#'        iterlim=50,verbose=FALSE)
#'
#' @param y n-vector of responses or (nx2)-matrix/data.frame when interval-censoring (with the 2 elements giving the interval bounds) or right-censoring (when the element in the 2nd column equals +Inf).
#' @param formula1 model formula for location (i.e. for the conditional mean).
#' @param formula2 model formula for dispersion (i.e. for the log of the conditional standard deviation).
#' @param w vector of length \code{n} containing (non-negative) weights (default: rep(1,n)).
#' @param data data frame containing the model covariates.
#' @param K1 (optional) number of B-splines to describe a given additive term in the location submodel (default: 10).
#' @param K2 (optional) number of B-splines to describe a given additive term in the dispersion submodel (default: 10).
#' @param pen.order1 (optional) penalty order for the additive terms in the location submodel (default: 2).
#' @param pen.order2 (optional) penalty order for the additive terms in the dispersion submodel (default: 2).
#' @param a.pen (optional) prior on penalty parameter \eqn{\tau} is Gamma(a.pen=1,b.pen) (for additive terms) (default: 1).
#' @param b.pen (optional) prior on penalty parameter \eqn{\tau} is Gamma(a.pen,b.pen=1e-4) (for additive terms) (default: 1e-4).
#' @param lambda1.min (optional) minimal value for the penalty parameters in the additive model for location (default: 1.0).
#' @param lambda2.min (optional) minimal value for penalty parameters in the additive model for log-dispersion (default: 1.0).
#' @param phi.0 (optional) initial values for the spline parameters in the log-hazard of the standardized error distribution.
#' @param psi1.0 (optional) initial values for the location submodel parameters.
#' @param psi2.0 (optional) initial values for the dispersion submodel parameters.
#' @param psi.method (optional) algorithm used for the computation of the conditional posterior mode of the regression and splines parameters. Possible choices are Levenberg-Marquardt ("LM"), Newton-Raphson ("NR"), Broyden–Fletcher–Goldfarb–Shanno ("BFGS") and Nelder-Mead ("NM") (default: "LM").
#' @param lambda0.0 (optional) initial value for the penalty parameter used during the estimation of the error density.
#' @param lambda1.0 (optional) initial value for the J1 penalty parameters of the additive terms in the location submodel.
#' @param lambda2.0 (optional) initial value for the J2 penalty parameters of the additive terms in the dispersion submodel.
#' @param lambda.method (optional) algorithm used for the computation of the marginal posterior mode of the penalty parameters tuning the smoothness of the estimated additive terms. Possible choices are Laplace P-splines ("LPS" or "LPS2") or "none" if fixed quantities given by initial values (by default set at 100) are desired.
#' @param logscale this parameter is only relevant when "LPS2" is chosen for \code{lambda.method}. Setting \code{logscale=TRUE} selects the penalty parameter on the log scale, maximizing the marginal posterior of its log-transformed value (default: FALSE).
#' @param density.method algorithm used for the selection of the penalty parameter during density estimation. Possible choices are Laplace- P-splines ("LPS") or Fellner-Schall ("Schall").
#' @param REML (optional) logical indicating if a REML correction is desired to estimate the dispersion parameters (default: FALSE).
#' @param diag.only (optional) logical indicating if only the diagonal of the Hessian needs to be corrected during REML (default: TRUE).
#' @param Normality (optional) logical indicating if Normality is assumed for the error term (default: FALSE).
#' @param sandwich (optional) logical indicating if sandwich variance estimators are needed for the regression parameters for a NP error distribution (when Normality=FALSE) ; it is forced to be TRUE when Normality=TRUE.
#' @param joint.computation indicates whether the computation of the effective degrees of freedom for the additive terms should be made jointly for the location and dispersion additive terms (default: FALSE).
#' @param K.error (optional) number of B-splines to approximate the log of the error hazard on (rmin,rmax) (default: 20).
#' @param rmin (optional) minimum value for the support of the standardized error distribution.
#' @param rmax (optional) maximum value for the support of the standardized error distribution.
#' @param ci.level (optional) nominal level for the reported credible intervals (default: .95)
#' @param grad.tol tolerance threshold for the L2-norm of the gradient when monitoring convergence. (default: 1e-2).
#' @param RDM.tol tolerance thershold for the Relative Damping Measure (= RDM) when monitoring convergence (default: 1e-4).
#' @param fun.tol tolerance threshold for variations in the maximized function during the final iterations of posterior mode computation and convergence monitoring (default: 1e-3).
#' @param iterlim (optional) maximum number of iterations (after which the algorithm is interrupted with a non-convergence diagnostic) (default: 50).
#' @param verbose (optional) logical indicating whether estimation step details should be displayed (default: FALSE).
#'
#' @return An object of class \code{DASLM} containing several components from the fit,
#' see \code{\link{DALSM.object}} for details. A summary can be printed using \code{\link{print.DALSM}}
#' or plotted using \code{\link{plot.DALSM}}.
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references Lambert, P. (2021). Fast Bayesian inference using Laplace approximations
#' in nonparametric double additive location-scale models with right- and
#' interval-censored data.
#' \emph{Computational Statistics and Data Analysis}, 161: 107250.
#' <doi:10.1016/j.csda.2021.107250>
#'
#' @seealso \code{\link{DALSM.object}}, \code{\link{print.DALSM}}, \code{\link{plot.DALSM}}.
#'
#' @examples
#' require(DALSM)
#' data(DALSM_IncomeData)
#' resp = DALSM_IncomeData[,1:2]
#' fit = DALSM(y=resp,
#'             formula1 = ~twoincomes+s(age)+s(eduyrs),
#'             formula2 = ~twoincomes+s(age)+s(eduyrs),
#'             data = DALSM_IncomeData)
#' print(fit)
#' plot(fit, select=0) ## Estimated error density
#' plot(fit, select=1:4, pages=1) ## Estimated additive terms
#'
#' @export

DALSM <- function(y, formula1,formula2, w, data,
                  K1=10, K2=10,pen.order1=2, pen.order2=2,
                  a.pen=1, b.pen=1e-4,
                  lambda1.min=1, lambda2.min=1,
                  phi.0=NULL,psi1.0=NULL,psi2.0=NULL,
                  psi.method=c("LM","NR","BFGS","NM"), ## "LevMarq","LevMarq1","nloptr","nlm"),
                  lambda0.0=NULL,lambda1.0=NULL,lambda2.0=NULL,
                  lambda.method=c("LPS","LPS2","none"), ## Penalty selection method for additive terms
                  logscale=FALSE, ## Maximize p(log(lambda)|D) when TRUE,  p(lambda|D) otherwise
                  density.method=c("LPS","Schall"),
                  REML=FALSE, diag.only=TRUE, Normality=FALSE, sandwich=TRUE,
                  joint.computation=FALSE,
                  K.error=20, rmin=NULL,rmax=NULL,
                  ci.level=.95,
                  grad.tol=1e-2, RDM.tol=1e-4, fun.tol=1e-3,
                  iterlim=50,verbose=FALSE){
  cl <- match.call()
  lambda.method = match.arg(lambda.method)
  ##
  density.method = match.arg(density.method)
  psi.method = match.arg(psi.method)
  if (missing(y)) {message("Missing <y> with the response vector (when exactly observed data) or matrix (when RC or IC data) !") ; return(NULL)}
  if (missing(formula1)) {message("Missing model formula <formula1> for location !") ; return(NULL)}
  if (missing(formula2)) {message("Missing model formula <formula2> for dispersion !") ; return(NULL)}
  if (missing(data)) {message("Missing data frame <data> with the covariate values !") ; return(NULL)}
  fixed.support = (!is.null(rmin)) & (!is.null(rmax)) ## Indicates if a fixed pre-specified support for the error term is desired
  alpha = 1-ci.level ## Significance level
  ## Sample size n, ylow, yup, ymid
  ## ------------------------------
  if (is.data.frame(y)) y = as.matrix(y)
  if (is.matrix(y)){ ## Response matrix (--> IC or RC data involved)
    ylow = y[,1] ; yup = y[,2]
  } else { ## Response vector (--> NO IC data)
    ylow = yup = y
    y = cbind(ylow,yup) # Turn the response vector to a matrix with 2 identical columns
  }
  n = nrow(y)
  if (n != nrow(data)){message("<y> and <data> should share the same number of units !") ; return(NULL)}
  if (missing(w)) w = rep(1,n) ## Default 1.0 weight for all units
  if (!all(w>=0)) stop("Weights in vector <w> should all be non negative !")
  if (length(w) != n) stop("Weight vector <w> should be of length <n> !")
  event = ifelse(is.finite(yup),1,0) ## Event = 1 if non-RC, 0 otherwise
  resp = y ; resp[,2] = ifelse(is.finite(yup),yup,ylow) ## Working (nx2) response matrix with identical columns if non-IC (in particular when RC)
  ymid = rowMeans(resp) ## apply(resp,1,mean) ## Midpoint of IC response, reported exactly observed or right-censored response otherwise
  is.IC = (ylow != yup) & (is.finite(yup)) ## Indicated whether Interval Censored (with, then, event=1)
  n.IC = sum(w * is.IC) ## Weighted number of interval-censored data
  ##
  ## Numbers of RC & IC data
  is.uncensored = (ylow==yup) & (event==1)
  n.uncensored = sum(w * is.uncensored) ## Weighted number of uncensored
  is.RC = (event==0)
  n.RC = sum(w * is.RC) ## Weighted number of right-censored
  sw = sum(w) ## Total weighted sample size
  ##
  ## Index of units with a reported event
  obs = (event==1)
  ##
  ## Regression models for location(formula1) & dispersion (formula2)
  ## ----------------------------------------------------------------
  regr1 = DesignFormula(formula1, data=data, K=K1, pen.order=pen.order1, n=n)
  regr2 = DesignFormula(formula2, data=data, K=K2, pen.order=pen.order2, n=n)
  ##
  ## Extract key stuff from regr1 & regr2
  J1 = regr1$J ; J2 = regr2$J
  K1 = regr1$K ; K2 = regr2$K
  Bx1 = regr1$Bx ; Bx2 = regr2$Bx
  pen.order.1 = regr1$pen.order ; pen.order.2 = regr2$pen.order
  nfixed1 = regr1$nfixed ; nfixed2 = regr2$nfixed
  lambda1.lab = regr1$lambda.lab ; lambda2.lab = regr2$lambda.lab
  Pd1.x = regr1$Pd.x ; Pd2.x = regr2$Pd.x
  ## Pd1.x = Pd1.x + diag(1e-4,ncol(Pd1.x)); Pd2.x = Pd2.x + diag(1e-4,ncol(Pd2.x))
  Z1 = regr1$Z ; Z2 = regr2$Z
  loc.lab = colnames(regr1$Xcal) ; disp.lab = colnames(regr2$Xcal) ## Labels of the regression parms
  addloc.lab = regr1$additive.lab ; adddisp.lab = regr2$additive.lab ## Labels of the additive terms
  # addloc.lab = colnames(regr1$X) ; adddisp.lab = colnames(regr2$X) ## Labels of the additive terms
  Xcal.1 = regr1$Xcal ; Xcal.2 = regr2$Xcal
  q1 = ncol(Xcal.1) ; q2 = ncol(Xcal.2) ## Total number of regression and spline parameters
  z.alpha = qnorm(1-.5*alpha)
  ##
  ## Initial values
  ## ... for the spline parameters defining the standardized error density
  if (is.null(lambda0.0)) lambda0.0 = eval(formals(densityLPS)$tau0) ## Default value for the Penalty parameter
  lambda0.cur = lambda0.0
  if (!is.null(phi.0)){
    K.error = length(phi.0)
    phi.cur =  phi.0
    if (Normality) cat("Normality=TRUE ignored as values for <phi.0> are provided !\n")
  }
  ## ... for location stuff
  if (J1 > 0){
    Bcal.1 = regr1$Bcal
    if (is.null(lambda1.0)) lambda1.0 = rep(100,J1)
    names(lambda1.0) = lambda1.lab
    lambda1.cur = lambda1.0 ; P1.cur = Pcal.fun(nfixed1,lambda1.cur,Pd1.x)
    psi1.cur = c(solve(t(Xcal.1[obs,,drop=FALSE])%*%(w[obs]*Xcal.1[obs,,drop=FALSE]) + P1.cur, t(w[obs]*Xcal.1[obs,,drop=FALSE])%*%ymid[obs]))
  } else {
    psi1.cur = c(solve(t(Xcal.1[obs,,drop=FALSE])%*%(w[obs]*Xcal.1[obs,,drop=FALSE]), t(w[obs]*Xcal.1[obs,,drop=FALSE])%*%ymid[obs]))
    lambda1.cur = NULL
  }
  names(psi1.cur) = loc.lab
  ##
  ## ... for dispersion stuff
  fun <- function(y,w){ ## Weighted mean and sd
    sw = sum(w)
    mean.y = sum((w/sw) * y) ; sd.y = sqrt(sum((w/sw) * (y - mean.y)^2) * (sw / max(1,(sw-1))))
    return(list(mu=mean.y,sd=sd.y))
  }
  sd0 = 1.1 * fun(y=ymid[obs]-c(Xcal.1[obs,,drop=FALSE]%*%psi1.cur), w=w[obs])$sd
  if (J2 > 0){
    Bcal.2 = regr2$Bcal
    if (is.null(lambda2.0)) lambda2.0 = rep(100,J2) ## PHL: 10 instead of 100
    names(lambda2.0) = lambda2.lab
    lambda2.cur = lambda2.0 ; P2.cur = Pcal.fun(nfixed2,lambda2.cur,Pd2.x)
    psi2.cur = c(log(sd0),rep(0,ncol(Z2)-1),rep(0,ncol(Bcal.2))) ## PHL
  } else {
    psi2.cur = c(log(sd0),rep(0,ncol(Z2)-1))
    lambda2.cur = NULL
  }
  names(psi2.cur) = disp.lab
  ##
  ## Substitute initial values for <psi1> and <psi2> if provided
  if (!is.null(psi1.0)) psi1.cur[1:length(psi1.0)] = psi1.0
  if (!is.null(psi2.0)) psi2.cur[1:length(psi2.0)] = psi2.0
  ##
  psi.cur = c(psi1.cur,psi2.cur)
  ##
  ## Effective dimensions / quadratic forms
  if (J1 > 0) ED1.cur = quad1.cur = 0*lambda1.cur
  if (J2 > 0) ED2.cur = quad2.cur = 0*lambda2.cur
  ##
  ## ####################
  ## ESTIMATION PROCEDURE
  ## ####################
  res.fun <- function(resp,mu,sd){
    res = (resp-mu)/sd
    if ((!fixed.support) & (!is.null(obj1d))){
      rmin = min(obj1d$knots) ; rmax = max(obj1d$knots)
    }
    res[res<=rmin] = .99*rmin ; res[res>=rmax] = .99*rmax ## If necessary, "Huber-like" M-estimation strategy
    return(res)
  }
  ##
  res.psi.fun <- function(resp,psi){
    psi1 = psi[1:q1] ; psi2 = psi[-(1:q1)]
    mu = c(Xcal.1 %*% psi1)
    eta = c(Xcal.2 %*% psi2) ; sd = exp(eta)
    res = res.fun(resp,mu,sd)
    return(res)
  }
  ## EDF.fun: Effective dimensions of additive terms
  ## -------
  parms1 = with(regr1, list(q=ncol(Xcal), nfixed=nfixed, K=K, J=J, addnames=colnames(X), Pd.x=Pd.x))
  parms2 = with(regr2, list(q=ncol(Xcal), nfixed=nfixed, K=K, J=J, addnames=colnames(X), Pd.x=Pd.x))
  ##
  EDF.fun <- function(Hes0.psi, lambda1, lambda2, parms1, parms2, joint.computation=TRUE){
      ## parms1 = with(regr1, list(q=ncol(Xcal), nfixed=nfixed, K=K, J=J, addnames=colnames(X), Pd.x=Pd.x))
      ## parms2 = with(regr2, list(q=ncol(Xcal), nfixed=nfixed, K=K, J=J, addnames=colnames(X), Pd.x=Pd.x))
      J1 = parms1$J ; J2 = parms2$J
      q1 = parms1$q ; q2 = parms2$q
      ED1 = rep(0,J1) ; ED2 = rep(0,J2)
      idx1 = 1:q1 ; idx2 = parms1$q + (1:q2)
      P.cur = diag(0,q1+q2)
      ##
      grp1 = with(parms1, rep(0,nfixed))
      if (J1 > 0){
          P1.cur = with(parms1, Pcal.fun(nfixed, lambda1, Pd.x))
          P.cur[idx1,idx1] = P1.cur
          grp1 = with(parms1, c(grp1, as.numeric(gl(J,K))))
      }
      ##
      grp2 = with(parms2, rep(0,nfixed))
      if (J2 > 0){
          P2.cur = with(parms2, Pcal.fun(nfixed, lambda2, Pd.x))
          P.cur[idx2,idx2] = P2.cur
          grp2 = with(parms2, c(grp2, as.numeric(gl(J,K))))
      }
      Hes.psi = Hes0.psi - P.cur
      ##
      ans = list()
      if (joint.computation){
          grp = c(grp1, max(grp1) + grp2)
          ##
          cond_number <- kappa(Hes.psi)  # Condition number
          threshold <- 1e12              # Threshold for singularity
          if (cond_number > threshold) {
              Cov.psi <- MASS::ginv(-Hes.psi)
          } else {
              Cov.psi <- solve(-Hes.psi)
          }
          ##
          temp = diag(Cov.psi %*% (-Hes0.psi))
          temp1 = temp[1:q1] ; temp2 = temp[-(1:q1)]
      } else {
          Cov.psi1 = solve(-Hes.psi[idx1,idx1])
          temp1 = diag(Cov.psi1 %*% (-Hes0.psi[idx1,idx1]))
          ##
          Cov.psi2 = solve(-Hes.psi[idx2,idx2])
          temp2 = diag(Cov.psi2 %*% (-Hes0.psi[idx2,idx2]))
      }
      if (J1 > 0) {
          ED1 = tapply(temp1, grp1, sum)
          ED1 = ED1[-1] # remove fixed part
          names(ED1) = parms1$addnames
          ans$ED1 = ED1
      }
      ans$ED1.tot = parms1$nfixed + ifelse(J1==0, 0, sum(ED1))
      if (J2 > 0) {
          ED2 = tapply(temp2, grp2, sum)
          ED2 = ED2[-1] # remove fixed part
          names(ED2) = parms2$addnames
          ans$ED2 = ED2
      }
      ans$ED2.tot = parms2$nfixed + ifelse(J2==0, 0, sum(ED2))
      ans$ED.tot = with(ans, ED1.tot + ED2.tot)
      return(ans)
  }
  ## End EDF.fun
  ##
  ## -----------------------------------------------------------------
  ## Key function:  Computation of Log-posterior, Score, Hessian, etc.
  ## -----------------------------------------------------------------
  ff <- function(psi, lambda1,lambda2,
                gradient=TRUE, hessian=TRUE,REML,diag.only=TRUE){
    if (REML) hessian=TRUE  ## Computation of the Hessian required for REML !!
    if (hessian) gradient = TRUE ## Gradient required to get Hessian
    ##
    psi1 = psi1.cur = psi[1:q1] ; psi2 = psi2.cur = psi[-(1:q1)]
    lambda1.cur = lambda1 ; lambda2.cur = lambda2
    ## Function to interrupt computation when found necessary
    stopthis <- function(){
      ans = list(lpost=NA,U.psi=NA,Hes.psi=NA)
      return(ans)
    }
    ## Evaluate the standardized residuals for given <psi>
    ## ---------------------------------------------------
    mu.cur = c(Xcal.1 %*% psi1.cur)
    sd.cur = exp(c(Xcal.2 %*% psi2.cur))
    res = res.fun(resp,mu.cur,sd.cur)
    ##
    # if (any(is.na(res)) | !all((res >= rmin) & (res <= rmax))) return(NA)
    ##
    ## Preliminary computation
    ## -----------------------
    dB1 = dB2 = D1Bsplines(res[,1], fit1d$knots) ; ## 1st derivative of the B-spline basis
    dB2[is.IC,] = D1Bsplines(res[is.IC,2], fit1d$knots) ## (also at RH limit of IC data)
    d2B = D2Bsplines(res[,1], fit1d$knots) ## 2nd derivative of the B-spline basis
    ## (Only required for non-IC data for which res[,1]=res[,2] !!)
    h1.cur = h2.cur = fit1d$hdist(res[,1]) ; h2.cur[is.IC] = fit1d$hdist(res[is.IC,2]) ## Hazard value
    H1.cur = H2.cur = fit1d$Hdist(res[,1]) ; H2.cur[is.IC] = fit1d$Hdist(res[is.IC,2]) ## Cumulative hazard value
    f1.cur = f2.cur = fit1d$ddist(res[,1]) ; f2.cur[is.IC] = fit1d$ddist(res[is.IC,2]) ## Density
    S1.cur = exp(-H1.cur) ; S2.cur = exp(-H2.cur) ## Survival
    ##
    idx1 = 1:q1 ; idx2 = parms1$q + (1:q2)
    if (J1 > 0) P1.cur = Pcal.fun(nfixed1,lambda1.cur,Pd1.x) ## Update penalty matrix for location
    else P1.cur = diag(0,length(psi1))
    ##
    if (J2 > 0) P2.cur = Pcal.fun(nfixed2,lambda2.cur,Pd2.x) ## Update penalty matrix for dispersion
    else P2.cur = diag(0,length(psi2))
    ## Global penalty matrix
    P.cur = diag(0,q1+q2)
    if (J1 > 0) P.cur[idx1,idx1] = P1.cur
    if (J2 > 0) P.cur[idx2,idx2] = P2.cur
    ##
    ## lpost
    ## -----
    Dif.S = pmax(1e-6,S1.cur-S2.cur)
    lpost.cur = llik.cur = sum(w * ifelse(is.IC,
                                          log(Dif.S),
                                          event*(log(h1.cur)-log(sd.cur)) - H1.cur))
    eps.ev = 1e-10
    if (J1 > 0){
        ## Prior on <lambda1>
        lpost.cur = lpost.cur + sum(dgamma(lambda1.cur,a.pen,b.pen,log=TRUE)) ## Prior on <lambda1>
        ## Prior on (psi1 | lambda1)
        ev.psi1 = svd(P1.cur,nu=0,nv=0)$d
        lpost.cur = lpost.cur + .5*sum(log(ev.psi1[ev.psi1>eps.ev])) - .5*sum(psi1.cur*c(P1.cur%*%psi1.cur))
    }
    if (J2 > 0){
        ## Prior on <lambda2>
        lpost.cur = lpost.cur + sum(dgamma(lambda2.cur,a.pen,b.pen,log=TRUE)) ## Prior on <lambda2>
        ## Prior on (psi2 | lambda2)
        ev.psi2 = svd(P2.cur,nu=0,nv=0)$d
        lpost.cur = lpost.cur + .5*sum(log(ev.psi2[ev.psi2>eps.ev])) - .5*sum(psi2.cur*c(P2.cur%*%psi2.cur))
    }
    ## Score
    ## -----
    U.psi = U.psi1 = U.psi2 = NULL
    if (gradient){
        ## <psi1>
        temp.Upsi1 = ifelse(is.IC,
                            -(1/sd.cur)*(f1.cur-f2.cur)/(Dif.S),
                            c(((event/sd.cur)*dB1)%*%phi.cur -1/sd.cur*h1.cur))
        omega.psi1 = -temp.Upsi1
        ## U.psi1 should be multiplied by <k> if (w -> k*w)  (hence w * omega.psi1 !!) :
        U.psi1 = c(t(Xcal.1) %*% (w * omega.psi1)) ## Score.psi1
        if (J1 > 0) U.psi1 = U.psi1 - c(P1.cur%*%psi1.cur) ## Score.psi1 if additive component(s)
        ##
        ## <psi2>
        temp.Upsi2 = ifelse(is.IC,
                            -(res[,1]*f1.cur-res[,2]*f2.cur)/(Dif.S),
                            c((event*res[,1]*dB1)%*%phi.cur + event-res[,1]*h1.cur))
        omega.psi2 = -temp.Upsi2
        ## U.psi2 should be multiplied by <k> if (w -> k*w)  (hence w * omega.psi2 !!) :
        U.psi2 = c(t(Xcal.2) %*% (w * omega.psi2)) ## Score.psi2
        if (J2 > 0) U.psi2 = U.psi2 - c(P2.cur%*%psi2.cur) ## Score.psi2 if additive component(s)
        ## ALTERNATIVE to get var-cov from score
        Ui.psi1 = Xcal.1 * omega.psi1 ## Xcal.1 * (w * omega.psi1)
        Ui.psi2 = Xcal.2 * omega.psi2 ## Xcal.2 * (w * omega.psi2)
        Ui.psi = cbind(Ui.psi1, Ui.psi2)
        U0.psi = colSums(w * Ui.psi)
        Ub.psi = U0.psi - c(P.cur %*% psi)
        Vn.psi = crossprod(sqrt(w) * Ui.psi)- (1/sw) * outer(U0.psi,U0.psi)
        idx1 = 1:q1 ; idx2 = (q1+1):(q1+q2)
        if (J1 > 0) Vn.psi[idx1,idx1] = Vn.psi[idx1,idx1] + P1.cur
        if (J2 > 0) Vn.psi[idx2,idx2] = Vn.psi[idx2,idx2] + P2.cur
        ##
        ## temp = cbind(Xcal.1 * (sqrt(w) * omega.psi1), Xcal.2 * (sqrt(w) * omega.psi2))
        ## Vn.psi = crossprod(temp) ## Vn.psi should be multiplied by <k> if (w -> k*w)  (hence sqrt(w) !!)
        ## idx1 = 1:q1 ; idx2 = (q1+1):(q1+q2)
        ## if (J1 > 0) Vn.psi[idx1,idx1] = Vn.psi[idx1,idx1] + P1
        ## if (J2 > 0) Vn.psi[idx2,idx2] = Vn.psi[idx2,idx2] + P2
        ## End alternative
        ##
        U.psi = c(U.psi1,U.psi2) ## Score for <psi> = <psi1,psi2>
        if (any(is.na(U.psi))) return(stopthis())
    }
    ##
    ## Computation of Hessian
    ## ----------------------
    Mcal.1 = Mcal.2 = NULL ## Quantities required to update <lambda1> & <lambda2>
    if (hessian){
      ## <psi1>
      g1.cur = c(dB1%*%phi.cur) - h1.cur ; g2.cur = c(dB2%*%phi.cur) - h2.cur
      W.psi1 = w * (1/sd.cur^2)*
        ifelse(is.IC,
               (f1.cur*g1.cur-f2.cur*g2.cur)/(Dif.S) + (f1.cur-f2.cur)^2/(Dif.S)^2,
               -c((event*d2B-h1.cur*dB1)%*%phi.cur)) ## Weight.psi1
      Hes.psi1 = Hes0.psi1 = -t(Xcal.1)%*%(W.psi1*Xcal.1) ## Hessian.psi1
      if (J1 > 0) {
        Hes.psi1 = Hes.psi1 - P1.cur ## Hessian.psi1 if additive component(s)
        tZWB.1 = t(Z1)%*%(W.psi1*Bcal.1) ; tZWZ.1 = t(Z1)%*%(W.psi1*Z1)
        inv.tZWZ.1 = try(solve(tZWZ.1),TRUE)
        if (!is.matrix(inv.tZWZ.1)) return(stopthis()) ## stopthis() ## PHL
        Mcal.1 = t(Bcal.1)%*%(W.psi1*Bcal.1) -t(tZWB.1)%*%inv.tZWZ.1%*%tZWB.1 ## Required to update <lambda1>
      }
      ## <psi2>
      m1.cur = 1 + res[,1]*(c(dB1%*%phi.cur)-h1.cur) ; m2.cur = 1 + res[,2]*(c(dB2%*%phi.cur)-h2.cur)
      W.psi2 = w * ifelse(is.IC,
                      (res[,1]*f1.cur*m1.cur-res[,2]*f2.cur*m2.cur)/(Dif.S) +
                        (res[,1]*f1.cur-res[,2]*f2.cur)^2/(Dif.S)^2,
                      c(-event*(res[,1]^2)*c(d2B%*%phi.cur)
                        -res[,1]*(event-h1.cur*res[,1])*c(dB1%*%phi.cur)
                        +h1.cur*res[,1])) ## Weight.psi2
      Hes.psi2 = Hes0.psi2 = -t(Xcal.2)%*%(W.psi2*Xcal.2) ## Hessian.psi2
      if (J2 > 0){
        Hes.psi2 = Hes.psi2 - P2.cur ## Hessian.psi2 if additive component(s)
        tZWB.2 = t(Z2)%*%(W.psi2*Bcal.2) ; tZWZ.2 = t(Z2)%*%(W.psi2*Z2)
        inv.tZWZ.2 = try(solve(tZWZ.2),TRUE)
        if (!is.matrix(inv.tZWZ.2)) return(stopthis()) ## stopthis() ## PHL
        Mcal.2 = t(Bcal.2)%*%(W.psi2*Bcal.2) -t(tZWB.2)%*%inv.tZWZ.2%*%tZWB.2 ## Required to update <lambda2>
      }
      ## Cross-derivatives
      W.cross = w * (1/sd.cur)*
        ifelse(is.IC,
               (f1.cur-f2.cur)/(Dif.S) +
                 (res[,1]*f1.cur*g1.cur-res[,2]*f2.cur*g2.cur)/(Dif.S) +
                 (f1.cur-f2.cur)*(res[,1]*f1.cur-res[,2]*f2.cur)/(Dif.S)^2,
               -event*res[,1]*c(d2B%*%phi.cur) - (event-h1.cur*res[,1])*c(dB1%*%phi.cur) + h1.cur)
        ## Cross-Weight.psi1.psi2
      Hes.21 = -t(Xcal.2)%*%(W.cross*Xcal.1) ## Hessian for cross derivatives
      Hes0.psi = cbind(rbind(Hes0.psi1, Hes.21),
                       rbind(t(Hes.21), Hes0.psi2))
      Hes.psi = cbind(rbind(Hes.psi1, Hes.21),
                      rbind(t(Hes.21), Hes.psi2))
      if(any(is.nan(Hes.psi))) return(stopthis()) ## stopthis() ## PHL
      Cov.psi = try(solve(-Hes.psi),TRUE) ## First computation of Cov.psi before REML
      if (!is.matrix(Cov.psi)){
        Hes.psi = cbind(rbind(Hes.psi1, 0*Hes.21),
                        rbind(0*t(Hes.21), Hes.psi2))
        Cov.psi = try(solve(-Hes.psi),TRUE)
      }
      if(!is.matrix(Cov.psi)) return(stopthis()) ## stopthis() ## PHL
    }
    ## REML-type correction to Hes.psi2
    ## --------------------------------
    dHes.REML = ddiagHes.REML = NULL
    if (REML){
      if(!is.matrix(Cov.psi)) return(stopthis()) ## PHL
      Cov.psi1 = Cov.psi[1:ncol(Xcal.1),1:ncol(Xcal.1)] ## Will be recomputed further using global Hessian...
      ##
      Mat.k = vector("list",q2) ## Prepare list with <q2> NULL elements
      dU.REML = rep(0,q2)  ## Change in Gradient U.psi2 due to REML
      dHes.REML = 0*diag(q2) ## Change in Hessian Hes.psi2 due to REML
      ddiagHes.REML = rep(0,q2) ## Change in diag(Hessian) diag(Hes.psi2) due to REML
      ##
      for (kk in 1:q2){
        if (is.null(Mat.k[[kk]])) Mat.k[[kk]] = Cov.psi1%*% (t(Xcal.1) %*% (-2*Xcal.2[,kk]*W.psi1*Xcal.1))
        Mat.k2 = Cov.psi1%*% (t(Xcal.1) %*% (4*Xcal.2[,kk]*Xcal.2[,kk]*W.psi1*Xcal.1))
        dU.REML[kk] = -.5*sum(diag(Mat.k[[kk]])) ##-.5*sum(diag(Mat.k))
        if (diag.only) ddiagHes.REML = -.5*sum(diag(Mat.k2)) +.5*sum(t(Mat.k[[kk]])*Mat.k[[kk]])
        else {
          for (ll in kk:q2){
            if (is.null(Mat.k[[ll]])) Mat.k[[ll]] = Cov.psi1%*% (t(Xcal.1) %*% (-2*Xcal.2[,ll]*W.psi1*Xcal.1))
            Mat.kl = Cov.psi1%*% (t(Xcal.1) %*% (4*Xcal.2[,kk]*Xcal.2[,ll]*W.psi1*Xcal.1))
            dHes.REML[kk,ll] = dHes.REML[ll,kk] = -.5*sum(diag(Mat.kl)) +.5*sum(t(Mat.k[[kk]])*Mat.k[[ll]])
          }
        }
      }
      U.psi2 = U.psi2 + dU.REML ##[1:nfixed2]
      if (diag.only) diag(Hes.psi2) = diag(Hes.psi2) + ddiagHes.REML ##[1:nfixed2]
      else Hes.psi2 = Hes.psi2 + dHes.REML
      ##
      ## Recomputation of the Score after REML
      U.psi = c(U.psi1,U.psi2)
      ## Recomputation of Hes.psi for <psi> = <psi1,psi2> after REML
      Hes.psi = cbind(rbind(Hes.psi1, Hes.21),
                      rbind(t(Hes.21), Hes.psi2))
      ## Hes.psi = -adjustPD(-Hes.psi)$PD
      Cov.psi = try(solve(-Hes.psi),TRUE) ## First computation of Cov.psi before REML
      if (!is.matrix(Cov.psi)){
        Hes.psi = cbind(rbind(Hes.psi1, 0*Hes.21),
                        rbind(0*t(Hes.21), Hes.psi2))
        Cov.psi = try(solve(-Hes.psi),TRUE)
      }
    }
    if (any(is.na(U.psi))) return(stopthis()) ## stopthis() ## PHL
    ## REML-type correction to log-posterior
    ldet.cur = NULL
    if (REML){
      if(any(is.nan(Hes.psi1))) return(stopthis()) ## stopthis() ## PHL
      ldet.cur = logDet.fun(-Hes.psi1) ## logdet(X1'W1 X1 + K1)
      lpost.cur = lpost.cur -.5*ldet.cur  ## REML: add -.5*logdet(X1'W1 X1 + K1)
    }
    ##
    ## Effective dimensions for additive terms
    ## ---------------------------------------
    if (hessian){
        EDF = EDF.fun(Hes0.psi, lambda1.cur, lambda2.cur, parms1, parms2, joint.computation=joint.computation)
        ED1 = EDF$ED1 ; ED2 = EDF$ED2 ; ED.tot = EDF$ED.tot
    } else {
        ED1 = NULL ; ED2 = NULL ; ED.tot = NULL
    }
    ##
    ## Output
    ## ------
    if (any(is.na(U.psi))) return(stopthis())
    ans = list(y=y, w=w, data=data,
               formula1=formula1, formula2=formula2,
               psi1=psi1,psi2=psi2, lambda1=lambda1,lambda2=lambda2,
               psi=c(psi1,psi2),
               res=res, mu=mu.cur, sd=sd.cur,
               lpost=lpost.cur, llik=llik.cur, dev=-2*llik.cur,
               h1=h1.cur, h2=h2.cur, H1=H1.cur, H2=H2.cur,
               f1=f1.cur, f2=f2.cur, S1=S1.cur, S2=S2.cur,
               P1=P1.cur, P2=P2.cur,
               ED1=ED1, ED2=ED2, ED.tot=ED.tot,
               Mcal.1=Mcal.1, Mcal.2=Mcal.2,
               ldet=ldet.cur,
               REML=REML,hessian=hessian)
    if (gradient){
        ans = c(ans,
               list(U.psi1=U.psi1, U.psi2=U.psi2, U.psi=U.psi,
                    omega.psi1=omega.psi1, omega.psi2=omega.psi2,
                    Vn.psi=Vn.psi)
               )
    }
    if (hessian){
      ans = c(ans,
              list(Hes.psi1=Hes.psi1, Hes.psi2=Hes.psi2,
                   Hes0.psi1=Hes0.psi1, Hes0.psi2=Hes0.psi2,
                   Hes.psi=Hes.psi, Hes0.psi=Hes0.psi,
                   Cov.psi=Cov.psi,
                   dHes.REML=dHes.REML, ddiagHes.REML=ddiagHes.REML)
      )
    }
    ##
    return(ans)
  } ## End <ff> function
  ##
  ##
  ## -------------------------------------------------
  ## Tools for the selection of the penalty parameters
  ## -------------------------------------------------
  ##
  ## Log-determinant of a positive definite matrix based on determinant(.)
  ##  Note: faster than methods based on functions cholevski, svd or qr
  ## ---------------------------------------------------------------------
  ldet.fun = function(x) c(determinant(x,logarithm=TRUE)$mod)
  ## ldet.fun = function(x) 2*sum(log(diag(chol(x))))
  ## ldet.fun = function(x) sum(log(svd(A,nu=0,nv=0)$d))
  ## ldet.fun  = function(x) sum(log(abs(diag(qr(x)$qr))))
  ##
  ## Fast calculation of log |B'WB + lambda*Pd| using pre-computed eigenvalues:
  ##  log |B'WB + lambda*Pd| = ldet0 + sum(log(lambda1)) + sum_j(log(tau+dj))
  ## --------------------------------------------------------------------------
  ## Starting from B'WB and Pd
  ev.fun <- function(BwB, Pd){  ## t(B) %*% diag(w) %*% B
      K = ncol(Pd) ; rk = qr(Pd)$rank ; id1 = 1:rk ; id2 = (1:K)[-id1]
      temp2 = svd(Pd) ; U = temp2$u ; lambda1 = temp2$d[id1]
      BwB = t(U) %*% (BwB %*% U)
      M = BwB[id1,id1] - BwB[id1,id2,drop=FALSE]%*%solve(BwB[id2,id2,drop=FALSE])%*%BwB[id2,id1,drop=FALSE]
      MM = sqrt(1/lambda1) * t(sqrt(1/lambda1)*t(M))
      ## MM2 = diag(1/sqrt(lambda1))%*%M%*%diag(1/sqrt(lambda1))
      dj = svd(MM,nu=0,nv=0)$d
      ldet0 = ldet.fun(BwB[id2,id2,drop=FALSE])
      ## log |B'WB + lambda*Pd| = ldet0 + sum(log(lambda1)) + sum_j(log(tau+dj))
      return(list(ldet0=ldet0,lambda1=lambda1,dj=dj))
  }
  ##
  ## Function 1: Penalty parameter selection using Laplace approximation to p(lambda|data)
  ##    with Newton-Raphson method
  ## --------------------------------------------------------------------------------------
  select.lambda.LPS <- function(psi, lambda1, lambda2, itermax = 50) {
      psi1.cur = psi[1:q1] ; psi2.cur = psi[-(1:q1)]
      lambda1.mat <- lambda2.mat <- NULL
      ## Euclidean norm function
      L2norm <- function(x) sqrt(sum(x^2))
      ##
      ## Function to update penalty matrix and calculate score and Hessian
      update.penalty <- function(psi.cur, lambda, nfixed, Pd.x, J, K, pen.order, Mcal, lambda.min, add.lab, verbose) {
          iter.lam = 0
          lambda.mat = NULL
          ## Home-made Newton-Raphson
          ok.xi = ifelse(J > 0,FALSE,TRUE) ; iter.lam = 0
          while(!ok.xi){
              iter.lam = iter.lam + 1
              P.cur = Pcal.fun(nfixed,lambda,Pd.x) ## Update penalty matrix
              iMpen = solve(Mcal + P.cur[-(1:nfixed),-(1:nfixed)])
              ## Score.lambda
              U.lam = U.xi = rep(0,J)
              Rj = list()
              ##
              quad.cur = rep(0,J)
              for (j in 1:J){
                  lam.j = rep(0,J) ; lam.j[j] = lambda[j]
                  Pj = Pcal.fun(nfixed,lam.j,Pd.x)[-(1:nfixed),-(1:nfixed)] ## Update penalty matrix
                  Rj[[j]] = iMpen%*%(Pj/lam.j[j])
                  idx = nfixed + (j-1)*K + (1:K)
                  theta.j = psi.cur[idx]
                  quad.cur[j] = sum(theta.j*c(Pd.x%*%theta.j)) + b.pen  ## <-----
                  U.lam[j] = .5*(K-pen.order)/lam.j[j] -.5*quad.cur[j] -.5*sum(diag(Rj[[j]]))
              }
              ## Score.xi  where  lambda = lambda.min + exp(xi)
              U.xi = U.lam * (lambda - lambda.min)
              ##
              ## Hessian.lambda
              Hes.lam = diag(0,J)
              Only.diagHes = FALSE ## PHL
              if (!Only.diagHes){
                  for (j in 1:J){
                      for (k in j:J){
                          temp = 0
                          if (j==k) temp = temp -.5*(K-pen.order)/lambda[j]^2
                          temp = temp + .5*sum(t(Rj[[j]])*Rj[[k]]) ## tr(XY) = sum(t(X)*Y)
                          Hes.lam[j,k] = Hes.lam[k,j] = temp
                      }
                  }
              } else {
                  for (j in 1:J){
                      Hes.lam[j,j] = -.5*(K-pen.order)/lambda[j]^2 +.5*sum(t(Rj[[j]])*Rj[[j]])
                  }
              }
              ## Hessian.xi  where  lambda = lambda.min + exp(xi)
              Hes.xi = t((lambda-lambda.min) * t((lambda-lambda.min) * Hes.lam))
              diag(Hes.xi) = diag(Hes.xi) + U.lam * (lambda-lambda.min)
              ## Update xi  where  lambda = lambda.min + exp(xi)
              xi.cur = log(lambda-lambda.min) ; names(xi.cur) = add.lab
              dxi = -c(MASS::ginv(Hes.xi) %*% U.xi) ## Increment.xi
              ##
              step = 1
              accept = FALSE
              while(!accept){
                  xi.prop = xi.cur + step*dxi
                  accept = all(xi.prop < 10) ## all(is.finite(exp(xi.prop)))
                  if (!accept) step = .5*step
              }
              if ((step !=1)&verbose) cat("Step-halving for lambda: step=",step,"\n")
              xi.cur = xi.prop
              lambda = lambda.min + exp(xi.cur)
              ## Convergence ?
              ok.xi = (L2norm(U.xi) < grad.tol) | (iter.lam > itermax)
          }
          xi.se = sqrt(pmax(1e-6,diag(MASS::ginv(-Hes.xi)))) ## sqrt(diag(solve(-Hes.xi)))
          xi.low = xi.cur - z.alpha*xi.se ; xi.up = xi.cur + z.alpha*xi.se
          lambda.low = lambda.min + exp(xi.low) ; lambda.up = lambda.min + exp(xi.up)
          lambda.mat = cbind(lambda=lambda, low=lambda.low, up=lambda.up)
          rownames(lambda.mat) = add.lab
          ##
          return(cbind(lambda = lambda, low = lambda.low, up = lambda.up))
      } ## End <update.penalty>
      ##
      ## -- Location: <lambda1> update --
      if (J1 > 0) {
          lambda1.mat <- update.penalty(psi1.cur, lambda1,
                                        nfixed1, Pd1.x, J1, K1, pen.order.1,
                                        obj.cur$Mcal.1, lambda1.min, addloc.lab,
                                        verbose)
          rownames(lambda1.mat) <- addloc.lab
      }
      ## -- Dispersion: <lambda2> update --
      if (J2 > 0) {
          lambda2.mat <- update.penalty(psi2.cur, lambda2,
                                        nfixed2, Pd2.x, J2, K2, pen.order.2,
                                        obj.cur$Mcal.2, lambda2.min, adddisp.lab,
                                        verbose)
          rownames(lambda2.mat) <- adddisp.lab
      }
      ##
      return(list(lambda1 = lambda1.mat[,1], lambda2 = lambda2.mat[,1],
                  lambda1.mat = lambda1.mat, lambda2.mat = lambda2.mat))
  } ## End select.lambda.LPS
  ##
  ## Function 2: Penalty parameter selection using Laplace approximation to  p(lambda|data)
  ##    with an underlying computation of |B'WB + lamdba Pd| using precomputed eigenvalues
  ##    and the fixed point method to find lambda.MAP
  ## --------------------------------------------------------------------------------------
  select.lambda.LPS2 <- function(psi, lambda1,lambda2, itermax=50){
      ## logscale = TRUE ## Maximize p(log(lambda)|D) when TRUE,  p(lambda|D) otherwise
      ## Generic update of a single penalty parameter
      ##
      psi1.cur = psi[1:q1] ; psi2.cur[-(1:q1)]
      ##
      ## -- Update of a single penalty parameter --
      update.lambda.fun <- function(lambda,quad,ev,rk){
          ok.lam = FALSE ; iter.lam = 0
          while(!ok.lam){
              iter.lam = iter.lam + 1
              lambda.old = lambda
              ttr = sum(1 / (lambda+ev))
              if (!logscale){
                  ## Posterior mode of lambda
                  lambda = ((a.pen-1) + .5*rk) / (b.pen + .5*(quad + ttr))
              } else {
                  ## Exponential of the posterior mode of log(lambda)
                  lambda = (a.pen + .5*rk) / (b.pen + .5*(quad + ttr))
              }
              lam.dif = abs((lambda - lambda.old)/lambda)
              ok.lam = (lam.dif < .001) | (iter.lam >= itermax)
          }
          curv = (b.pen+.5*quad)/lambda + .5*sum(ev/lambda/(lambda+ev)^2)
          se.lambda = 1 / sqrt(curv)
          ## se.lambda = 1 / sqrt(.5 * sum(ev/(lambda*(lambda+ev)^2)))
          se.loglambda = 1 / sqrt(lambda^2 * curv) ## <-- HERE
          ##
          return(list(lambda=lambda, se.lambda=se.lambda, se.loglambda=se.loglambda))
      }
      ##
      lambda1.mat = se.lambda1 = se.loglambda1 = NULL
      lambda2.mat = se.lambda2 = se.loglambda2 = NULL
      ## -- Update <lambda1> if (J1 > 0) --
      if (J1 > 0){
          rk1 = qr(regr1$Pd.x)$rank
          Hes0.psi1 = obj.cur$Hes0.psi[1:q1,1:q1]
          se.lambda1 = se.loglambda1 = 0*lambda1
          ev1.lst = list()
          for (j in 1:J1){ ## Loop over additive terms in the location sub-model
              ## Update lambda1
              idx = nfixed1 + (j-1)*K1 + (1:K1)
              theta.j = psi1.cur[idx]
              quad.j = sum(theta.j*c(Pd1.x%*%theta.j))
              Hes0.psi1j = Hes0.psi1[idx,idx,drop=FALSE]
              ev1.lst[[j]] = ev.fun(BwB=Hes0.psi1j,Pd=regr1$Pd.x)$dj
              temp = update.lambda.fun(lambda1[j],quad.j,ev1.lst[[j]],rk1)
              lambda1[j] = temp$lambda
              se.lambda1[j] = temp$se.lambda
              se.loglambda1[j] = temp$se.loglambda
          }
          lambda.low = exp(log(lambda1) - z.alpha*se.loglambda1)
          lambda.up  = exp(log(lambda1) + z.alpha*se.loglambda1)
          lambda1.mat = cbind(lambda=lambda1, low=lambda.low, up=lambda.up)
          rownames(lambda1.mat) = addloc.lab
      }
      ## -- Update <lambda2> if (J2 > 0) --
      if (J2 > 0){
          rk2 = qr(regr2$Pd.x)$rank
          Hes0.psi2 = obj.cur$Hes0.psi[q1+(1:q2), q1+(1:q2)]
          se.lambda2 = se.loglambda2 = 0*lambda2
          ev2.lst = list()
          for (j in 1:J2){ ## Loop over additive terms in the log-dispersion sub-model
              ## Update lambda2.cur
              idx = nfixed2 + (j-1)*K2 + (1:K2)
              theta.j = psi2.cur[idx] ## psi.cur[idx] ## <---- Correction 2025.02.21 !!
              quad.j = sum(theta.j*c(Pd2.x%*%theta.j))
              Hes0.psi2j = Hes0.psi2[idx,idx,drop=FALSE]
              ev2.lst[[j]] = ev.fun(BwB=Hes0.psi2j,Pd=regr2$Pd.x)$dj
              temp = update.lambda.fun(lambda2[j],quad.j,ev2.lst[[j]],rk2)
              lambda2[j] = temp$lambda
              se.lambda2[j] = temp$se.lambda
              se.loglambda2[j] = temp$se.loglambda
          }
          lambda.low = exp(log(lambda2) - z.alpha*se.loglambda2)
          lambda.up  = exp(log(lambda2) + z.alpha*se.loglambda2)
          lambda2.mat = cbind(lambda=lambda2, low=lambda.low, up=lambda.up)
          rownames(lambda2.mat) = adddisp.lab
      }
      ##
      return(list(lambda1=lambda1, lambda2=lambda2,
                  lambda1.mat=lambda1.mat, lambda2.mat=lambda2.mat,
                  se.lambda1=se.lambda1, se.lambda2=se.lambda2,
                  se.loglambda1=se.loglambda1, se.loglambda2=se.loglambda2))
  } ## End select.lambda.LPS2
  ##
  ##
  ## ## Function 1b: Penalty parameter selection using Laplace approximation to p(lambda|data)
  ## ##    with Newton-Raphson method
  ## ##    Note: same as Function 1, but with less structured code.
  ## ## --------------------------------------------------------------------------------------
  ## select.lambda.LPS1b <- function(psi, lambda1,lambda2, itermax=50){
  ##     psi1.cur = psi[1:q1] ;  psi2.cur = psi[-(1:q1)]
  ##     iter.lam1 = iter.lam2 = 0
  ##     lambda1.mat = lambda2.mat = NULL
  ##     ## Euclidean norm function
  ##     L2norm <- function(x) sqrt(sum(x^2))
  ##     ## -3a- LOCATION
  ##     ## -------------
  ##     if (J1 > 0){
  ##         ## Home-made Newton-Raphson
  ##         ok.xi1 = ifelse(J1 > 0,FALSE,TRUE) ; iter.lam1 = 0
  ##         Mcal.1 = obj.cur$Mcal.1
  ##         while(!ok.xi1){
  ##             iter.lam1 = iter.lam1 + 1
  ##             P1.cur = Pcal.fun(nfixed1,lambda1,Pd1.x) ## Update penalty matrix
  ##             iMpen.1 = solve(Mcal.1 + P1.cur[-(1:nfixed1),-(1:nfixed1)])
  ##             ## Score.lambda1
  ##             U.lam1 = U.xi1 = rep(0,J1)
  ##             Rj.1 = list()
  ##             ##
  ##             quad1.cur = rep(0,J1)
  ##             for (j in 1:J1){
  ##                 lam1.j = rep(0,J1) ; lam1.j[j] = lambda1[j]
  ##                 Pj.1 = Pcal.fun(nfixed1,lam1.j,Pd1.x)[-(1:nfixed1),-(1:nfixed1)] ## Update penalty matrix
  ##                 Rj.1[[j]] = iMpen.1%*%(Pj.1/lam1.j[j])
  ##                 idx = nfixed1 + (j-1)*K1 + (1:K1)
  ##                 theta.j = psi1.cur[idx]
  ##                 quad1.cur[j] = sum(theta.j*c(Pd1.x%*%theta.j)) + b.pen  ## <-----
  ##                 U.lam1[j] = .5*(K1-pen.order.1)/lam1.j[j] -.5*quad1.cur[j] -.5*sum(diag(Rj.1[[j]]))
  ##             }
  ##             ## Score.xi1  where  lambda1 = lambda1.min + exp(xi1)
  ##             U.xi1 = U.lam1 * (lambda1 - lambda1.min)
  ##             ##
  ##             ## Hessian.lambda1
  ##             Hes.lam1 = diag(0,J1)
  ##             Only.diagHes1 = FALSE ## PHL
  ##             if (!Only.diagHes1){
  ##                 for (j in 1:J1){
  ##                     for (k in j:J1){
  ##                         temp = 0
  ##                         if (j==k) temp = temp -.5*(K1-pen.order.1)/lambda1[j]^2
  ##                         temp = temp + .5*sum(t(Rj.1[[j]])*Rj.1[[k]]) ## tr(XY) = sum(t(X)*Y)
  ##                         Hes.lam1[j,k] = Hes.lam1[k,j] = temp
  ##                     }
  ##                 }
  ##             } else {
  ##                 for (j in 1:J1){
  ##                     Hes.lam1[j,j] = -.5*(K1-pen.order.1)/lambda1[j]^2 +.5*sum(t(Rj.1[[j]])*Rj.1[[j]])
  ##                 }
  ##             }
  ##             ## Hessian.xi1  where  lambda1 = lambda1.min + exp(xi1)
  ##             Hes.xi1 = t((lambda1-lambda1.min) * t((lambda1-lambda1.min) * Hes.lam1))
  ##             diag(Hes.xi1) = diag(Hes.xi1) + U.lam1 * (lambda1-lambda1.min)
  ##             ## Update xi1  where  lambda1 = lambda1.min + exp(xi1)
  ##             xi1.cur = log(lambda1-lambda1.min) ; names(xi1.cur) = addloc.lab
  ##             dxi1 = -c(MASS::ginv(Hes.xi1) %*% U.xi1) ## Increment.xi1
  ##             ##
  ##             step = 1
  ##             accept = FALSE
  ##             while(!accept){
  ##                 xi1.prop = xi1.cur + step*dxi1
  ##                 accept = all(xi1.prop < 10) ## all(is.finite(exp(xi1.prop)))
  ##                 if (!accept) step = .5*step
  ##             }
  ##             if ((step !=1)&verbose) cat("Step-halving for lambda1: step=",step,"\n")
  ##             xi1.cur = xi1.prop
  ##             lambda1 = lambda1.min + exp(xi1.cur)
  ##             ##
  ##             mat = Mcal.1 + P1.cur[-(1:nfixed1),-(1:nfixed1)]
  ##             ## Convergence ?
  ##             ok.xi1 = (L2norm(U.xi1) < grad.tol) | (iter.lam1 > itermax)
  ##         }
  ##         xi1.se = sqrt(pmax(1e-6,diag(MASS::ginv(-Hes.xi1)))) ## sqrt(diag(solve(-Hes.xi1)))
  ##         xi1.low = xi1.cur - z.alpha*xi1.se ; xi1.up = xi1.cur + z.alpha*xi1.se
  ##         lambda1.low = lambda1.min + exp(xi1.low) ; lambda1.up = lambda1.min + exp(xi1.up)
  ##         lambda1.mat = cbind(lambda=lambda1, low=lambda1.low, up=lambda1.up)
  ##         rownames(lambda1.mat) = addloc.lab
  ##     }
  ##     ##
  ##     ## -3b- DISPERSION
  ##     ## ---------------
  ##     ##
  ##     if (J2 > 0){
  ##         ## Home-made Newton-Raphson
  ##         ok.xi2 = ifelse(J2 > 0,FALSE,TRUE) ; iter.lam2 = 0
  ##         Mcal.2 = obj.cur$Mcal.2
  ##         while(!ok.xi2){
  ##             iter.lam2 = iter.lam2 + 1
  ##             P2.cur = Pcal.fun(nfixed2,lambda2,Pd2.x) ## Update penalty matrix
  ##             iMpen.2 = MASS::ginv(Mcal.2 + P2.cur[-(1:nfixed2),-(1:nfixed2)])
  ##             ## Score.lambda2
  ##             U.lam2 = U.xi2 = rep(0,J2)
  ##             Rj.2 = list()
  ##             quad2.cur = rep(0,J2)
  ##             for (j in 1:J2){
  ##                 lam2.j = rep(0,J2) ; lam2.j[j] = lambda2[j]
  ##                 Pj.2 = Pcal.fun(nfixed2,lam2.j,Pd2.x)[-(1:nfixed2),-(1:nfixed2)] ## Update penalty matrix
  ##                 Rj.2[[j]] = iMpen.2%*%(Pj.2/lam2.j[j])
  ##                 idx = nfixed2 + (j-1)*K2 + (1:K2)
  ##                 theta2.j = psi2.cur[idx]
  ##                 quad2.cur[j] = sum(theta2.j*c(Pd2.x%*%theta2.j)) + b.pen  ## <-----
  ##                 U.lam2[j] = .5*(K2-pen.order.2)/lam2.j[j] -.5*quad2.cur[j] -.5*sum(diag(Rj.2[[j]]))
  ##             }
  ##             ## Score.xi2  where  lambda2 = lambda2.min + exp(xi2)
  ##             U.xi2 = U.lam2 * (lambda2-lambda2.min)
  ##             ##
  ##             ## Hessian.lambda2
  ##             Hes.lam2 = diag(0,J2)
  ##             Only.diagHes2 = FALSE
  ##             if (!Only.diagHes2){
  ##                 for (j in 1:J2){
  ##                     for (k in j:J2){
  ##                         temp = 0
  ##                         if (j==k) temp = temp -.5*(K2-pen.order.2)/lambda2[j]^2
  ##                         temp = temp +.5*sum(t(Rj.2[[j]])*Rj.2[[k]]) ## tr(XY) = sum(t(X)*Y)
  ##                         Hes.lam2[j,k] = Hes.lam2[k,j] = temp
  ##                     }
  ##                 }
  ##             } else {
  ##                 for (j in 1:J2){
  ##                     Hes.lam2[j,j] = -.5*(K2-pen.order.2)/lambda2[j]^2 +.5*sum(t(Rj.2[[j]])*Rj.2[[j]])
  ##                 }
  ##             }
  ##             ## Hessian.xi2  where  lambda2 = lambda2.min + exp(xi2)
  ##             Hes.xi2 = t((lambda2-lambda2.min) * t((lambda2-lambda2.min) * Hes.lam2))
  ##             diag(Hes.xi2) = diag(Hes.xi2) + U.lam2 * (lambda2-lambda2.min)
  ##             ## Update xi2  where  lambda2 = lambda2.min + exp(xi2)
  ##             xi2.cur = log(lambda2-lambda2.min) ; names(xi2.cur) = adddisp.lab
  ##             dxi2 = -c(MASS::ginv(Hes.xi2) %*% U.xi2) ## Increment.xi2
  ##             ##
  ##             step = 1
  ##             accept = FALSE
  ##             while(!accept){
  ##                 xi2.prop = xi2.cur + step*dxi2
  ##                 accept = all(xi2.prop < 10)
  ##                 if (!accept) step = .5*step
  ##             }
  ##             if ((step !=1)&verbose) cat("Step-halving for lambda2: step=",step,"\n")
  ##             xi2.cur = xi2.prop
  ##             lambda2 = lambda2.min + exp(xi2.cur)
  ##             ## Convergence ?
  ##             ok.xi2 = (L2norm(U.xi2) < grad.tol) | (iter.lam2 > itermax)
  ##         }
  ##         xi2.se = sqrt(pmax(1e-6,diag(MASS::ginv(-Hes.xi2)))) ## sqrt(diag(solve(-Hes.xi2)))
  ##         xi2.low = xi2.cur - z.alpha*xi2.se ; xi2.up = xi2.cur + z.alpha*xi2.se
  ##         lambda2.low = lambda2.min + exp(xi2.low) ; lambda2.up = lambda2.min + exp(xi2.up)
  ##         lambda2.mat = cbind(lambda=lambda2, low=lambda2.low, up=lambda2.up)
  ##         rownames(lambda2.mat) = adddisp.lab
  ##     }
  ##     ##
  ##     ## print(lambda1.mat)
  ##     ## print(U.xi1)
  ##     ## cat("\n")
  ##     return(list(lambda1=lambda1, lambda2=lambda2,
  ##                 lambda1.mat=lambda1.mat, lambda2.mat=lambda2.mat))
  ## } ## End select.lambda.LPS1b
  ##
  ##
  ## ----------------------------------
  ##  Generic Newton-Raphson algorithm
  ## ---------------------------------
  ## GOAL: maximize function using its gradient and Hessian with Newton-Raphson
  ## INPUT:
  ##    * g: function returning a list:
  ##       - g: function value to be maximized
  ##       - grad: gradient of the function
  ##       - RDM: grad' (-Hes)^-1 grad / ntheta
  ##       - dtheta: solve(-Hes) %*% grad
  ##       - RDM: sum(grad * dtheta) / ntheta
  ##       - dtheta2: solve(Vn, grad) with Vn = mean of gradient crossproducts
  ##       - RDM2: sum(grad * dtheta2) / ntheta
  ##    * theta: starting value
  ##    * grad.tol: tolerance value for L2norm(grad)
  ##    * RDM.tol: tolerance value for the RDM (= Relative Damping Measure)
  ##    * fun.tol: tolerance value for changes in function g to declare convergence
  ##    * itermax: maximum number of iterations
  ##
  ## OUTPUT:
  ##    List containing:
  ##       * val : function value at <theta>
  ##       * val.start: starting value of the function at theta0
  ##       * theta: value of theta maximizing the function (if convergence occured)
  ##       * grad: gradient at <theta>
  ##       * RDM: Relative Damping Measure at <theta>
  ##       * iter: number of iterations
  ##       * convergence: convergence indicator based on L2(grad), RDM and last changes in g(.)
  NewtonRaphson <- function(g, theta, grad.tol=1e-2, RDM.tol=1e-4, fun.tol=1e-3,
                            itermax=100, step_factor=10, verbose=FALSE){
      ntheta = length(theta)
      theta.cur = theta
      obj.cur = g(theta.cur,Dtheta=TRUE)
      g.start = obj.cur$g ## Function at the iteration start
      ## Euclidean norm function
      L2norm <- function(x) sqrt(sum(x^2))
      ## Newton-Raphson (NR) or Gradient descent (GD) step
      NRGDstep <- function(obj.cur, step, method=c("NR","NR2","GD")){
          method = match.arg(method)
          if (method == "NR") theta.prop = with(obj.cur, theta + step * dtheta)
          if (method == "NR2") theta.prop = with(obj.cur, theta + step * dtheta2)
          if (method == "GD"){
              grad.norm = with(obj.cur, grad / sqrt(L2norm(grad)))
              theta.prop = with(obj.cur, theta + step * grad.norm)
          }
          obj.prop = tryCatch(expr=g(theta.prop,Dtheta=TRUE), error=function(e) e)
          if (inherits(obj.prop, "error") || !("g" %in% names(obj.prop)) || is.na(obj.prop$g)) {
              return(list(better = FALSE, gain = 0, obj=obj.cur, step=step/step_factor, method=method))
          }
          gain = obj.prop$g - obj.cur$g
          if (gain > 0){ ## Accept obj.prop
              step = step_factor * step ## 1.1*step
              obj = obj.prop
          } else { ## Reject obj.prop
              step = step / step_factor
              obj = obj.cur
          }
          return(list(better = gain > 0, gain = gain, obj=obj, step=step, method=method))
      }
      ## Convergence criterion using RDM (see e.g. Prague et al, 2013)
      ##    RDM < eps.RMD ?   where  RDM = grad' (-H)^-1 grad / ntheta
      RDM = obj.cur$RDM ; RDM2 = obj.cur$RDM2
      ## grad = obj.cur$grad
      ## dtheta = obj.cur$dtheta
      ## RDM = sum(grad * dtheta) / ntheta
      ##
      iter = 0
      grad.desc = TRUE ## Also test gradient descent step
      converged = (L2norm(obj.cur$grad) < grad.tol)
      ## cat("\nBefore: range(grad):",range(obj.cur$grad),"\n")
      ## converged = (RDM < RDM.tol) && all(abs(grad) < grad.tol)
      MaxStepHalving = 15
      gain = 0 ; nrep = 0
      ## Main loop
      ## ---------
      step1 = step2 = 1
      step.GD = .1
      while(!(converged || (iter >= itermax) || (nrep >= MaxStepHalving))) {
          iter = iter + 1
          ##
          action = "Reject" ;
          g0 = obj.cur$g ## Starting value of the function
          gain = 0 ## Monitoring function change
          nrep = 0
          ## -- Step-halving loop --
          for (nrep in 1:MaxStepHalving){
              ## -- Newton-Raphson step --
              if (obj.cur$RDM > 0){
                  temp = NRGDstep(obj.cur, step1, "NR")
                  obj.cur = temp$obj ; step1 = temp$step
                  if (temp$better){
                      gain = temp$gain
                      action = "NR"
                      break
                  }
              } else {
                  ## -- Newton-Raphson step using score-based Vn matrix --
                  temp = NRGDstep(obj.cur, step2, "NR2")
                  obj.cur = temp$obj ; step2 = temp$step
                  if (temp$better){
                      gain = temp$gain
                      action = "NR2"
                      break
                  } else {
                      ## -- Gradient descent step --
                      temp = NRGDstep(obj.cur, step.GD, "GD")
                      obj.cur = temp$obj ; step.GD = temp$step
                      if (temp$better){
                          gain = temp$gain
                          action = "GD"
                          break
                      }
                  }
              }
          } ## End Step-halving
          ##
          ## Convergence criterion using RDM (see e.g. Prague et al, 2013)
          ##    RDM < eps.RMD ?   where  RDM = grad' (-H)^-1 grad / ntheta
          grad = c(obj.cur$grad)
          RDM = obj.cur$RDM ; RDM2 = obj.cur$RDM2
          ## RDM = with(obj.cur, sum(grad * dtheta)) / ntheta
          conv.RDM = (0 < RDM & RDM < RDM.tol) || (RDM2 < RDM.tol)
          conv.grad = (L2norm(grad) < grad.tol) ## all(abs(grad) < grad.tol)
          conv.gain = (gain > 0 & gain < fun.tol)
          converged = (conv.RDM && conv.gain) && conv.grad
          if (verbose){
              if (action == "NR") cat(action,"step =",step1," RDM=",RDM,g0,obj.cur$g,gain,"\n")
              if (action == "NR2") cat(action,"step =",step2," RDM2=",RDM2,g0,obj.cur$g,gain,"\n")
              if (action == "GD") cat(action,"step =",step.GD,g0,obj.cur$g,gain,"\n")
          }
      } ## End <while> (main loop)
      ## cat("\nAfter: range(grad):",range(obj.cur$grad),"\n")
      ans = list(val=obj.cur$g, val.start=g.start,
                 theta=obj.cur$theta,
                 grad=obj.cur$grad, RDM=RDM,
                 iter=iter, converged=converged)
      return(ans)
  } ## End NewtonRaphson
  ##
  ## ---------------------------------------
  ##  Generic Levenberg-Marquardt algorithm
  ## ---------------------------------------
  ## GOAL: maximize function using its gradient and Hessian with Levenberg-Marquardt
  ## INPUT:
  ##    * g: function returning a list:
  ##       - g: function value to be maximized
  ##       - grad: gradient of the function
  ##       - Hes: Hessian matrix of the function
  ##       - RDM: grad' (-H)^-1 grad / ntheta
  ##    * theta0: starting value
  ##    * grad.tol: tolerance value for the L2-norm of the gradient
  ##    * RDM.tol: tolerance value for the RDM (= Relative Damping Measure)
  ##    * max_iter: maximum number of iterations
  ##    * lambda_init: initial value of the damping parameter in Levenberg-Marquardt
  ##    * lambda_factor: adaptative factor for the damping parameter
  ##
  ## OUTPUT:
  ##    List containing:
  ##       * val : function value at <theta>
  ##       * val.start: starting value of the function at theta0
  ##       * theta: value of theta maximizing the function (if convergence occured)
  ##       * grad: gradient at <theta>
  ##       * RDM: Relative Damping Measure at <theta>
  ##       * iter: number of iterations
  ##       * convergence: logical indicating if convergence occured: L2(grad) < grad.tol
  myLevMarq <- function(g, theta0, grad.tol=1e-2, RDM.tol=1e-4,
                        max_iter=25, lambda_init=1e-3, lambda_factor=10,
                        verbose=FALSE){
      theta <- theta0
      ntheta <- length(theta)
      lambda <- lambda_init
      L2norm <- function(x) sqrt(sum(x^2)) ## Euclidean norm function
      obj.cur <- g(theta0, Dtheta=TRUE)
      g.start = obj.cur$g
      ##
      nreject = 0
      for (i in 1:max_iter) {
          grad <- obj.cur$grad
          hess <- -obj.cur$Hes ## hess = NEGATIVE Hessian !!
          ## Hessian regularisation
          hess_reg <- hess ; diag(hess_reg) <- diag(hess_reg) + lambda
          ##
          ## Check positive definite
          is_positive_definite <- FALSE
          while (!is_positive_definite) {
              tryCatch({
                  chol(hess_reg)
                  is_positive_definite <- TRUE
              }, error = function(e) {
                  ## Increase lambda if not positive-definite
                  lambda <<- lambda * lambda_factor
                  hess_reg <<- hess ; diag(hess_reg) <<- diag(hess_reg) + lambda
              })
          }
          ## Compute dtheta and update theta
          dtheta <- solve(hess_reg, grad)
          theta.prop <- theta + dtheta ## Update theta
          ##
          ## Accept update ?
          obj.prop <- g(theta.prop, Dtheta=TRUE)
          g.cur = obj.cur$g ; g.prop = obj.prop$g
          g.dif = g.prop - g.cur
          if (verbose) cat(i,"> ",lambda, g.cur, g.prop, "--> g.dif =",g.dif,"\n")
          if (!is.na(g.dif) && (g.dif > 0)) {
              nreject = 0
              theta <- theta.prop
              lambda <- lambda / lambda_factor
              obj.cur <- obj.prop
          } else {
              nreject = nreject + 1
              if (nreject > 10) break
              lambda <- lambda * lambda_factor
          }
          if (verbose) cat("L2(grad):",L2norm(obj.cur$grad),"\n")
          ## Compute RDM
          grad = obj.cur$grad
          RDM = sum(grad * dtheta) / ntheta
          ## Check convergence
          converged = (L2norm(obj.cur$grad) < grad.tol)
          if (converged) break
          ## if ((converged) || (abs(g.dif) < fun.tol)) break
      } ## End iteration loop
      ans = list(val=obj.cur$g, val.start=g.start,
                 theta=obj.cur$theta,
                 grad=obj.cur$grad, RDM=RDM,
                 iter=i, converged=converged)
      return(ans)
  } ## End myLevMarq
  ##
  ## ///////////////////////////////////////
  ## ----> BEGIN of the Estimation algorithm
  ## ///////////////////////////////////////
  ##
  converged = FALSE
  iter = 0
  ptm <- proc.time()
  ##
  KK = K.error ## Number of spline parameters for the (log-hazard of the) stdzd error density
  ##
  ## Spline parameters for (the log-hazard of) a N(0,1) density
  phi.Norm <- function(knots){
    xg = seq(min(knots),max(knots),length=201)
    Bg = Bsplines(xg, knots) ## B-spline basis on a fine grid
    lhaz.norm = dnorm(xg,log=TRUE) - pnorm(xg,log.p=TRUE,lower.tail=FALSE)
    phi.cur = c(solve(t(Bg)%*%Bg+diag(1e-6,ncol(Bg)))%*%t(Bg)%*%lhaz.norm)
    return(phi=phi.cur)
  }
  ##
  final.iteration = FALSE ## Idea: when the EDs of the additive terms stabilize,
  ## perform a final iteration to get the final estimates for the regression parms <psi1> & <psi2>
  while (!converged){
    iter = iter + 1
    if (verbose){
        cat("\n--------------\n")
        cat("Iteration ",iter,"\n")
        cat("--------------\n")
    }
    psi1.old = psi1.cur ; psi2.old = psi2.cur ; psi.old = c(psi1.old,psi2.old)
    ## Evaluate fitted location and dispersion
    obj1d = NULL ; res = res.psi.fun(resp,c(psi1.cur,psi2.cur))
    ##
    ## -1- Update B-spline estimation of the error density
    ## ---------------------------------------------------
    if (!final.iteration){
        ## Estimate error density
        if (!fixed.support) {
            rmin = -6 ; rmax = 6
        }
        obj1d = Dens1d(res, event=event, w=w,
                       ymin=rmin, ymax=rmax, K=K.error,
                       equid.knots=TRUE)
        if (!is.null(phi.0)){
            phi.cur = phi.0
            phi.ref = NULL
            fixed.phi = FALSE
        } else {
            if ((iter <= 2) | Normality) { ## Start with Normal error...
                phi.ref = phi.Norm(obj1d$knots)
                phi.cur = phi.ref
                fixed.phi = TRUE
            } else { ## ... leave it free later on
                phi.ref = NULL
                fixed.phi = FALSE
            }
        }
        ## Start by not constraining moments for the error density...
        fit1d = fit1d.0 = densityLPS(obj1d,
                                     tau0=lambda0.cur, phi0=phi.cur,
                                     fixed.phi=fixed.phi,phi.ref=phi.ref,
                                     is.density=TRUE,Mean0=NULL,Var0=NULL,
                                     method=density.method,verbose=FALSE)
        ## ... then constrain moments if necessary
        moments.cur = with(fit1d.0, c(mean.dist,var.dist))
        if (!fixed.phi & (max(moments.cur-c(0,1)) > 1e-3)){
            fit1d = densityLPS(obj1d,
                               tau0=fit1d.0$tau,phi0=fit1d.0$phi,
                               fixed.phi=fixed.phi, phi.ref=phi.ref,
                               is.density=TRUE,Mean0=0,Var0=1,
                               method=density.method,verbose=FALSE)
        }
        ## Sample size by censoring status
        fit1d$n=n ; fit1d$n.uncensored=n.uncensored ; fit1d$n.IC=n.IC ; fit1d$n.RC=n.RC
        ## Current penalty parameter for the spline parameters in the error density
        lambda0.cur = fit1d$tau
        ## Current spline parameter estimates in error density estimation
        phi.cur = fit1d$phi
    }
    ## -2- Update <psi1> & <psi2> (i.e. regression coefs related to location)
    ## --------------------------
    ## Initialize obj.cur (lpost, Score, Hessian, etc.)
    obj.cur = ff(c(psi1.cur,psi2.cur),lambda1.cur,lambda2.cur,REML=REML,diag.only=diag.only)
    ED1.old = obj.cur$ED1 ; ED2.old = obj.cur$ED2
    ##
    ## -- Estimation of (psi1,psi2) given (lambda1,lambda2) --
    if (verbose) cat("\n-- Update regression and spline parameters --\n\n")
    if (psi.method == "LM"){
        ## -- LM: Levenberg-Marquardt (DEFAULT) --
        g.regr <- function(theta,Dtheta=TRUE){
            ntheta = length(theta)
            obj.cur = ff(theta,lambda1.cur,lambda2.cur,REML=REML,diag.only=diag.only,
                         hessian=Dtheta)
            grad = Hes = dtheta = NULL
            if (Dtheta){
                grad = obj.cur$U.psi
                Hes = obj.cur$Hes.psi
            }
            ##
            ans = list(g=obj.cur$lpost, theta=theta, grad=grad, Hes=Hes)
            return(ans)
        } ## End g.regr
        ##
        regr.LM = myLevMarq(g=g.regr, theta=psi.cur,
                            grad.tol=grad.tol, RDM.tol=RDM.tol,
                            verbose=verbose)
        if (verbose){
            with(regr.LM,
                 cat("\nLevenberg-Marquardt: niter=",iter,
                     " ; RDM=",RDM,
                     " ; L2(grad)=",L2norm(grad),
                     " ; converged=",converged,
                     "\n",sep=""))
        }
        iter.psi = regr.LM$iter
        psi.cur = regr.LM$theta
        psi1.cur = psi.cur[1:q1] ; psi2.cur = psi.cur[-(1:q1)]
        obj.cur = ff(psi.cur,lambda1.cur,lambda2.cur,REML=REML,diag.only=diag.only)
        ok.psi = ok.psi1 = ok.psi2 = regr.LM$converged
        ## End LM
    } else if (psi.method == "NR"){
        ## -- NR: Newton-Raphson --
        g.regr <- function(theta,Dtheta=TRUE){
            ntheta = length(theta)
            obj.cur = ff(theta,lambda1.cur,lambda2.cur,REML=REML,diag.only=diag.only,
                         hessian=Dtheta)
            grad = Hes = dtheta = RDM = NULL
            if (Dtheta){
                grad = obj.cur$U.psi
                Hes = obj.cur$Hes.psi
                A = -obj.cur$Hes.psi + diag(1e-6,length(grad))
                dtheta = solve(A, grad)
                ##
                A2 = obj.cur$Vn.psi
                dtheta2 = solve(A2, grad)
                ##
                RDM = sum(grad * dtheta) / ntheta
                RDM2 = sum(grad * dtheta2) / ntheta
            }
            ##
            ans = list(g=obj.cur$lpost, theta=theta,
                       grad=grad, ## Hes=Hes,
                       dtheta=dtheta, RDM=RDM,
                       dtheta2=dtheta2, RDM2=RDM2
                       )
            return(ans)
        } ## End g.regr
        ##
        regr.NR = NewtonRaphson(g=g.regr,theta=psi.cur,
                                grad.tol=grad.tol,RDM.tol=RDM.tol,fun.tol=fun.tol,
                                verbose=verbose)
        if (verbose){
            cat("\nNewton-Raphson: niter=",regr.NR$iter,
                " ; RDM=",regr.NR$RDM,
                " ; converged=",regr.NR$converged,"\n",sep="")
        }
        iter.psi = regr.NR$iter
        psi.cur = regr.NR$theta
        psi1.cur = psi.cur[1:q1] ; psi2.cur = psi.cur[-(1:q1)]
        obj.cur = ff(psi.cur,lambda1.cur,lambda2.cur,REML=REML,diag.only=diag.only)
        ok.psi = ok.psi1 = ok.psi2 = regr.NR$converged
        ## End NR
    } else if (psi.method == "BFGS"){
        ## -- BFGS --
        fn <- function(psi) -ff(psi, lambda1.cur,lambda2.cur, hessian=FALSE,
                               REML=REML, diag.only=diag.only)$lpost
        gr <- function(psi) -ff(psi, lambda1.cur,lambda2.cur, hessian=FALSE,
                               REML=REML, diag.only=diag.only)$U.psi
        result = optim(psi.cur, fn=fn, gr=gr, method="BFGS")
        iter.psi = NA
        psi.cur = result$par  ; psi1.cur = psi.cur[1:q1] ; psi2.cur = psi.cur[-(1:q1)]
        obj.cur = ff(psi.cur,lambda1.cur,lambda2.cur,REML=REML,diag.only=diag.only)
        obj.cur$U.psi1 = obj.cur$U.psi[1:q1] ; obj.cur$U.psi2 = obj.cur$U.psi[-(1:q1)]
        ok.psi = ok.psi1 = ok.psi2 = (result$convergence == 0)
        ## End BFGS
    } else if (psi.method == "NM"){
        ## -- Nelder-Mead --
        fn <- function(psi) -ff(psi, lambda1.cur,lambda2.cur,
                               gradient=FALSE, hessian=FALSE,
                               REML=REML, diag.only=diag.only)$lpost
        gr <- function(psi) -ff(psi, lambda1.cur,lambda2.cur,
                               gradient=TRUE, hessian=FALSE,
                               REML=REML, diag.only=diag.only)$U.psi
        result = optim(psi.cur, method="Nelder-Mead", fn=fn, gr=gr)
        iter.psi = result$iter
        psi.cur = result$par  ; psi1.cur = psi.cur[1:q1] ; psi2.cur = psi.cur[-(1:q1)]
        obj.cur = ff(psi.cur,lambda1.cur,lambda2.cur,REML=REML,diag.only=diag.only)
        obj.cur$U.psi1 = obj.cur$U.psi[1:q1] ; obj.cur$U.psi2 = obj.cur$U.psi[-(1:q1)]
        ok.psi = ok.psi1 = ok.psi2 = (result$convergence > 0)
        ## End Nelder-Mead
    }
    ## -----------------------------------------------
    ##  Estimation methods based on external packages
    ## -----------------------------------------------
    ## } else if (psi.method == "LevMarq"){ ## Requires <MarqLevAlg> R-package
    ##     ## -- Levenberg-Marquardt (fn + gradient + Hessian) --
    ##     fn <- function(psi) -ff(psi, lambda1.cur,lambda2.cur, hessian=FALSE,
    ##                            REML=REML, diag.only=diag.only)$lpost
    ##     gr <- function(psi) -ff(psi, lambda1.cur,lambda2.cur, hessian=FALSE,
    ##                            REML=REML, diag.only=diag.only)$U.psi
    ##     hess <- function(psi) -ff(psi, lambda1.cur,lambda2.cur, hessian=TRUE,
    ##                              REML=REML, diag.only=diag.only)$Hes.psi
    ##     result = marqLevAlg(b=psi.cur, fn=fn, gr=gr, hess=hess, epsb=fun.tol, epsd=RDM.tol)
    ##     iter.psi = result$ni
    ##     psi.cur = result$b  ; psi1.cur = psi.cur[1:q1] ; psi2.cur = psi.cur[-(1:q1)]
    ##     obj.cur = ff(psi.cur,lambda1.cur,lambda2.cur,REML=REML,diag.only=diag.only)
    ##     obj.cur$v = result$v ## Upper diagonal elements of Variance matrix
    ##     obj.cur$U.psi = result$grad
    ##     obj.cur$U.psi1 = obj.cur$U.psi[1:q1] ; obj.cur$U.psi2 = obj.cur$U.psi[-(1:q1)]
    ##     ok.psi = ok.psi1 = ok.psi2 = (result$rdm <= RDM.tol)
    ##     ## End LevMarq
    ## }
    ## else if (psi.method == "nloptr"){ ## Requires <nloptr> R-package
    ##     ## -- nloptr --
    ##     ff.opt.lst <- function(psi){
    ##         temp = ff(psi, lambda1.cur,lambda2.cur, hessian=FALSE, REML=REML, diag.only=diag.only)
    ##         list("objective" = -temp$lpost,
    ##              "gradient" = -temp$U.psi)
    ##     }
    ##     opts = list("algorithm" = "NLOPT_LD_LBFGS", "xtol_rel" = 1e-6, "maxeval" = 500, print.level=3)
    ##     result = nloptr(x0=psi.cur, eval_f=ff.opt.lst, opts=opts)
    ##     iter.psi = result$iterations
    ##     Psi.cur = result$solution ; psi1.cur = psi.cur[1:q1] ; psi2.cur = psi.cur[-(1:q1)]
    ##     obj.cur = ff(psi.cur,lambda1.cur,lambda2.cur,REML=REML,diag.only=diag.only)
    ##     obj.cur$U.psi = -nl.grad(psi.cur, function(x) ff.opt.lst(x)$objective)
    ##     obj.cur$U.psi1 = obj.cur$U.psi[1:q1] ; obj.cur$U.psi2 = obj.cur$U.psi[-(1:q1)]
    ##     ok.psi = ok.psi1 = ok.psi2 = (obj.cur$status == 0) ## TRUE
    ##     ## End nloptr
    ## } else if (psi.method == "LevMarq1"){ ## Requires <MarqLevAlg> R-package
    ##     ## -- Levenberg-Marquardt (fn + gradient) --
    ##     fn <- function(psi) -ff(psi, lambda1.cur,lambda2.cur, hessian=FALSE,
    ##                            REML=REML, diag.only=diag.only)$lpost
    ##     gr <- function(psi) -ff(psi, lambda1.cur,lambda2.cur, hessian=FALSE,
    ##                            REML=REML, diag.only=diag.only)$U.psi
    ##     result = marqLevAlg(b=psi.cur, fn=fn, gr=gr, epsb=fun.tol, epsd=RDM.tol)
    ##     iter.psi = result$ni
    ##     psi.cur = result$b  ; psi1.cur = psi.cur[1:q1] ; psi2.cur = psi.cur[-(1:q1)]
    ##     obj.cur = ff(psi.cur,lambda1.cur,lambda2.cur,REML=REML,diag.only=diag.only)
    ##     obj.cur$v = result$v ## Upper diagonal elements of Variance matrix
    ##     obj.cur$U.psi = result$grad
    ##     obj.cur$U.psi1 = obj.cur$U.psi[1:q1] ; obj.cur$U.psi2 = obj.cur$U.psi[-(1:q1)]
    ##     ok.psi = ok.psi1 = ok.psi2 = (result$rdm <= RDM.tol)
    ##     ## End LevMarq1
    ## } else if (psi.method == "nlm"){
    ##     ## -- nlm (fn + gradient) --
    ##     ff.opt <- function(psi){
    ##         temp = ff(psi, lambda1.cur,lambda2.cur, hessian=FALSE, REML=REML, diag.only=diag.only)
    ##         ans = -temp$lpost
    ##         attr(ans,'gradient') = -temp$U.psi
    ##         return(ans)
    ##     }
    ##     result = nlm(ff.opt, psi.cur, iterlim=1000,hessian=TRUE)
    ##     iter.psi = result$iterations
    ##     psi.cur = result$est  ; psi1.cur = psi.cur[1:q1] ; psi2.cur = psi.cur[-(1:q1)]
    ##     obj.cur = ff(psi.cur,lambda1.cur,lambda2.cur,REML=REML,diag.only=diag.only)
    ##     obj.cur$U.psi = result$gradient
    ##     obj.cur$U.psi1 = obj.cur$U.psi[1:q1] ; obj.cur$U.psi2 = obj.cur$U.psi[-(1:q1)]
    ##     ok.psi = ok.psi1 = ok.psi2 = (result$code <= 3)
    ## End nlm
    ## } else if (psi.method == "nlm2"){
    ##     ## -- nlm2 (fn + gradient + Hessian) --
    ##     ff.opt <- function(psi){
    ##         temp = ff(psi, lambda1.cur,lambda2.cur, hessian=TRUE, REML=REML, diag.only=diag.only)
    ##         ans = -temp$lpost
    ##         attr(ans,'gradient') = -temp$U.psi
    ##         attr(ans,'hessian') = -temp$Hes.psi
    ##         return(ans)
    ##     }
    ##     result = nlm(ff.opt, psi.cur, iterlim=1000, hessian=TRUE, check.analyticals=FALSE)
    ##     psi.cur = result$est  ; psi1.cur = psi.cur[1:q1] ; psi2.cur = psi.cur[-(1:q1)]
    ##     obj.cur = ff(psi.cur,lambda1.cur,lambda2.cur,REML=REML,diag.only=diag.only)
    ##     obj.cur$U.psi = result$gradient
    ##     obj.cur$U.psi1 = obj.cur$U.psi[1:q1] ; obj.cur$U.psi2 = obj.cur$U.psi[-(1:q1)]
    ##     ok.psi = ok.psi1 = ok.psi2 = TRUE
    ##     ## End nlm2
    ## }
    ##
    ## -3- Update <lambda1> & <lambda2> (i.e. penalty vectors for additives terms in loc and disp)
    ## --------------------------------
    itermin = 3 ## Only start updating lambda's after <itermin> iterations
    update.lambda = ifelse(iter <= itermin, FALSE, (J1 > 0)|(J2 > 0))
    ## res = obj.cur$res
    ## browser()
    ##
    if (update.lambda & !final.iteration){
        if (verbose) cat("\n-- Update penalty parameters --\n")
        ## -- Method 1 (DEFAULT) --
        if (lambda.method == "LPS"){
            temp = select.lambda.LPS(psi=psi.cur, lambda1=lambda1.cur,lambda2=lambda2.cur)
        }
        ## -- Method 2 --
        if (lambda.method == "LPS2"){
            temp = select.lambda.LPS2(psi=psi.cur, lambda1=lambda1.cur,lambda2=lambda2.cur)
        }
        ## ## -- Method 2b --
        ## if (lambda.method == "LPS2b"){
        ##     temp = select.lambda.LPS2b(psi=psi.cur, lambda1=lambda1.cur,lambda2=lambda2.cur)
        ## }
        ##
        ## -- Method 3: fixed penalty parameters --
        if (lambda.method == "none"){
            lambda1.mat = cbind(lambda=lambda1.cur,low=lambda1.cur,up=lambda1.cur)
            lambda2.mat = cbind(lambda=lambda2.cur,low=lambda2.cur,up=lambda2.cur)
        } else {
            ## Location
            lambda1.cur = temp$lambda1
            lambda1.mat = temp$lambda1.mat
            ## se.lambda1 = temp$se.lambda1 ; se.loglambda1 = temp$se.loglambda1
            ##
            ## Dispersion
            lambda2.cur = temp$lambda2
            lambda2.mat = temp$lambda2.mat
            ##se.lambda2 = temp$se.lambda2 ; se.loglambda2 = temp$se.loglambda2
        }
        ## -- Reevaluate the log-posterior and associated quantities --
        obj.cur = ff(psi.cur,lambda1.cur,lambda2.cur,REML=REML,diag.only=diag.only)
    }
    ## Global Stopping rule
    ## --------------------
    L2 <- function(x) sqrt(sum(x^2))
    L1 <- function(x) sum(abs(x))
    ##
    eps = grad.tol
    ##
    if (!final.iteration){
      converged.loc = (L2(obj.cur$U.psi1) < eps)   ## L2 criterion
      converged.disp = (L2(obj.cur$U.psi2) < eps)  ## L2 criterion
      ## converged.loc = (max(abs(obj.cur$U.psi1)) < eps)   ## Linf criterion
      ## converged.disp = (max(abs(obj.cur$U.psi2)) < eps)  ## Linf criterion
      ##
      ## Check if the EDs of the additive terms stabilize
      if (J1 > 0) converged.ED1 = ifelse(iter<=itermin,
                                         FALSE,
                                         ((L1(obj.cur$ED1-ED1.old) / L1(ED1.old)) < eps))
      else converged.ED1 = TRUE
      if (J2 > 0) converged.ED2 = ifelse(iter<=itermin,
                                         FALSE,
                                         ((L1(obj.cur$ED2-ED2.old) / L1(ED2.old)) < eps))
      else converged.ED2 = TRUE
    }
    converged.detail = c(psi1=ok.psi1,psi2=ok.psi2,ED1=converged.ED1,ED2=converged.ED2)
    ok.psi = (ok.psi1 & ok.psi2)
    ## When the EDs of the additive terms stabilize,
    ##  perform a final iteration to get the final estimates for the regression parms <psi1> & <psi2>
    green = all(c(ok.psi,converged.ED1,converged.ED2)) ## Necessary condition to complete the estimation process
    ## Necessary & sufficient condition to complete the estimation process: green light followed by one last iteration
    if (final.iteration & green) converged = TRUE
    final.iteration = green ## One final iteration required after the 1st green light to complete the estimation process
    ##
    if (verbose){
      if (J1 > 0){
        cat("\nlambda1.cur:",lambda1.cur,"\n")
        ## cat("ED1 before & after: ",ED1.old," ---> ",obj.cur$ED1,"\n")
        cat("ED1\n")
        print(cbind(before=ED1.old,after=obj.cur$ED1))
      }
      if (J2 > 0){
        cat("\nlambda2.cur:",lambda2.cur,"\n")
        ## cat("ED2 before & after: ",ED2.old," ---> ",obj.cur$ED2,"\n")
        cat("ED2:\n")
        print(cbind(before=ED2.old,after=obj.cur$ED2))
      }
      ## cat("\niter:",iter," (",iter.psi,iter.lam1,iter.lam2,") ; convergence:",converged.detail,"\n")
      cat("\nScore for <psi1> in",range(obj.cur$U.psi1),"\n")
      cat("Score for <psi2> in",range(obj.cur$U.psi2),"\n")
      ##
      cat("\nConvergence details:\n")
      print(converged.detail)
    }
    ##
    if (iter > iterlim){
        cat("Algorithm did not converge after",iter,"iterations\n")
        break
        fit = list(converged=FALSE)
        return(fit) ## Stop the algorithm before convergence
    }
  }
  elapsed.time <- (proc.time()-ptm)[1] ## Elapsed time
  ## ///////////////////////////////////////
  ## <---- End of the Estimation algorithm
  ## ///////////////////////////////////////
  ##
  ## Sandwich variance-covariance for <psi>
  ## --------------------------------------
  temp = cbind(Xcal.1 * (sqrt(w) * obj.cur$omega.psi1), Xcal.2 * (sqrt(w) * obj.cur$omega.psi2))
  Vn.psi = crossprod(temp) ## Vn.psi should be multiplied by <k> if (w -> k*w)  (hence sqrt(w) !!)
  ## ## Meaning of crossprod:
  ## Vn.psi1 = temp[1,]%o%temp[1,]
  ## for (i in 2:n) Vn.psi1 = Vn.psi1 + temp[i,]%o%temp[i,]
  ## ## End meaning
  idx1 = 1:q1 ; idx2 = (q1+1):(q1+q2)
  if (J1 > 0) Vn.psi[idx1,idx1] = Vn.psi[idx1,idx1] + obj.cur$P1
  if (J2 > 0) Vn.psi[idx2,idx2] = Vn.psi[idx2,idx2] + obj.cur$P2
  ##
  if (REML){  ## ATTENTION: sign correction needed:  MINUS dHes.REML & ddiagHes.REML !!
    if (diag.only) diag(Vn.psi)[idx2] = diag(Vn.psi)[idx2] - obj.cur$ddiagHes.REML
    else Vn.psi[idx2,idx2] = Vn.psi[idx2,idx2] - obj.cur$dHes.REML
  }
  ##
  bread.psi = obj.cur$Cov.psi ## Minus-Inverse Hessian for <psi>
  Sand.psi = bread.psi %*% (Vn.psi %*% bread.psi) ## Sandwich estimator
  ##
  ## Compute variance-covariance matrices for location and dispersion parameters
  Sand.psi1  = Sand.psi[idx1,idx1,drop=FALSE] ## Sandwich estimator for Var-Cov of location regr coef
  Sand.psi2  = Sand.psi[idx2,idx2,drop=FALSE] ## Sandwich estimator for Var-Cov of dispersion regr coef
  bread.psi1 = bread.psi[idx1,idx1,drop=FALSE] ## Minus-inverse Hessian part for location regr coef
  bread.psi2 = bread.psi[idx2,idx2,drop=FALSE] ## Minus-inverse Hessian part for ion regr coef
  if (!Normality & !sandwich){
    Cov.psi1 = bread.psi1 ## Model-based estimation for Var-Cov of location regr coef
    Cov.psi2 = bread.psi2 ## Model-based estimation for Var-Cov of dispersion regr coef
  } else {
    Cov.psi1 = Sand.psi1 ## Sandwich estimator for Var-Cov of location regr coef
    Cov.psi2 = Sand.psi2 ## Sandwich estimator for Var-Cov of dispersion regr coef
  }
  ## Report on fixed effects
  ## -----------------------
  ## ... for Location
  beta.hat = psi1.cur[1:nfixed1] ; names(beta.hat)=colnames(Xcal.1)[1:nfixed1]
  beta.se = sqrt(pmax(0,diag(Cov.psi1[1:nfixed1,1:nfixed1,drop=FALSE])))
  low = beta.hat-z.alpha*beta.se ; up = beta.hat+z.alpha*beta.se
  z.beta = beta.hat / beta.se ; Pval.beta = 1-pchisq(z.beta^2,1)
  mat.beta = cbind(est=beta.hat,se=beta.se,low=low,up=up,Z=z.beta,Pval=Pval.beta)
  ## ... for Dispersion
  delta.hat = psi2.cur[1:nfixed2] ; names(delta.hat)=colnames(Xcal.2)[1:nfixed2]
  delta.se = sqrt(diag(Cov.psi2[1:nfixed2,1:nfixed2,drop=FALSE]))
  low = delta.hat-z.alpha*delta.se ; up = delta.hat+z.alpha*delta.se
  z.delta = delta.hat / delta.se ; Pval.delta = 1-pchisq(z.delta^2,1)
  mat.delta = cbind(est=delta.hat,se=delta.se,low=low,up=up,Z=z.delta,Pval=Pval.delta)
  ##
  ## Report on effective dimension of additive components / Test Horiz / Test Linear
  ## -------------------------------------------------------------------------------
  if (J1 > 0){
      ## xi1.se = sqrt(pmax(1e-6,diag(MASS::ginv(-Hes.xi1)))) ## sqrt(diag(solve(-Hes.xi1)))
      ## xi1.low = xi1.cur - z.alpha*xi1.se ; xi1.up = xi1.cur + z.alpha*xi1.se
      lambda1.low = lambda1.mat[,"low"] ; lambda1.up = lambda1.mat[,"up"]
      ## lambda1.low = lambda1.min + exp(xi1.low) ; lambda1.up = lambda1.min + exp(xi1.up)
      ##
      ## Elements for testing absence of effect (horizontality) or linear effect for additive terms
      D1=diff(diag(K1),dif=1) ; D2=diff(diag(K1),dif=2) ## Diff matrices of order 1 & 2
      ## D1=diff(diag(K1+1),dif=1)[,-(K1+1)] ; D2=diff(diag(K1+1),dif=2)[,-(K1+1)] ## Diff matrices of order 1 & 2 for centered B-spline matrices without their last basis component for identification purposes in additive models
      psi = psi1.cur ## Concerned spline coefs
      Cov = Cov.psi1 ## solve(-Hes.psi1) ## Var-Cov matrix of location regr coef  ## HERE HERE
      Chi2.horiz1 = Chi2.linear1 = Pval.horiz1 = Pval.linear1 = 0*ED1.cur ## Test for Horizontality and Linearity of additive terms
      ED1.low = ED1.up = 0*ED1.cur
      BWB.1 = with(obj.cur, -(Hes.psi1+P1)) ## gives  t(Xcal.1)%*%(W.psi1*Xcal.1)
      for (j in 1:J1){
          ## Horizontal ?
          ## ----------
          idx = nfixed1 + (j-1)*K1 + (1:K1)
          Chi2.horiz1[j] = sum(psi[idx] * c(solve(Cov[idx,idx]) %*% c(psi[idx]))) ## Test whether can set all spline parms to 0 ; Indeed, with a centered B-spline basis, the last spline coef is not estimated and set to zero
          Pval.horiz1[j] = 1-pchisq(Chi2.horiz1[j],obj.cur$ED1[j])
          ##
          ## Linear ?
          ## ------
          Chi2.linear1[j] = sum(c(D1%*%psi[idx]) * c(solve(D1%*%Cov[idx,idx]%*%t(D1)) %*% c(D1%*%psi[idx]))) ##
          Pval.linear1[j] = 1-pchisq(Chi2.linear1[j],obj.cur$ED1[j]-1)
          ##
          ## Lower and upper value for effective degrees of freedom
          ## ------------------------------------------------------
          lam1.low = lam1.up = lambda1.cur
          ## Only set the jth component of lambda1 to lower or upper CI value
          lam1.low[j] = lambda1.low[j] ; lam1.up[j] = min(lambda1.up[j],1e5)
          ##
          ## !! Lower (upper) bound for lambda corresponds to upper (lower) bound for EDF !!
          ED1.up[j]= EDF.fun(obj.cur$Hes0.psi, lam1.low, lambda2.cur,
                             parms1, parms2, joint.computation=joint.computation)$ED1[j]
          ED1.low[j]= EDF.fun(obj.cur$Hes0.psi, lam1.up, lambda2.cur,
                              parms1, parms2, joint.computation=joint.computation)$ED1[j]
          ## ## Begin Old code -->
          ## P.low = Pcal.fun(nfixed1,lam1.low,Pd1.x) ; P.up = Pcal.fun(nfixed1,lam1.up,Pd1.x)
          ## temp.low = diag(solve(BWB.1 + P.low) %*% BWB.1) ## Effective dim for each regression parameter in location
          ## temp.up  = diag(solve(BWB.1 + P.up ) %*% BWB.1) ##  with only jth penalty set at its lower or upper CI value
          ## idx = nfixed1 + (j-1)*K1 + (1:K1) ## Concerned diagonal elements for jth additive term
          ## ## Effective dim of the jth additive term in location with lambda.j at its lower CI value
          ## ED1.up[j]  = sum(temp.low[idx])
          ## ## Effective dim of the jth additive term in location with lambda.j at its upper CI value
          ## ED1.low[j] = sum(temp.up[idx])
          ## ## <-- End Old code
      }
  }
  ##
  if (J2 > 0){
      ## xi2.se = sqrt(pmax(1e-6,diag(MASS::ginv(-Hes.xi2)))) ## sqrt(diag(solve(-Hes.xi1)))
      ## xi2.low = xi2.cur - z.alpha*xi2.se ; xi2.up = xi2.cur + z.alpha*xi2.se
      lambda2.low = lambda2.mat[,"low"] ; lambda2.up = lambda2.mat[,"up"]
      ## lambda2.low = lambda2.min + exp(xi2.low) ; lambda2.up = lambda2.min + exp(xi2.up)
      ##
      ## Elements for testing absence of effect (horizontality) or linear effect for additive terms
      D1=diff(diag(K2),dif=1) ; D2=diff(diag(K2),dif=2) ## Diff matrices of order 1 & 2 for centered B-spline matrices without their last basis component for identification purposes in additve models
      ## D1=diff(diag(K2+1),dif=1)[,-(K2+1)] ; D2=diff(diag(K2+1),dif=2)[,-(K2+1)] ## Diff matrices of order 1 & 2 for centered B-spline matrices without their last basis component for identification purposes in additve models
      psi = psi2.cur ## Concerned spline coefs
      Cov = Cov.psi2 ## solve(-Hes.psi2) ## Var-Cov matrix of dispersion regr coef
      Chi2.horiz2 = Chi2.linear2 = Pval.horiz2 = Pval.linear2 = 0*ED2.cur ## Test for Horizontality and Linearity of additive terms
      ED2.low = ED2.up = 0*ED2.cur
      BWB.2 = with(obj.cur, -(Hes.psi2+P2)) ## gives  t(Xcal.2)%*%(W.psi2*Xcal.2)
      for (j in 1:J2){
          ## Horizontal ?
          ## ----------
          idx = nfixed2 + (j-1)*K2 + (1:K2)
          Chi2.horiz2[j] = sum(psi[idx] * c(solve(Cov[idx,idx]) %*% c(psi[idx]))) ## Test whether can set all spline parms to 0 ; Indeed, with a centered B-spline basis, the last spline coef is not estimated and set to zero
          Pval.horiz2[j] = 1-pchisq(Chi2.horiz2[j],obj.cur$ED2[j])
          ##
          ## Linear ?
          ## ------
          P1 = t(D1)%*%D1
          Chi2.linear2[j] = sum(c(D1%*%psi[idx]) * c(solve(D1%*%Cov[idx,idx]%*%t(D1)) %*% c(D1%*%psi[idx])))
          ## Chi2.linear2[j] = sum(c(D2%*%psi[idx]) * c(solve(D2%*%Cov[idx,idx]%*%t(D2)) %*% c(D2%*%psi[idx])))
          Pval.linear2[j] = 1-pchisq(Chi2.linear2[j],obj.cur$ED2[j]-1)
          ##
          ## Lower and upper value for effective degrees of freedom
          ## ------------------------------------------------------
          lam2.low = lam2.up = lambda2.cur
          ## Only set the jth component of lambda2 to lower or upper CI value
          lam2.low[j] = lambda2.low[j] ; lam2.up[j] = min(lambda2.up[j],1e5)
          ## !! Lower (upper) bound for lambda corresponds to upper (lower) bound for EDF !!
          ED2.up[j]= EDF.fun(obj.cur$Hes0.psi, lambda1.cur, lam2.low,
                             parms1, parms2, joint.computation=joint.computation)$ED2[j]
          ED2.low[j]= EDF.fun(obj.cur$Hes0.psi, lambda1.cur, lam2.up,
                              parms1, parms2, joint.computation=joint.computation)$ED2[j]
          ##
          ## ## Begin Old code -->
          ## P.low = Pcal.fun(nfixed2,lam2.low,Pd2.x) ; P.up = Pcal.fun(nfixed2,lam2.up,Pd2.x)
          ## temp.low = diag(solve(BWB.2 + P.low) %*% BWB.2) ## Effective dim for each regression parameter in dispersion
          ## temp.up  = diag(solve(BWB.2 + P.up ) %*% BWB.2) ##   with only jth penalty set at its lower or upper CI value
          ## idx = nfixed2 + (j-1)*K2 + (1:K2) ## Concerned diagonal elements for jth additive term
          ## ## Effective dim of the jth additive term in dispersion with lambda.j at its lower CI value
          ## ED2.up[j]  = sum(temp.low[idx])
          ## ## Effective dim of the jth additive term in dispersion with lambda.j at its upper CI value
          ## ED2.low[j] = sum(temp.up[idx])
          ## ## <-- End Old code
      }
  }
  ## Percentage of exactly observed, interval-censored and right-censored data
  ## -------------------------------------------------------------------------
  perc.obs = round(100*n.uncensored/sw,1) ## Percentage exactly observed
  perc.IC = round(100*n.IC/sw,1) ## Percentage Interval Censored
  perc.RC = round(100*n.RC/sw,1) ## Percentage Right Censored
  ##
  ## Information criteria
  ## --------------------
  AIC = obj.cur$dev + 2*obj.cur$ED.tot
  BIC = obj.cur$dev + log(sw - n.RC) * obj.cur$ED.tot
  logEvid = with(obj.cur, lpost - .5 * (ldet.fun(-Hes.psi1) + ldet.fun(-Hes.psi2)))
  ## Output
  ## ------
  ans = list(converged=converged,
             y=y, data=data,
             formula1=formula1, formula2=formula2,
             n=n,sw=sw,resp=resp,is.IC=is.IC,
             n.uncensored=n.uncensored, n.IC=n.IC, n.RC=n.RC,
             perc.obs=perc.obs,perc.IC=perc.IC,perc.RC=perc.RC,
             regr1=regr1,regr2=regr2,
             is.IC=is.IC,event=event, ## IC and event indicators
             Normality=Normality, ## Indicates whether normality requested
             derr=fit1d, phi=phi.cur, rmin=min(fit1d$knots),rmax=max(fit1d$knots), ## Error density
             K.error=K.error, knots.error=fit1d$knots,
             density.method=density.method, ## Algorithm to select the penalty parameter for the density
             fixed.loc=mat.beta, fixed.disp=mat.delta, ## Summary on the estimation of the 'fixed' effects <beta> and <delta>
             ci.level=ci.level, ## Level for credible intervals
             alpha=alpha, ## Significance level
             sandwich=sandwich, ## Indicates whether sandwich estimator for variance-covariance
             psi.method=psi.method, ## Algorithm used to estimate <psi> for given <lambda1,lambda2>
             lambda.method=lambda.method, ## Algorithm to select the penalty parameter for the additive terms
             psi1=psi1.cur, ## Location regr coefs
             bread.psi=bread.psi, Sand.psi=Sand.psi, ## Var-Cov for psi
             bread.psi1=bread.psi1, Sand.psi1=Sand.psi1, Cov.psi1=Cov.psi1, ## Var-Cov for psi1
             Vn.psi=Vn.psi,
             U.psi1=obj.cur$U.psi1, ## Gradient psi1
             psi2=psi2.cur, ## Dispersion regr coefs
             bread.psi2=bread.psi2, Sand.psi2=Sand.psi2, Cov.psi2=Cov.psi2, ## Var-Cov for psi2
             U.psi2=obj.cur$U.psi2, ## Gradient psi2
             U.psi=obj.cur$U.psi,   ## Gradient for <psi1,psi2>
             Cov.psi=obj.cur$Cov.psi,
             Hes.psi=obj.cur$Hes.psi, Hes0.psi=obj.cur$Hes0.psi,
             ## Hesb.psi=obj.cur$Hesb.psi,
             joint.computation=joint.computation)
  if (with(obj.cur, exists("v"))) ans$v = obj.cur$v
  if (J1 > 0){
    ED1 = cbind(EDF=obj.cur$ED1,Chi2=Chi2.horiz1,Pval=Pval.horiz1)
    ##    ED1 = cbind(ED.hat=obj.cur$ED1,Chi2.H0=Chi2.horiz1,Pval.H0=Pval.horiz1,Chi2.lin=Chi2.linear1,Pval.lin=Pval.linear1)
    if (is.null(addloc.lab)) rownames(ED1) = paste("f.mu.",1:J1,sep="")
    else rownames(ED1) = addloc.lab
    ans = c(ans,
            list(K1=K1, ## knots.loc=knots.1, Bgrid.loc=Bgrid.loc,
                 ## xi1=cbind(xi1.hat=xi1.cur,xi1.se=xi1.se,xi1.low=xi1.low,xi1.up=xi1.up),
                 ## U.xi1=U.xi1, U.lambda1=U.lam1,
                 ## Cov.xi1=MASS::ginv(-Hes.xi1),
                 lambda1.min=lambda1.min,
                 lambda1=lambda1.mat,
                 ## lambda1=cbind(lambda1.hat=lambda1.min+exp(xi1.cur),
                 ##               lambda1.low=lambda1.min+exp(xi1.low),lambda1.up=lambda1.min+exp(xi1.up)),
                 ED1=cbind(ED1,low=ED1.low,up=ED1.up)))
  }
  if (J2 > 0){
    ED2 = cbind(EDF=obj.cur$ED2,Chi2=Chi2.horiz2,Pval=Pval.horiz2)
##    ED2 = cbind(ED2.hat=obj.cur$ED2,Chi2.H0=Chi2.horiz2,Pval.H0=Pval.horiz2,Chi2.lin=Chi2.linear2,Pval.lin=Pval.linear2)
    if (is.null(adddisp.lab)) rownames(ED2) = paste("f.sig.",1:J2,sep="")
    else rownames(ED2) = adddisp.lab
    ##
    ans = c(ans,
            list(K2=K2, ## knots.disp=knots.2, Bgrid.disp=Bgrid.disp,
                 ## xi2=cbind(xi2.hat=xi2.cur,xi2.se=xi2.se,xi2.low=xi2.low,xi2.up=xi2.up),
                 ## U.xi2=U.xi2, U.lambda2=U.lam2,
                 ## Cov.xi2=MASS::ginv(-Hes.xi2),
                 lambda2.min=lambda2.min,
                 lambda2=lambda2.mat,
                 ## lambda2=cbind(lambda2.hat=lambda2.min+exp(xi2.cur),
                 ##               lambda2.low=lambda2.min+exp(xi2.low),lambda2.up=lambda2.min+exp(xi2.up)),
                 ED2=cbind(ED2,low=ED2.low,up=ED2.up)))
  }
  ##
  ans = c(ans,
          list(
              res=res,mu=obj.cur$mu,sd=obj.cur$sd,
              REML=REML,diag.only=diag.only,
              dev=obj.cur$dev, llik=obj.cur$llik, logEvid=logEvid,
              AIC=AIC, BIC=BIC, ED.tot=obj.cur$ED.tot,
              iter=iter,
              call=cl,
              elapsed.time=elapsed.time))
  ##
  ## Numerical computation of gradient and Hessian
  ## g.regr <- function(psi) ff(psi,lambda1.cur,lambda2.cur,REML=REML,diag.only=diag.only,
  ##                           gradient=FALSE,hessian=FALSE)$lpost
  ## Dg.regr <- function(psi) ff(psi,lambda1.cur,lambda2.cur,REML=REML,diag.only=diag.only,
  ##                            gradient=TRUE,hessian=FALSE)$U.psi

  ## HesBis = numDeriv::jacobian(func=Dg.regr, x=psi.cur)
  ## gradBis = numDeriv::grad(func=g.regr, x=psi.cur)
  ## ans = c(ans,list(gradBis=gradBis, HesBis=HesBis))
  ##
  class(ans) = "DALSM"
  ##
  if (verbose) print(ans)
  return(ans)
}
