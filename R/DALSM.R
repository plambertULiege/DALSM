#' Fit a double additive location-scale model (DALSM) with a flexible error distribution
#' @description Fit a location-scale regression model with a flexible error distribution and
#' additive terms in location (=mean) and dispersion (= log(sd)) using Laplace P-splines
#' from potentially right- and interval-censored response data.
#' @usage DALSM(y, formula1,formula2, data,
#'        K1=10, K2=10, pen.order1=2, pen.order2=2,
#'        b.tau=1e-4, lambda1.min=1, lambda2.min=1,
#'        phi.0=NULL,psi1.0=NULL,psi2.0=NULL,lambda1.0=NULL,lambda2.0=NULL,
#'        REML=FALSE, diag.only=TRUE, Normality=FALSE, sandwich=TRUE,
#'        K.error=20, rmin=NULL,rmax=NULL,
#'        ci.level=.95, iterlim=50,verbose=FALSE)
#'
#' @param y n-vector of responses or (nx2)-matrix/data.frame when interval-censoring (with the 2 elements giving the interval bounds) or right-censoring (when the element in the 2nd column equals Inf).
#' @param formula1 model formula for location (i.e. for the conditional mean).
#' @param formula2 model formula for dispersion (i.e. for the log of the conditional standard deviation).
#' @param data data frame containing the model covariates.
#' @param K1 (optional) number of B-splines to describe a given additive term in the location submodel (default: 10).
#' @param K2 (optional) number of B-splines to describe a given additive term in the dispersion submodel (default: 10).
#' @param pen.order1 (optional) penalty order for the additive terms in the location submodel (default: 2).
#' @param pen.order2 (optional) penalty order for the additive terms in the dispersion submodel (default: 2).
#' @param b.tau (optional) prior on penalty parameter \eqn{\tau} is Gamma(1,b.tau=1e-4) (for additive terms) (default: 1e-4).
#' @param lambda1.min (optional) minimal value for the penalty parameters in the additive model for location (default: 1.0).
#' @param lambda2.min (optional) minimal value for penalty parameters in the additive model for log-dispersion (default: 1.0).
#' @param phi.0 (optional) initial values for the spline parameters in the log-hazard of the standardized error distribution.
#' @param psi1.0 (optional) initial values for the location submodel parameters.
#' @param psi2.0 (optional) initial values for the dispersion submodel parameters.
#' @param lambda1.0 (optional) initial value for the J1 penalty parameters of the additive terms in the location submodel.
#' @param lambda2.0 (optional) initial value for the J2 penalty parameters of the additive terms in the dispersion submodel.
#' @param REML (optional) logical indicating if a REML correction is desired to estimate the dispersion parameters (default: FALSE).
#' @param diag.only (optional) logical indicating if only the diagonal of the Hessian needs to be corrected during REML (default: TRUE).
#' @param Normality (optional) logical indicating if Normality is assumed for the error term (default: FALSE).
#' @param sandwich (optional) logical indicating if sandwich variance estimators are needed for the regression parameters for a NP error distribution (when Normality=FALSE) ; it is forced to be TRUE when Normality=TRUE.
#' @param K.error (optional) number of B-splines to approximate the log of the error hazard on (rmin,rmax) (default: 20).
#' @param rmin (optional) minimum value for the support of the standardized error distribution (default: min(y)-sd(y)).
#' @param rmax (optional) maximum value for the support of the standardized error distribution (default: max(y)+sd(y)).
#' @param ci.level (optional) nominal level for the reported credible intervals (default: .95).
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
#' plot(fit)
#'
#' @export

DALSM = function(y, formula1,formula2, data,
                 K1=10, K2=10,pen.order1=2, pen.order2=2, b.tau=1e-4,
                 lambda1.min=1, lambda2.min=1,
                 phi.0=NULL,psi1.0=NULL,psi2.0=NULL,lambda1.0=NULL,lambda2.0=NULL,
                 REML=FALSE, diag.only=TRUE, Normality=FALSE, sandwich=TRUE,
                 K.error=20, rmin=NULL,rmax=NULL,
                 ci.level=.95,
                 iterlim=50,verbose=FALSE){
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
  event = ifelse(is.finite(yup),1,0) ## Event = 1 if non-RC, 0 otherwise
  resp = y ; resp[,2] = ifelse(is.finite(yup),yup,ylow) ## Working (nx2) response matrix with identical columns if non-IC (in particular when RC)
  ymid = apply(resp,1,mean) ## Midpoint of IC response, reported exactly observed or right-censored response otherwise
  is.IC = (ylow != yup) & (is.finite(yup)) ## Indicated whether Interval Censored (with, then, event=1)
  n.IC = sum(is.IC)
  ##
  ## Numbers of RC & IC data
  is.uncensored = (ylow==yup) & (event==1)
  n.uncensored = sum(is.uncensored)
  is.RC = (event==0)
  n.RC = sum(is.RC)
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
  Z1 = regr1$Z ; Z2 = regr2$Z
  loc.lab = colnames(regr1$Xcal) ; disp.lab = colnames(regr2$Xcal) ## Labels of the regression parms
  addloc.lab = regr1$additive.lab ; adddisp.lab = regr2$additive.lab ## Labels of the additive terms
  # addloc.lab = colnames(regr1$X) ; adddisp.lab = colnames(regr2$X) ## Labels of the additive terms
  Xcal.1 = regr1$Xcal ; Xcal.2 = regr2$Xcal
  q1 = ncol(Xcal.1) ; q2 = ncol(Xcal.2) ## Total number of regression and spline parameters
  z.alpha = qnorm(1-.5*alpha)
  ##
  ## Initial values
  ## ... for location stuff
  if (J1 > 0){
    Bcal.1 = regr1$Bcal
    if (is.null(lambda1.0)) lambda1.0 = rep(100,J1)
    names(lambda1.0) = lambda1.lab
    lambda1.cur = lambda1.0 ; P1.cur = Pcal.fun(nfixed1,lambda1.cur,Pd1.x)
    psi1.cur = c(solve(t(Xcal.1[obs,,drop=FALSE])%*%Xcal.1[obs,,drop=FALSE] + P1.cur, t(Xcal.1[obs,,drop=FALSE])%*%ymid[obs]))
  } else {
    psi1.cur = c(solve(t(Xcal.1[obs,,drop=FALSE])%*%Xcal.1[obs,,drop=FALSE], t(Xcal.1[obs,,drop=FALSE])%*%ymid[obs]))
    lambda1.cur = NULL
  }
  names(psi1.cur) = loc.lab
  ##
  ## ... for dispersion stuff
  sd0 = 1.1*sd(ymid[obs]-c(Xcal.1[obs,,drop=FALSE]%*%psi1.cur))
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
  res.fun = function(resp,mu,sd){
    res = (resp-mu)/sd
    if ((!fixed.support) & (!is.null(obj1d))){
      rmin = min(obj1d$knots) ; rmax = max(obj1d$knots)
    }
    res[res<=rmin] = .99*rmin ; res[res>=rmax] = .99*rmax ## If necessary, "Huber-like" M-estimation strategy
    return(res)
  }
  ##
  res.psi.fun = function(resp,psi){
    psi1 = psi[1:q1] ; psi2 = psi[-(1:q1)]
    mu = c(Xcal.1 %*% psi1)
    eta = c(Xcal.2 %*% psi2) ; sd = exp(eta)
    res = res.fun(resp,mu,sd)
    return(res)
  }
  ## Key function:  Computation of Log-posterior, Score, Hessian, etc.
  ## ------------
  ff = function(psi, lambda1,lambda2,
                hessian=TRUE,REML,diag.only=TRUE){
    if (REML) hessian=TRUE  ## Computation of the Hessian required for REML !!
    psi1 = psi1.cur = psi[1:q1] ; psi2 = psi2.cur = psi[-(1:q1)]
    lambda1.cur = lambda1 ; lambda2.cur = lambda2
    ## Function to interrupt computation when found necessary
    stopthis = function(){
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
    if (J1 > 0) P1.cur = Pcal.fun(nfixed1,lambda1.cur,Pd1.x) ## Update penalty matrix for location
    else P1.cur = 0
    ##
    if (J2 > 0) P2.cur = Pcal.fun(nfixed2,lambda2.cur,Pd2.x) ## Update penalty matrix for dispersion
    else P2.cur = 0
    ## lpost
    ## -----
    Dif.S = pmax(1e-6,S1.cur-S2.cur)
    lpost.cur = llik.cur = sum(ifelse(is.IC,
                                      log(Dif.S),
                                      event*(log(h1.cur)-log(sd.cur)) - H1.cur))
    if (J1 > 0) lpost.cur = lpost.cur - .5*sum(psi1.cur*c(P1.cur%*%psi1.cur))
    if (J2 > 0) lpost.cur = lpost.cur - .5*sum(psi2.cur*c(P2.cur%*%psi2.cur))
    ## Score
    ## -----
    ## <psi1>
    temp.Upsi1 = ifelse(is.IC,
                        -(1/sd.cur)*(f1.cur-f2.cur)/(Dif.S),
                        c(((event/sd.cur)*dB1)%*%phi.cur -1/sd.cur*h1.cur))
    omega.psi1 = -temp.Upsi1
    U.psi1 = c(t(Xcal.1) %*% omega.psi1) ## Score.psi1
    if (J1 > 0) U.psi1 = U.psi1 - c(P1.cur%*%psi1.cur) ## Score.psi1 if additive component(s)
    ##
    ## <psi2>
    temp.Upsi2 = ifelse(is.IC,
                        -(res[,1]*f1.cur-res[,2]*f2.cur)/(Dif.S),
                        c((event*res[,1]*dB1)%*%phi.cur + event-res[,1]*h1.cur))
    omega.psi2 = -temp.Upsi2
    U.psi2 = c(t(Xcal.2) %*% omega.psi2) ## Score.psi2
    if (J2 > 0) U.psi2 = U.psi2 - c(P2.cur%*%psi2.cur) ## Score.psi2 if additive component(s)
    ##
    U.psi = c(U.psi1,U.psi2) ## Score for <psi> = <psi1,psi2>
    if (any(is.na(U.psi))) return(stopthis())
    ##
    ## Computation of Hessian
    ## ----------------------
    Mcal.1 = Mcal.2 = NULL ## Quantities required to update <lambda1> & <lambda2>
    if (hessian){
      ## <psi1>
      g1.cur = c(dB1%*%phi.cur) - h1.cur ; g2.cur = c(dB2%*%phi.cur) - h2.cur
      W.psi1 = (1/sd.cur^2)*
        ifelse(is.IC,
               (f1.cur*g1.cur-f2.cur*g2.cur)/(Dif.S) + (f1.cur-f2.cur)^2/(Dif.S)^2,
               -c((event*d2B-h1.cur*dB1)%*%phi.cur)) ## Weight.psi1
      Hes.psi1 = -t(Xcal.1)%*%(W.psi1*Xcal.1) ## Hessian.psi1
      if (J1 > 0) {
        Hes.psi1 = Hes.psi1 - P1.cur ## Hessian.psi1 if additive component(s)
        tZWB.1 = t(Z1)%*%(W.psi1*Bcal.1) ; tZWZ.1 = t(Z1)%*%(W.psi1*Z1)
        inv.tZWZ.1 = try(solve(tZWZ.1),TRUE)
        if (!is.matrix(inv.tZWZ.1)) return(stopthis()) ## stopthis() ## PHL
        Mcal.1 = t(Bcal.1)%*%(W.psi1*Bcal.1) -t(tZWB.1)%*%inv.tZWZ.1%*%tZWB.1 ## Required to update <lambda1>
      }
      ## <psi2>
      m1.cur = 1 + res[,1]*(c(dB1%*%phi.cur)-h1.cur) ; m2.cur = 1 + res[,2]*(c(dB2%*%phi.cur)-h2.cur)
      W.psi2 = ifelse(is.IC,
                      (res[,1]*f1.cur*m1.cur-res[,2]*f2.cur*m2.cur)/(Dif.S) +
                        (res[,1]*f1.cur-res[,2]*f2.cur)^2/(Dif.S)^2,
                      c(-event*(res[,1]^2)*c(d2B%*%phi.cur)
                        -res[,1]*(event-h1.cur*res[,1])*c(dB1%*%phi.cur)
                        +h1.cur*res[,1])) ## Weight.psi2
      Hes.psi2 = -t(Xcal.2)%*%(W.psi2*Xcal.2) ## Hessian.psi2
      if (J2 > 0){
        Hes.psi2 = Hes.psi2 - P2.cur ## Hessian.psi2 if additive component(s)
        tZWB.2 = t(Z2)%*%(W.psi2*Bcal.2) ; tZWZ.2 = t(Z2)%*%(W.psi2*Z2)
        inv.tZWZ.2 = try(solve(tZWZ.2),TRUE)
        if (!is.matrix(inv.tZWZ.2)) return(stopthis()) ## stopthis() ## PHL
        Mcal.2 = t(Bcal.2)%*%(W.psi2*Bcal.2) -t(tZWB.2)%*%inv.tZWZ.2%*%tZWB.2 ## Required to update <lambda2>
      }
      ## Cross-derivatives
      W.cross = (1/sd.cur)*
        ifelse(is.IC,
               (f1.cur-f2.cur)/(Dif.S) +
                 (res[,1]*f1.cur*g1.cur-res[,2]*f2.cur*g2.cur)/(Dif.S) +
                 (f1.cur-f2.cur)*(res[,1]*f1.cur-res[,2]*f2.cur)/(Dif.S)^2,
               -event*res[,1]*c(d2B%*%phi.cur) -
                 (event-h1.cur*res[,1])*c(dB1%*%phi.cur) +
                 h1.cur) ## Cross-Weight.psi1.psi2
      Hes.21 = -t(Xcal.2)%*%(W.cross*Xcal.1) ## Hessian for cross derivatives
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
      ldet.cur = logDet.fun(-Hes.psi1) ## logdet(X1'W1 X1 + K1
      lpost.cur = lpost.cur -.5*ldet.cur  ## REML: add -.5*logdet(X1'W1 X1 + K1)
    }
    ##
    ## Effective dimensions for additive terms
    ## ---------------------------------------
    ED1 = rep(0,J1) ; ED2 = rep(0,J2)
    if ((J1 > 0)&(is.matrix(Cov.psi))){
      idx1 = 1:q1
      temp.1 = diag(Cov.psi[idx1,idx1]%*%(-Hes.psi1-P1.cur)) ## Effective dim for each regression parameter in location
      for (j in 1:J1){
        idx = nfixed1 + (j-1)*K1 + (1:K1)
        ED1[j] = sum(temp.1[idx]) ## Effective dim of the jth additive term in location
      }
      names(ED1) = colnames(regr1$X)
    }
    if ((J2 > 0)&(is.matrix(Cov.psi))){
      idx2 = (q1+1):(q1+q2)
      temp.2 = diag(Cov.psi[idx2,idx2]%*%(-Hes.psi2-P2.cur)) ## Effective dim for each regression parameter in location
      for (j in 1:J2){
        idx = nfixed2 + (j-1)*K2 + (1:K2)
        ED2[j] = sum(temp.2[idx]) ## Effective dim of the jth additive term in location
      }
      names(ED2) = colnames(regr2$X)
    }
    ## Output
    ## ------
    if (any(is.na(U.psi))) return(stopthis())
    ans = list(y=y, data=data,
               formula1=formula1, formula2=formula2,
               psi1=psi1,psi2=psi2, lambda1=lambda1,lambda2=lambda2,
               psi=c(psi1,psi2),
               res=res, mu=mu.cur, sd=sd.cur,
               lpost=lpost.cur, llik=llik.cur,
               U.psi1=U.psi1, U.psi2=U.psi2, U.psi=U.psi,
               omega.psi1=omega.psi1, omega.psi2=omega.psi2,
               h1=h1.cur, h2=h2.cur, H1=H1.cur, H2=H2.cur,
               f1=f1.cur, f2=f2.cur, S1=S1.cur, S2=S2.cur,
               ldet=ldet.cur,
               P1=P1.cur, P2=P2.cur,
               ED1=ED1, ED2=ED2,
               Mcal.1=Mcal.1, Mcal.2=Mcal.2,
               REML=REML,hessian=hessian)
    if (hessian){
      ans = c(ans,
              list(Hes.psi1=Hes.psi1, Hes.psi2=Hes.psi2,
                   Hes.psi=Hes.psi,
                   Cov.psi=Cov.psi,
                   dHes.REML=dHes.REML, ddiagHes.REML=ddiagHes.REML)
      )
    }
    ##
    return(ans)
  }
  ##
  ## ESTIMATION
  ## ----------
  converged = FALSE
  eps.Score = 1e-1
  iter = 0
  ptm <- proc.time()
  ## ///////////////////////////////////////
  ## ----> BEGIN of the Estimation algorithm
  ## ///////////////////////////////////////
  ##
  KK = K.error ## Number of spline parameters for the (log-hazard of the) stdzd error density
  ##
  ## Spline parameters for (the log-hazard of) a N(0,1) density
  phi.Norm = function(knots){
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
    psi1.old = psi1.cur ; psi2.old = psi2.cur ; psi.old = c(psi1.old,psi2.old)
    ## Evaluate fitted location and dispersion
    obj1d = NULL ; res = res.psi.fun(resp,c(psi1.cur,psi2.cur))
    ##
    ## -1- Update B-spline estimation of the error density
    ## *****************************************************
    if (!final.iteration){
      KK = K.error
      ## Estimate error density
      obj1d = Dens1d(res, event=event,
                     ymin=rmin, ymax=rmax, K=K.error,
                     equid.knots=TRUE)
      ## Force a N(0,1) density during a couple of iterations OR till the end if Normality=TRUE
      fixed.phi = ifelse((iter <= 2) | Normality ,TRUE,FALSE)
      if (iter <=2){
        phi.cur = phi.ref = phi.Norm(obj1d$knots)
      }
      ##
      fit1d.0 = densityLPS(obj1d,
                       phi0=phi.cur,
                       fixed.phi=fixed.phi,phi.ref=phi.ref,
                       is.density=TRUE,Mean0=NULL,Var0=NULL,
                       method="LPS",verbose=FALSE)
      fit1d = densityLPS(obj1d,
                     phi0=fit1d.0$phi,
                     fixed.phi=fixed.phi, phi.ref=phi.ref,
                     is.density=TRUE,Mean0=0,Var0=1,
                     method="LPS",verbose=FALSE)
      fit1d$n=n ; fit1d$n.uncensored=n.uncensored ; fit1d$n.IC=n.IC ; fit1d$n.RC=n.RC
      phi.cur = fit1d$phi ## Current spline parameter estimates in error density estimation
    }
    ## -2- Update <psi1> & <psi2> (i.e. regression coefs related to location)
    ## **************************
    ## Computation of the Newton-Raphson step
    dpsi.fun = function(obj,hessian=TRUE){
      if (hessian){
        ## Newton-Raphson step using classical minus inverse-Hessian for Var-Cov
        dpsi = with(obj, Cov.psi%*%U.psi) ## Increment.psi
      } else {
        ## Newton-Raphson step using Sandwich Var-Cov
        temp = cbind(Xcal.1 * obj$omega.psi1, Xcal.2 * obj$omega.psi2)
        Vn.psi = crossprod(temp) ## PHL
        idx1 = 1:q1 ; idx2 = (q1+1):(q1+q2)
        if (J1 > 0) Vn.psi[idx1,idx1] = Vn.psi[idx1,idx1] + obj$P1
        if (J2 > 0) Vn.psi[idx2,idx2] = Vn.psi[idx2,idx2] + obj$P2
        bread.psi = obj$Cov.psi ## Minus-Inverse Hessian for <psi>
        Sand.psi = bread.psi %*% Vn.psi %*% bread.psi ## Sandwich estimator
        dpsi = Sand.psi %*% obj$U.psi
      }
      return(dpsi)
    }
    ##
    ok.psi1 = ok.psi2 = FALSE ## Monitor convergence in estimation of <psi1> & <psi2>
    ok.psi = (ok.psi1 & ok.psi2)
    iter.psi = 0
    while(!ok.psi){
      iter.psi = iter.psi + 1
      obj.cur = ff(c(psi1.cur,psi2.cur),lambda1.cur,lambda2.cur,REML=REML,diag.only=diag.only) ## Current lpost, Score, Hessian...
      ED1.old = obj.cur$ED1 ; ED2.old = obj.cur$ED2
      ok.psi = (max(abs(obj.cur$U.psi)) < eps.Score) | (iter.psi > 10) ## Check if can already stop (or if algorithm stuck)
      if (is.na(ok.psi)) ok.psi = FALSE
      if (ok.psi){
        ok.psi1 = (max(abs(obj.cur$U.psi1)) < eps.Score)
        ok.psi2 = (max(abs(obj.cur$U.psi2)) < eps.Score)
        break
      }
      ## Newton-Raphson step using classical minus inverse-Hessian for Var-Cov
      dpsi = dpsi.fun(obj.cur,hessian=TRUE) ## Increment.psi
      accept = FALSE ; step = 1
      while(!accept){ ## Make proposals till acceptance
        psi.cur = obj.cur$psi ## Current values for the regression parameters
        nrep = 0
        repeat { ## Repeat step-halving directly till sensible proposal for the std error
          nrep = nrep + 1
          if (nrep == 3) dpsi = dpsi.fun(obj.cur,hessian=FALSE) ## Recompute Increment.psi
          psi.prop = psi.cur + step*dpsi ## Update.psi
          psi1.prop = psi.prop[1:q1] ; psi2.prop = psi.prop[-(1:q1)] ## Separate into <psi1.prop> & <psi2.prop>
          eta.prop = c(Xcal.2 %*% psi2.prop) ## log(sd.prop)
          if (max(eta.prop) < 40 | (nrep > 10)) break ## Make sure that does not generate irrealistic std error
          step = .5*step
        }
        if (verbose & (nrep > 1)) cat("nrep = ",nrep,"\n")
        ## Evaluate lpost, Score, Hessian, etc. at <psi.prop>
        obj.prop = ff(c(psi1.prop,psi2.prop), lambda1.cur,lambda2.cur,REML=REML,diag.only=diag.only)
        ## Check if increases <lpost>
        if (is.na(obj.prop$lpost)) accept = FALSE ## Reject if <lpost> is NA
        else {
          if (!is.finite(obj.prop$lpost)) accept = FALSE ## Reject if <lpost> is Inf
          else accept = ((obj.prop$lpost-obj.cur$lpost)>= -1) ## Accept if increases <lpost>
        }
        ## ... otherwise, step-halving
        if (!accept) step = .5*step ## Step-halving if proposal rejected
        if (step < 1e-3){
          break
        }
      }
      if ((step !=1)&verbose) cat("Step-halving for regression parameters: iter.psi=",iter.psi," ; step=",step,"\n")
      ##
      ## Record the accepted update for <psi1> and <res>
      ## -----------------------------------------------
      if (accept){
        obj.cur = obj.prop
        psi.cur = obj.cur$psi
        psi1.cur = obj.cur$psi1 ; psi2.cur = obj.cur$psi2
      }
      ## Stopping criteria based on the Score
      ## ------------------------------------
      ok.psi1 = (max(abs(obj.cur$U.psi1)) < eps.Score) ; if (is.na(ok.psi1)) ok.psi1 = FALSE
      ok.psi2 = (max(abs(obj.cur$U.psi2)) < eps.Score) ; if (is.na(ok.psi2)) ok.psi2 = FALSE
      ok.psi = (ok.psi1 & ok.psi2)
      if (iter.psi > 10) break ## Give up trying to update if that doesn't work !!
      if (!is.logical(ok.psi)) ok.psi = FALSE
    }
    ##
    ## -3- Update <lambda1> & <lambda2> (i.e. penalty vectors for location and dispersion additive terms)
    ## ********************************
    res = obj.cur$res
    itermin = 3 ## Only start updating after <itermin> iterations
    update.lambda = ifelse(iter <= itermin, FALSE, (J1 > 0)|(J2 > 0))
    iter.lam1 = iter.lam2 = 0
    if (update.lambda & !final.iteration){
      ## -4a- LOCATION
      ## -------------
      ok.xi1 = ifelse(J1 > 0,FALSE,TRUE) ; iter.lam1 = 0
      ##
      ## Home-made Newton-Raphson
      Mcal.1 = obj.cur$Mcal.1
      while(!ok.xi1){
        iter.lam1 = iter.lam1 + 1
        P1.cur = Pcal.fun(nfixed1,lambda1.cur,Pd1.x) ## Update penalty matrix
        iMpen.1 = solve(Mcal.1 + P1.cur[-(1:nfixed1),-(1:nfixed1)])
        ## Score.lambda1
        U.lam1 = U.xi1 = rep(0,J1)
        Rj.1 = list()
        ##
        for (j in 1:J1){
          lam1.j = rep(0,J1) ; lam1.j[j] = lambda1.cur[j]
          Pj.1 = Pcal.fun(nfixed1,lam1.j,Pd1.x)[-(1:nfixed1),-(1:nfixed1)] ## Update penalty matrix
          Rj.1[[j]] = iMpen.1%*%(Pj.1/lam1.j[j])
          idx = nfixed1 + (j-1)*K1 + (1:K1)
          theta.j = psi1.cur[idx]
          quad1.cur[j] = sum(theta.j*c(Pd1.x%*%theta.j)) + b.tau  ## <-----
          U.lam1[j] = .5*(K1-pen.order.1)/lam1.j[j] -.5*quad1.cur[j] -.5*sum(diag(Rj.1[[j]]))
        }
        ## Score.xi1  where  lambda1 = lambda1.min + exp(xi1)
        U.xi1 = U.lam1 * (lambda1.cur - lambda1.min)
        ##
        ## Hessian.lambda1
        Hes.lam1 = diag(0,J1)
        Only.diagHes1 = FALSE ## PHL
        if (!Only.diagHes1){
          for (j in 1:J1){
            for (k in j:J1){
              temp = 0
              if (j==k) temp = temp -.5*(K1-pen.order.1)/lambda1.cur[j]^2
              temp = temp + .5*sum(t(Rj.1[[j]])*Rj.1[[k]]) ## tr(XY) = sum(t(X)*Y)
              Hes.lam1[j,k] = Hes.lam1[k,j] = temp
            }
          }
        } else {
          for (j in 1:J1) Hes.lam1[j,j] = -.5*(K1-pen.order.1)/lambda1.cur[j]^2 +.5*sum(t(Rj.1[[j]])*Rj.1[[j]])
        }
        ## Hessian.xi1  where  lambda1 = lambda1.min + exp(xi1)
        Hes.xi1 = t((lambda1.cur-lambda1.min) * t((lambda1.cur-lambda1.min) * Hes.lam1))
        diag(Hes.xi1) = diag(Hes.xi1) + U.lam1 * (lambda1.cur-lambda1.min)
        ## Update xi1  where  lambda1 = lambda1.min + exp(xi1)
        xi1.cur = log(lambda1.cur-lambda1.min) ; names(xi1.cur) = addloc.lab ## paste("f.mu.",1:J1,sep="")
        dxi1 = -c(MASS::ginv(Hes.xi1) %*% U.xi1) ## Increment.xi1
        ##
        step = 1
        accept = FALSE
        while(!accept){
          xi1.prop = xi1.cur + step*dxi1
          accept = all(xi1.prop < 10) ## all(is.finite(exp(xi1.prop)))
          if (!accept) step = .5*step
        }
        if ((step !=1)&verbose) cat("Step-halving for lambda1: step=",step,"\n")
        xi1.cur = xi1.prop
        lambda1.cur = lambda1.min + exp(xi1.cur)
        ##
        mat = Mcal.1 + P1.cur[-(1:nfixed1),-(1:nfixed1)]
        ## Convergence ?
        ok.xi1 = (max(abs(U.xi1)) < eps.Score) | (iter.lam1 > 20)
      }
      ## -4b- DISPERSION
      ## ---------------
      ok.xi2 = ifelse(J2 > 0,FALSE,TRUE) ; iter.lam2 = 0
      ##
      ## Home-made Newton-Raphson
      Mcal.2 = obj.cur$Mcal.2
      while(!ok.xi2){
        iter.lam2 = iter.lam2 + 1
        P2.cur = Pcal.fun(nfixed2,lambda2.cur,Pd2.x) ## Update penalty matrix
        iMpen.2 = MASS::ginv(Mcal.2 + P2.cur[-(1:nfixed2),-(1:nfixed2)])
        ## Score.lambda2
        U.lam2 = U.xi2 = rep(0,J2)
        Rj.2 = list()
        for (j in 1:J2){
          lam2.j = rep(0,J2) ; lam2.j[j] = lambda2.cur[j]
          Pj.2 = Pcal.fun(nfixed2,lam2.j,Pd2.x)[-(1:nfixed2),-(1:nfixed2)] ## Update penalty matrix
          Rj.2[[j]] = iMpen.2%*%(Pj.2/lam2.j[j])
          idx = nfixed2 + (j-1)*K2 + (1:K2)
          theta2.j = psi2.cur[idx]
          quad2.cur[j] = sum(theta2.j*c(Pd2.x%*%theta2.j)) + b.tau  ## <-----
          U.lam2[j] = .5*(K2-pen.order.2)/lam2.j[j] -.5*quad2.cur[j] -.5*sum(diag(Rj.2[[j]]))
        }
        ## Score.xi2  where  lambda2 = lambda2.min + exp(xi2)
        U.xi2 = U.lam2 * (lambda2.cur-lambda2.min)
        ##
        ## Hessian.lambda2
        Hes.lam2 = diag(0,J2)
        Only.diagHes2 = FALSE
        if (!Only.diagHes2){
          for (j in 1:J2){
            for (k in j:J2){
              temp = 0
              if (j==k) temp = temp -.5*(K2-pen.order.2)/lambda2.cur[j]^2
              temp = temp +.5*sum(t(Rj.2[[j]])*Rj.2[[k]]) ## tr(XY) = sum(t(X)*Y)
              Hes.lam2[j,k] = Hes.lam2[k,j] = temp
            }
          }
        } else {
          for (j in 1:J2) Hes.lam2[j,j] = -.5*(K2-pen.order.2)/lambda2.cur[j]^2 +.5*sum(t(Rj.2[[j]])*Rj.2[[j]])
        }
        ## Hessian.xi2  where  lambda2 = lambda2.min + exp(xi2)
        Hes.xi2 = t((lambda2.cur-lambda2.min) * t((lambda2.cur-lambda2.min) * Hes.lam2))
        diag(Hes.xi2) = diag(Hes.xi2) + U.lam2 * (lambda2.cur-lambda2.min)
        ## Update xi2  where  lambda2 = lambda2.min + exp(xi2)
        xi2.cur = log(lambda2.cur-lambda2.min) ; names(xi2.cur) = adddisp.lab ##paste("f.sig.",1:J2,sep="")
        dxi2 = -c(MASS::ginv(Hes.xi2) %*% U.xi2) ## Increment.xi2
        xi2.cur = xi2.cur + dxi2
        lambda2.cur = lambda2.min + exp(xi2.cur)
        ## Convergence ?
        ok.xi2 = (max(abs(U.xi2)) < eps.Score) | (iter.lam2 > 20)
      }
    }
    ##
    ## Global Stopping rule
    ## ********************
    L2 = function(x) sqrt(sum(x^2))
    L1 = function(x) sum(abs(x))
    ##
    eps = eps.Score ## 1e-3
    ##
    if (!final.iteration){
      converged.loc = (max(abs(obj.cur$U.psi1)) < eps)   ## L2 criterion
      converged.disp = (max(abs(obj.cur$U.psi2)) < eps)  ## L2 criterion
      ## Check if the EDs of the additive terms stabilize
      if (J1 > 0) converged.ED1 = ifelse(iter<=itermin, FALSE, ((L1(obj.cur$ED1-ED1.old) / L1(ED1.old)) < eps))
      else converged.ED1 = TRUE
      if (J2 > 0) converged.ED2 = ifelse(iter<=itermin, FALSE, ((L1(obj.cur$ED2-ED2.old) / L1(ED2.old)) < eps))
      else converged.ED2 = TRUE
    }
    converged.detail = c(ok.psi1,ok.psi2,converged.ED1,converged.ED2)
    ok.psi = (ok.psi1 & ok.psi2)
    ## When the EDs of the additive terms stabilize,
    ##  perform a final iteration to get the final estimates for the regression parms <psi1> & <psi2>
    green = all(c(ok.psi,converged.ED1,converged.ED2)) ## Necessary condition to complete the estimation process
    ## Necessary & sufficient condition to complete the estimation process: green light followed by one last iteration
    if (final.iteration & green) converged = TRUE
    final.iteration = green ## One final iteration required after the 1st green light to complete the estimation process
    ##
    if (verbose){
      cat("\niter:",iter," (",iter.psi,iter.lam1,iter.lam2,") ; convergence:",converged.detail,"\n")
      cat("Score for <psi1> in",range(obj.cur$U.psi1),"\n")
      cat("Score for <psi2> in",range(obj.cur$U.psi2),"\n")
      if (J1 > 0){
        cat("lambda1.cur:",lambda1.cur,"\n")
        cat("ED1 before & after: ",ED1.old," ---> ",obj.cur$ED1,"\n")
      }
      if (J2 > 0){
        cat("lambda2.cur:",lambda2.cur,"\n")
        cat("ED2 before & after: ",ED2.old," ---> ",obj.cur$ED2,"\n")
      }
    }
    ##
    if (iter > iterlim){
      cat("Algorithm did not converge after",iter,"iterations\n")
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
  temp = cbind(Xcal.1 * obj.cur$omega.psi1, Xcal.2 * obj.cur$omega.psi2)
  Vn.psi = crossprod(temp)
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
  Sand.psi = bread.psi %*% Vn.psi %*% bread.psi ## Sandwich estimator
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
    xi1.se = sqrt(pmax(1e-6,diag(MASS::ginv(-Hes.xi1)))) ## sqrt(diag(solve(-Hes.xi1)))
    xi1.low = xi1.cur - z.alpha*xi1.se ; xi1.up = xi1.cur + z.alpha*xi1.se
    lambda1.low = lambda1.min + exp(xi1.low) ; lambda1.up = lambda1.min + exp(xi1.up)
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
      lam1.low = lam1.up = lambda1.cur
      lam1.low[j] = lambda1.low[j] ; lam1.up[j] = min(lambda1.up[j],1e5) ## Only set the jth component of lambda1 to lower or upper CI value
      P.low = Pcal.fun(nfixed1,lam1.low,Pd1.x) ; P.up = Pcal.fun(nfixed1,lam1.up,Pd1.x)
      temp.low = diag(solve(BWB.1 + P.low) %*% BWB.1) ## Effective dim for each regression parameter in location
      temp.up  = diag(solve(BWB.1 + P.up ) %*% BWB.1) ##   with only jth penalty set at its lower or upper CI value
      idx = nfixed1 + (j-1)*K1 + (1:K1) ## Concerned diagonal elements for jth additive term
      ED1.up[j]  = sum(temp.low[idx]) ## Effective dim of the jth additive term in location with lambda.j at its lower CI value
      ED1.low[j] = sum(temp.up[idx])  ## Effective dim of the jth additive term in location with lambda.j at its upper CI value
    }
  }
  ##
  if (J2 > 0){
    xi2.se = sqrt(pmax(1e-6,diag(MASS::ginv(-Hes.xi2)))) ## sqrt(diag(solve(-Hes.xi1)))
    xi2.low = xi2.cur - z.alpha*xi2.se ; xi2.up = xi2.cur + z.alpha*xi2.se
    lambda2.low = lambda2.min + exp(xi2.low) ; lambda2.up = lambda2.min + exp(xi2.up)
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
      ## Lower and upper value for effective degrees of freedom
      lam2.low = lam2.up = lambda2.cur
      lam2.low[j] = lambda2.low[j] ; lam2.up[j] = min(lambda2.up[j],1e5) ## Only set the jth component of lambda2 to lower or upper CI value
      P.low = Pcal.fun(nfixed2,lam2.low,Pd2.x) ; P.up = Pcal.fun(nfixed2,lam2.up,Pd2.x)
      temp.low = diag(solve(BWB.2 + P.low) %*% BWB.2) ## Effective dim for each regression parameter in dispersion
      temp.up  = diag(solve(BWB.2 + P.up ) %*% BWB.2) ##   with only jth penalty set at its lower or upper CI value
      idx = nfixed2 + (j-1)*K2 + (1:K2) ## Concerned diagonal elements for jth additive term
      ED2.up[j]  = sum(temp.low[idx]) ## Effective dim of the jth additive term in dispersion with lambda.j at its lower CI value
      ED2.low[j] = sum(temp.up[idx])  ## Effective dim of the jth additive term in dispersion with lambda.j at its upper CI value
    }
  }
  # ## Expected Stdzd residuals
  # expctd.res = res
  # ugrid = fit1d$ugrid ; du = fit1d$du ; dens.grid = fit1d$dens.grid ; bins = fit1d$bins
  # ## Conditional expectation of the error given that larger than bins[1:(length(bins)-1)]
  # temp = rev( cumsum(rev(ugrid*dens.grid*du)) / cumsum(rev(dens.grid*du)) )
  # if (sum(1-event)>0) expctd.res[event==0] = temp[(res[event==0]-bins[1])%/%du+1]
  ##
  ## Percentage of exactly observed, interval-censored and right-censored data
  perc.obs = round(100*n.uncensored/n,1) ## Percentage exactly observed
  perc.IC = round(100*n.IC/n,1) ## Percentage Interval Censored
  perc.RC = round(100*n.RC/n,1) ## Percentage Right Censored
  ##
  ## Output
  ans = list(converged=converged,
             y=y, data=data,
             formula1=formula1, formula2=formula2,
             n=n,resp=resp,is.IC=is.IC,
             n.uncensored=n.uncensored, n.IC=n.IC, n.RC=n.RC,
             perc.obs=perc.obs,perc.IC=perc.IC,perc.RC=perc.RC,
             regr1=regr1,regr2=regr2,
             is.IC=is.IC,event=event, ## IC and event indicators
             Normality=Normality, ## Indicates whether normality requested
             derr=fit1d, phi=phi.cur, rmin=min(fit1d$knots),rmax=max(fit1d$knots), K.error=K.error, knots.error=fit1d$knots, ## Error density
             fixed.loc=mat.beta, fixed.disp=mat.delta, ## Summary on the estimation of the 'fixed' effects <beta> and <delta>
             ci.level=ci.level, ## Level for credible intervals
             alpha=alpha, ## Significance level
             sandwich=sandwich, ## Indicates whether sandwich estimator for variance-covariance
             psi1=psi1.cur, ## Location regr coefs
             bread.psi1=bread.psi1, Sand.psi1=Sand.psi1, Cov.psi1=Cov.psi1, ## Var-Cov for psi1
             U.psi1=obj.cur$U.psi1, ## Gradient psi1
             psi2=psi2.cur, ## Dispersion regr coefs
             bread.psi2=bread.psi2, Sand.psi2=Sand.psi2, Cov.psi2=Cov.psi2, ## Var-Cov for psi2
             U.psi2=obj.cur$U.psi2, ## Gradient psi2
             U.psi=obj.cur$U.psi,   ## Gradient for <psi1,psi2>
             Cov.psi=obj.cur$Cov.psi)
  if (J1 > 0){
    ED1 = cbind(ED.hat=obj.cur$ED1,Chi2=Chi2.horiz1,Pval=Pval.horiz1)
##    ED1 = cbind(ED.hat=obj.cur$ED1,Chi2.H0=Chi2.horiz1,Pval.H0=Pval.horiz1,Chi2.lin=Chi2.linear1,Pval.lin=Pval.linear1)
    if (is.null(addloc.lab)) rownames(ED1) = paste("f.mu.",1:J1,sep="")
    else rownames(ED1) = addloc.lab
    ans = c(ans,
            list(K1=K1, ## knots.loc=knots.1, Bgrid.loc=Bgrid.loc,
                 xi1=cbind(xi1.hat=xi1.cur,xi1.se=xi1.se,xi1.low=xi1.low,xi1.up=xi1.up),
                 U.xi1=U.xi1, U.lambda1=U.lam1,
                 Cov.xi1=MASS::ginv(-Hes.xi1),
                 lambda1.min=lambda1.min,
                 lambda1=cbind(lambda1.hat=lambda1.min+exp(xi1.cur),
                               lambda1.low=lambda1.min+exp(xi1.low),lambda1.up=lambda1.min+exp(xi1.up)),
                 ED1=cbind(ED1,low=ED1.low,up=ED1.up)))
  }
  if (J2 > 0){
    ED2 = cbind(ED.hat=obj.cur$ED2,Chi2=Chi2.horiz2,Pval=Pval.horiz2)
##    ED2 = cbind(ED2.hat=obj.cur$ED2,Chi2.H0=Chi2.horiz2,Pval.H0=Pval.horiz2,Chi2.lin=Chi2.linear2,Pval.lin=Pval.linear2)
    if (is.null(adddisp.lab)) rownames(ED2) = paste("f.sig.",1:J2,sep="")
    else rownames(ED2) = adddisp.lab
    ##
    ans = c(ans,
            list(K2=K2, ## knots.disp=knots.2, Bgrid.disp=Bgrid.disp,
                 xi2=cbind(xi2.hat=xi2.cur,xi2.se=xi2.se,xi2.low=xi2.low,xi2.up=xi2.up),
                 U.xi2=U.xi2, U.lambda2=U.lam2,
                 Cov.xi2=MASS::ginv(-Hes.xi2),
                 lambda2.min=lambda2.min,
                 lambda2=cbind(lambda2.hat=lambda2.min+exp(xi2.cur),
                               lambda2.low=lambda2.min+exp(xi2.low),lambda2.up=lambda2.min+exp(xi2.up)),
                 ED2=cbind(ED2,low=ED2.low,up=ED2.up)))
  }
  ##
  ans = c(ans,
          list(
               res=res,mu=obj.cur$mu,sd=obj.cur$sd,
               # expctd.res=expctd.res,
               REML=REML,diag.only=diag.only,
               iter=iter,elapsed.time=elapsed.time))
  ##
  class(ans) = "DALSM"
  ##
  if (verbose) print(ans)
  return(ans)
}
