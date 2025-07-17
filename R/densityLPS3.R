## Author: Philippe LAMBERT (ULg, UCL, Belgium), Nov 2018 ;
##  Updated in October 2024 to handle weights.
##  Updated in May 2025: code cleaning with Newton-Raphson
###################################################################################
#' Constrained density estimation from censored data using Laplace P-splines
#' @description P-spline estimation of the density (pdf), cumulative distribution (cdf),
#' hazard and cumulative hazard functions from interval- or right-censored data under possible marginal
#' mean and/or variance constraints. The penalty parameter \eqn{\tau} tuning the smoothness of
#' the log-hazard can be selected using the Laplace P-splines (LPS) method maximizing an approximation to the marginal posterior
#' of \eqn{\tau} (also named the 'evidence')  or using Schall's method.
#' @usage densityLPS2(obj.data,
#'        is.density=TRUE, Mean0=NULL, Var0=NULL,
#'        fixed.penalty=FALSE, method=c("LPS","Schall"),
#'        fixed.phi=FALSE,phi.ref=NULL, phi0=NULL,tau0=exp(5),tau.min=.1,
#'        itermax=21, verbose=FALSE)
#' @param obj.data a list created from potentially right- or interval-censored data using \code{\link{Dens1d}}. It includes summary statistics, the assumed density support, the knots for the B-spline basis, etc.
#' @param is.density (optional) logical indicating whether the estimated density should integrate to 1.0 over the range of the knots in obj.data$knots (default: TRUE).
#' @param Mean0 (optional) constrained value for the mean of the fitted density (defaut: NULL).
#' @param Var0 (optional) constrained value for the variance of the fitted density (defaut: NULL).
#' @param fixed.penalty (optional) logical indicating whether the penalty parameter should be selected from the data (\code{fixed.penalty}=FALSE) or fixed (\code{fixed.penalty}=TRUE) at its initial value \eqn{\tau_0}.
#' @param method method used for the selection of the penalty parameter: "LPS" (by maximizing the marginal posterior for \eqn{\tau}, cf. "Laplace P-splines") or "Schall" (Schall's method).
#' @param fixed.phi (optional) logical indicating whether the spline parameters are fixed (\code{fixed.phi}=TRUE) or estimated from the data (default: \code{fixed.phi}=FALSE).
#' @param phi.ref (optional) reference value for the spline parameters with respect to which deviations are penalized (default: zero vector).
#' @param phi0 starting value for the spline parameters (default: spline parameters corresponding to a Student density with 5 DF).
#' @param tau0 (optional) initial value for the penalty parameter \eqn{\tau} (default: exp(5)).
#' @param tau.min (optional) minimal value for the penalty parameter \eqn{\tau} (default: .1).
#' @param Lnorm Default Lp-norm for vectors: "Linf" (default) or "L2".
#' @param grad.tol Tolerance value for the Lp-norm of the gradient at convergence (default: 1e-2).
#' @param RDM.tol Tolerance value for the RDM (= Relative Damping Measure) at convergence (default: 1e-2).
#' @param fun.tol Tolerance value for changes in function \code{g} to declare convergence (default: 1e-3).
#' @param itermax Maximum number of iterations (after which the algorithm is interrupted with a non-convergence diagnostic) (default: 101).
#' @param verbose (optional) logical indicating whether estimation step details should be displayed (default: FALSE).
#'
#' @return a \code{\link{densLPS.object}} containing the density estimation results.
#' @export

#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references Lambert, P. (2021). Fast Bayesian inference using Laplace approximations
#' in nonparametric double additive location-scale models with right- and
#' interval-censored data.
#' \emph{Computational Statistics and Data Analysis}, 161: 107250.
#' <doi:10.1016/j.csda.2021.107250>
#'
#' @seealso \code{\link{densLPS.object}}, \code{\link{print.densLPS}}, \code{\link{plot.densLPS}}, \code{\link{Dens1d.object}}, \code{\link{Dens1d}}.
#'
#' @examples
#' library(DALSM)
#'
#' ## Example 1: estimation of the error density in a DALSM model
#' require(DALSM)
#' data(DALSM_IncomeData)
#' resp = DALSM_IncomeData[,1:2]
#' fit = DALSM(y=resp,
#'             formula1 = ~twoincomes+s(age)+s(eduyrs),
#'             formula2 = ~twoincomes+s(age)+s(eduyrs),
#'             data = DALSM_IncomeData)
#' plot(fit$derr)  ## Plot the estimated error density
#' print(fit$derr) ## ... and provide summary statistics for it
#'
#' ## Example 2: density estimation from censored income data
#' data(DALSM_IncomeData)
#' resp = DALSM_IncomeData[,1:2]
#' head(resp,n=20)
#' temp = Dens1d(y=resp,ymin=0) ## Create Dens1d object from positive censored data
#' obj = densityLPS(temp) ## Density estimation from IC & RC data
#' print(obj) ## Summary information on the estimated density
#' plot(obj,hist=TRUE) ## Visualize the estimated density
#' legend("topright",col=c("black","grey"),lwd=c(2,20),
#'        legend=c("Fitted density","Pseudo-data"),bty="n")
#'
#' ## Example 3: density estimation from simulated RC and IC data
#' ## Data generation
#' set.seed(123)
#' n = 500 ## Sample size
#' x = rgamma(n,10,2) ## Exact (unobserved) data
#' width = runif(n,1,3) ## Width of the IC data (mean width = 2)
#' w = runif(n) ## Positioning of the exact data within the interval
#' xmat = cbind(pmax(0,x-w*width),x+(1-w)*width) ## Generated IC data
#' t.cens = rexp(n,1/15) ## Right-censoring values
#' idx.RC = (1:n)[t.cens<x] ## Id's of the right-censored units
#' xmat[idx.RC,] = cbind(t.cens[idx.RC],Inf) ## Data for RC units: (t.cens,Inf)
#' head(xmat,15)
#' ## Density estimation with mean and variance constraints
#' obj.data = Dens1d(xmat,ymin=0) ## Prepare the data for estimation
#' obj = densityLPS(obj.data,Mean0=10/2,Var0=10/4) ## Density estimation
#' print(obj)
#' plot(obj) ## Plot the estimated density
#' curve(dgamma(x,10,2), ## ... and compare it to the true density (in red)
#'       add=TRUE,col="red",lwd=2,lty=2)
#' legend("topright",col=c("black","red"),lwd=c(2,2),lty=c(1,2),
#'        legend=c("Estimated density","True density"),bty="n")
#' ## Same story for the cdf
#' with(obj, curve(pdist(x),ymin,ymax,lwd=2,xlab="",ylab="F(x)"))
#' curve(pgamma(x,10,2),add=TRUE,col="red",lwd=2,lty=2)
#' legend("right",col=c("black","red"),lwd=c(2,2),lty=c(1,2),
#'        legend=c("Estimated cdf","True cdf"),bty="n")
#'
#' ## Example 4: density estimation from weighted grouped data
#' set.seed(123)
#' n = 500 ## Sample size
#' y = rgamma(n,10,2) ## Generate raw data
#' brks = c(0,2,4,6,10,12) ## Limits of the grouped data
#' y.cat = cut(y,breaks=brks) ## Categorization of the raw data
#' ymat = cbind(low=head(brks,-1), up=brks[-1]) ## Different grouped data
#' w = table(y.cat) ## Weight associated to each group data
#' print(cbind(ymat,w)) ## Observed grouped data
#' obj.data = Dens1d(y=ymat,w=w,ymin=0,pen.order=2) ## Data object
#' obj = densityLPS(obj.data,Mean0=10/2,Var0=10/4)  ## Density estimation
#' print(obj)
#' plot(obj,lwd=2)
#' weighted.hist(rowMeans(ymat),breaks=brks,weight=w,freq=FALSE,
#'               main="Estimated density",border=TRUE,col="#D9D9D980",add=TRUE)
#' curve(dgamma(x,10,2),add=TRUE,col="red",lwd=2,lty=2)
#' legend("topright",col=c("black","red","grey"),lwd=c(2,2,10),lty=c(1,2,1),
#'        legend=c("Estimated density","True density","Observed grouped data"),bty="n")
#'
densityLPS3 = function(obj.data,
                    is.density=TRUE,Mean0=NULL, Var0=NULL,
                    fixed.penalty=FALSE,method=c("LPS","Schall"),
                    fixed.phi=FALSE,phi.ref=NULL,
                    phi0=NULL,tau0=exp(5),tau.min=.1,
                    Lnorm=c("Linf","L2"),
                    grad.tol=1e-2, RDM.tol=1e-2, fun.tol=1e-3, itermax=101,
                    verbose=FALSE, criterion=c("AIC","dev","L2"), crit.tol=1e-2){
  method = match.arg(method)
  Lnorm = match.arg(Lnorm)
  criterion = match.arg(criterion)
  if (fixed.penalty) method = "Fixed penalty"
  if (is.null(obj.data)){
        print("Data object prepared using <Dens1d> required !!")
        return(NA)
  }
  ##
  nconstraints = (is.density) + (!is.null(Mean0)) + (!is.null(Var0)) ## Number of active moment constraints
  ##
  y = obj.data$y ## Reported times ; a nx2 matrix if IC data
  event = obj.data$event
  w = obj.data$w ## Weights
  sw = sum(w) ## Sum of the weights
  ymid = obj.data$ymid ## y if non-IC or .5*(y[,1]+y[,2]) if IC
  ## n.IC = obj.data$n.IC ## Number of Interval Censored data
  K = obj.data$K ## Number of B-splines coefs
  Pd = obj.data$Pd ## Penalty matrix
  pen.order = obj.data$pen.order ## Penalty order
  Bg = obj.data$Bgrid ## B-spline basis at the bin midpoints
  Bb = obj.data$Bbins ## B-spline basis at the bin limits
  bins = obj.data$bins ## Small bin limits
  nbins = obj.data$nbins ## Number of small bins
  C = obj.data$C ## (n x ngrid) Matrix <C>: indicates in which bin the event or censoring occured ;
  library(Matrix)
  CC <- as(C, "dgCMatrix")
  ##  For RC data, all the possible bins are indicated.
  fg = obj.data$fgrid ## Number of events within a given bin
  rg = obj.data$rgrid ## Number of units at risk within a given bin
  du = obj.data$dbins ## Width of a small bin
  ugrid = obj.data$ugrid ## Mid of small bins
  ngrid = length(ugrid) ## Number of small bins
  n = obj.data$n ## Number of data entries
  ev = obj.data$ev ## Eigenvalues
  ##
  ## Initial value for <phi> & <tau>
  T = with(obj.data, sum(w * ymid)) ## Total weighted "follow-up" duration
  dd = with(obj.data, sum(w * event)) ## Total weighted number of events (i.e. of non right-censored units)
  ##
  ## Initial values for the spline parameters
  mean.y = sum((w/sw) * ymid) ; sd.y = sqrt(sum((w/sw) * (ymid - mean.y)^2) * (sw / (sw-1)))
  ##  mean.y = mean(ymid) ; sd.y = sd(ymid)
  haz.t = (1/sd.y)*dt((ugrid-mean.y)/sd.y,df=5)/(1-pt((ugrid-mean.y)/sd.y,df=5))
  ## haz.norm = dnorm(ugrid,mean.y,sd.y)/(1-pnorm(ugrid,mean.y,sd.y))
  if (is.null(phi0)) phi.cur = phi0 = c(solve(t(Bg)%*%Bg)%*%t(Bg)%*%log(haz.t)) ## phi.cur = rep(log(dd/T),K)+runif(K,-.1,.1) ##
  else phi.cur = phi0
  ##
  if (is.null(phi.ref)) phi.ref = 0*phi.cur
  nphi = length(phi.cur)
  tau = max(tau.min,abs(tau0))
  ##
  if (!fixed.penalty){
    if (verbose){
      if (method=="LPS") cat("Selection of <tau>: Evidence-based method\n")
      if (method=="Schall") cat("Selection of <tau>: Schall's-based method\n")
    }
  }
  ##
  eta.cur = c(Bg %*% phi.cur)
  h.cur = exp(eta.cur) ## Current hazard at grid midpoints
  expctd = rg*(h.cur*du) ## Fitted number of events per bin
  llik.cur = sum(fg*eta.cur-expctd) ## Log-lik
  quad.cur = sum((phi.cur-phi.ref)*c(Pd%*%(phi.cur-phi.ref)))
  lpost.cur = llik.cur +.5*(K-pen.order)*log(tau)-.5*tau*quad.cur ## Log-posterior
  ## dev = -2*sum(fg*log(h.cur)-expctd) ## Deviance
  ##
  ## Start estimation process
  ## ########################
  ptm <- proc.time() ## ---------->
  ##
  omega0 = omega1 = omega2 = NULL ## Lagrange multipliers for Mean and Variance constraints
  if (is.density) omega0 = 0
  if (!is.null(Mean0)) omega1 = 0
  if (!is.null(Var0)) omega2 = 0
  omega.cur = c(omega0,omega1,omega2) ## Current Lagrange multipliers
  L1 = function(x) sum(abs(x)) ## L1-norm
  L2 = function(x) sqrt(sum(x^2)) ## L2-norm
  ##
  ## Functions to compute a moment and its gradient for the density associated to spline parms <phi>
  ## -----------------------------------------------------------------------------------------------
  ## Moment of order <pow>
  Moment.fun = function(phi,pow){
    eta.grid = c(Bg %*% phi)
    h.grid = exp(eta.grid)
    H.grid = cumsum(h.grid*du) -.5*h.grid*du ## <---- remove .5*du !!
    ans = sum((ugrid^pow) * h.grid*exp(-H.grid) * du)
    return(ans)
  }
  Finfty.fun = function(phi) Moment.fun(phi,pow=0) ## \int f(x)dx
  Mean.fun = function(phi) Moment.fun(phi,pow=1) ## Mean function
  Var.fun = function(phi) return(Moment.fun(phi,pow=2)-(Moment.fun(phi,pow=1))^2) ## Variance function
  ## Analytical gradient of a moment constraint
  gradMoment.fun = function(phi,pow){
    eta.grid = c(Bg %*% phi)
    h.grid = exp(eta.grid)
    H.grid = cumsum(h.grid*du) -.5*h.grid*du ## <---- remove .5*du !!
    dens.grid = h.grid*exp(-H.grid)
    ##
    Hb.jk = matrixStats::colCumsums(h.grid * Bg * du)
    ## Hb.jk = apply(h.grid*Bg*du,2,cumsum)
    temp = (ugrid^pow)*dens.grid*du
    ans = c(t(Bg) %*% (temp*(1+.5*h.grid*du)) - t(Hb.jk) %*% temp) ## <---- remove .5*du !!
    return(ans)
  }
  ##
  ## Linearization of the mean and variance constraints
  ## --------------------------------------------------
  ##  F1(phi) = f1  <--> Mean(phi) = Mean0
  ##  F2(phi) = f2  <--> Var(phi) = Var0
  ##       become (linearly)
  ##  Hcal %*% phi = Const  with  Hcal = rbind(H1,H2)
  ##   with row vectors  H1=grad.F1(phi) and H2=grad.F2(phi)
  H0 = H1 = H2 = NULL ; F0 = F1 = F2 = NULL ; const0 = const1 = const2 = NULL
  if (is.density){
    H0 = gradMoment.fun(phi.cur,pow=0) ## Current grad.F0(phi.cur)
    F0 = Finfty.fun(phi.cur) ## Current F(infty) = F0(phi.cur)
    const0 = 1 - F0    ## Current constraint failure
  }
  if (!is.null(Mean0)){
    H1 = gradMoment.fun(phi.cur,pow=1) ## Current grad.F1(phi.cur)
    ## H1 = grad(Mean.fun,phi.cur) ## Current grad.F1(phi.cur)
    F1 = Mean.fun(phi.cur) ## Current mean = F1(phi.cur)
    const1 = Mean0 - F1    ## Current constraint failure
  }
  if (!is.null(Var0)){
    if (is.null(Mean0)){
      H1 = gradMoment.fun(phi.cur,pow=1) ## Current grad.F1(phi.cur)
      F1 = Mean.fun(phi.cur) ## Current mean = F1(phi.cur)
    }
    H2 = gradMoment.fun(phi.cur,pow=2) - 2*F1*H1 ## Current grad.F2(phi.cur)
    ## H2 = grad(Var.fun,phi.cur) ## Current grad.F2(phi.cur)
    F2 = Var.fun(phi.cur)  ## Current variance = F2(phi.cur)
    const2 = Var0 - F2     ## Current constraint failure
  }
  if (nconstraints > 0) {
    Hcal = rbind(H0,H1,H2) ## Linearizer of (F1,F2)(phi)
    Const = c(Hcal%*%phi.cur) + c(const0,const1,const2)
  } else Const = NULL
  ##
  ## Start iterations
  ## ################
  ## Necessary precomputation of some quantities when (fixed.phi == TRUE)
  if (fixed.phi){
    eta.cur = c(Bg %*% phi.cur) ; h.cur = exp(eta.cur) ; expctd = rg*(h.cur*du)
    BWB = t(Bg)%*%(expctd*Bg)
    U.phi = c(t(Bg) %*% (fg - expctd)) - tau*c(Pd%*%(phi.cur-phi.ref))
    quad = sum(diff(phi.cur,diff=pen.order)^2)
    denom = sum(sw*ev/(tau+sw*ev)) ## cf. Evidence based method
    ## denom = sum(n*ev/(tau+n*ev)) ## cf. Evidence based method
    ed = sum(t(solve(BWB+tau*Pd))*BWB) ## cf. Schall's method
  }
  ## --------------------------------------------------------
  ## -- Function: Lagrangian with its gradient and Hessian --
  ## --------------------------------------------------------
  E.step <- function(phi){
      ## -- Current estimate for hazard, cumulative hazard and density --
      eta.grid = c(Bg %*% phi)
      h.grid = exp(eta.grid)
      H.grid = cumsum(h.grid*du) -.5*h.grid*du
      dens.grid = h.grid*exp(-H.grid)
      ## -- E-step for IC data: revise <fgrid> and <rgrid> using current density estimate --
      ##    (fgrid = Nbr of events per bin & rgrid = Nbr at risk per bin)
      if (is.matrix(y)){ ## i.e. if IC data present
          ##   <fgrid>: grid frequencies (= weighted number of events within a given bin)
          ##   <rgrid>: weighted number of units 'still at risk' in a given bin defined by the grid
          temp = t(dens.grid*t(C)) ## (n x ngrid) matrix
          normCst = rowSums(temp)
          temp = ifelse(normCst==0,0,1/normCst) * temp
          temp = w * temp  ## Account for weights
          fg = colSums(temp[event==1,,drop=FALSE])
          temp2 = colSums(temp)
          rg = sw-c(0,cumsum(temp2)[-ngrid]) ## Weighted number of units 'still at risk' in a given bin
      }
      ans = list(fg=fg, rg=rg)
      return(ans)
  }
  ##
  ## E.step2 <- function(phi) {
  ##     eta.grid <- c(Bg %*% phi)
  ##     h.grid <- exp(eta.grid)
  ##     H.grid <- cumsum(h.grid * du) - 0.5 * h.grid * du
  ##     dens.grid <- h.grid * exp(-H.grid)
  ##     ##
  ##     start_idx = attr(C,"start")
  ##     end_idx = attr(C,"end")
  ##     n <- length(start_idx)
  ##     ngrid <- length(dens.grid)
  ##     temp <- matrix(0, n, ngrid)
  ##     ##
  ##     for (i in seq_len(n)) {
  ##         s <- start_idx[i]
  ##         e <- end_idx[i]
  ##         if (!is.na(s) && s <= e) {
  ##             row_vals <- dens.grid[s:e]
  ##             norm <- sum(row_vals)
  ##             if (norm > 0) {
  ##                 temp[i, s:e] <- (row_vals / norm) * w[i]
  ##             }
  ##         }
  ##     }
  ##     ##
  ##     fg <- colSums(temp[event == 1, , drop = FALSE])
  ##     temp2 <- colSums(temp)
  ##     rg <- sw - c(0, cumsum(temp2)[-ngrid])
  ##     ##
  ##     ans = list(fg = fg, rg = rg)
  ##     return(ans)
  ## }
  ##
  g.fun <- function(theta,Dtheta=TRUE){
      phi.cur = theta[1:nphi]
      ##
      ## -- Current estimate for hazard, cumulative hazard and density --
      eta.grid = c(Bg %*% phi.cur)
      h.grid = exp(eta.grid)
      H.grid = cumsum(h.grid*du) -.5*h.grid*du
      dens.grid = h.grid*exp(-H.grid)
      pi.grid = dens.grid*du ## / sum(dens.grid)
      ##
      ## -- Likelihood and log-posterior --
      expctd = rg*(h.grid*du)
      llik = sum(fg*eta.grid-expctd) ## LogLik for complete data
      ## llik.obs = ifelse(is.density, sum(w * log(c(C %*% pi.grid))), NULL) ## LogLik for observed data
      quad.cur = sum((phi.cur-phi.ref)*c(Pd%*%(phi.cur-phi.ref)))
      lcur = llik +.5*(K-pen.order)*log(tau) -.5*tau * quad.cur ## Log-posterior
      ##
      ## -- Gradient and Hessian (before constraints) --
      BWB = t(Bg)%*%(expctd*Bg) ## crossprod(sqrt(expctd)*Bg)
      ## U.phi = colSums(Bg * (fg - expctd)) - tau*c(Pd%*%(phi.cur-phi.ref))
      ## U.phi = c(t(Bg) %*% (fg - expctd)) - tau*c(Pd%*%(phi.cur-phi.ref))
      U.phi = crossprod(Bg, fg - expctd) - tau*c(Pd%*%(phi.cur-phi.ref))
      H.phi = -BWB - tau*Pd ## Hessian when no constraints
      Hes = H.phi
      grad = U.phi
      gcur = lcur
      ## -- Constraints --
      F0 = F1 = F2 = NULL
      ## Density constraint
      if (is.density){
        H0 = gradMoment.fun(phi.cur,pow=0) ## Current grad.F0(phi.cur)
        F0 = Finfty.fun(phi.cur) ## Current F(infty) = F0(phi.cur)
        const0 = 1 - F0    ## Current constraint failure
      }
      ## Mean constraint
      if (!is.null(Mean0)){
        H1 = gradMoment.fun(phi.cur,pow=1) ## Current grad.F1(phi.cur)
        ## H1 = grad(Mean.fun,phi.cur) ## Current grad.F1(phi.cur)
        F1 = Mean.fun(phi.cur) ## Current mean = F1(phi.cur)
        const1 = Mean0 - F1    ## Current constraint failure
      }
      ## Variance constraint
      if (!is.null(Var0)){
        if (is.null(Mean0)){
          H1 = gradMoment.fun(phi.cur,pow=1) ## Current grad.F1(phi.cur)
          F1 = Mean.fun(phi.cur) ## Current mean = F1(phi.cur)
        }
        H2 = gradMoment.fun(phi.cur,pow=2) - 2*F1*H1 ## Current grad.F2(phi.cur)
        ## H2 = grad(Var.fun,phi.cur) ## Current grad.F2(phi.cur)
        F2 = Var.fun(phi.cur)  ## Current variance = F2(phi.cur)
        const2 = Var0 - F2     ## Current constraint failure
      }
      ## -- Gradient and Hessian (after constraints) --
      if (nconstraints > 0) {
          omega.cur = theta[-(1:nphi)]
          Hcal = rbind(H0,H1,H2) ## Linearizer of (F1,F2)(phi)
          Const = c(Hcal%*%phi.cur) + c(const0,const1,const2)
          ## Extend Score & Hessian due to constraints (Lagrange multiplier method)
          Hes = rbind(cbind(H.phi,-t(Hcal)),
                      cbind(-Hcal,diag(0,nconstraints)))
          Hes = Hes + diag(1e-5,ncol(Hes)) ## Add ridge penalty to force invertibility of Hes
          grad = c(U.phi-c(t(Hcal)%*%omega.cur), c(-Hcal%*%phi.cur+Const))
          ## Lagrangian
          gcur = lcur - sum(omega.cur * c(Hcal%*%phi.cur - Const))
      } else Const = NULL
      ## -- Return: Lagrangian with its gradient and Hessian --
      edf = sum(solve(BWB + tau * Pd) * BWB) ##  edf = sum(diag(solve(BWB + tau * Pd, BWB)))
      dev = -2*llik ##-2*llik.obs
      AIC = dev + 2*edf
      grid = list(eta=eta.grid, h=h.grid, H=H.grid, dens=dens.grid, pi=pi.grid, expctd=expctd)
      ans = list(g=gcur, theta=theta, grad=grad, Hes=Hes, ## Strictly necessary for Newton-Raphson
                 llik=llik, dev=dev, AIC=AIC, ## Monitoring GoF
                 F0=F0, F1=F1, F2=F2, ## Monitoring constraints
                 BWB=BWB, edf=edf, grid=grid) ## BWB + EDF + Density related quantities at gridpoints
      return(ans)
  } ## End g.fun
  ## -----------------------------------------------------------------------------
  ## -- Algorithm for density estimation from right- and interval-censored data --
  ## -----------------------------------------------------------------------------
  theta.cur = c(phi.cur, omega.cur)
  g.cur = g.fun(theta.cur)
  niter = 0
  ## if (verbose) cat("Iteration",niter,
  ##                  ": tau =",round(tau,2),
  ##                  "; Deviance =",round(g.cur$dev,2),
  ##                  "; AIC =",round(g.cur$AIC,2),"\n\n")
  converged = FALSE
  while (!(converged | fixed.phi)){ ## Proceed with iterations if (!converged & !fixed.phi)
    niter = niter + 1
    ok.phi = FALSE ; niter.phi = 0
    phi.old = phi.cur
    AIC.old = g.cur$AIC
    dev.old = g.cur$dev
    ##
    ##
    ## -- (1) Iterative update of <phi> and <omega> (given penalty parm <tau>) --
    obj.NR = AFTcure::NewtonRaphson(g.fun, theta.cur, itermax=5,
                                    grad.tol=grad.tol, RDM.tol=RDM.tol, fun.tol=fun.tol,
                                    Lnorm=Lnorm, verbose=FALSE)
    ## obj.LM = LevMarq(g.fun, theta.cur, verbose=verbose)
                   ## grad.tol=grad.tol, RDM.tol=RDM.tol, Lnorm=Lnorm,
    ##
    ## g.fn = function(x) g.fun(x)$g
    ## gr.fn = function(x) g.fun(x)$grad
    ## Hes.fn = function(x) g.fun(x)$Hes
    ## bele = marqLevAlg(theta.cur, fn=g.fn, gr=gr.fn, hess=Hes.fn)
    niter.phi = obj.NR$iter
    theta.cur = obj.NR$theta
    phi.cur = theta.cur[1:nphi]
    U.phi = obj.NR$grad[1:nphi]
    ok.phi = obj.NR$converged
    ##
    g.cur = g.fun(theta.cur) ## Update g.cur after Newton-Raphson
    ## ed = g.cur$edf ## Current EDF
    ##
    ## -- (2) Iterative selection of penalty parameter <tau> (given <phi>) --
    if (!fixed.penalty){
        quad = sum(diff((phi.cur-phi.ref),diff=pen.order)^2)
        ok.tau = FALSE ; niter.tau = 0
        while(!ok.tau){
            niter.tau = niter.tau + 1
            tau.old = tau
            if (method=="LPS"){
                ## Evidence-based method
                denom = sum(sw*ev/(tau+sw*ev))
                tau = max(tau.min,denom/quad)
            }
            if (method=="Schall"){
                ## Schall's method
                ed = sum(t(solve(g.cur$BWB+tau*Pd))*g.cur$BWB) ## = sum(diag(solve(BWB+tau*Pd)%*%BWB))
                tau = max(tau.min,(ed-pen.order)/quad)
            }
            ok.tau = (abs(tau-tau.old)/tau) < 1e-2
        }
        ## Update g.cur after <tau> update
        g.cur = g.fun(theta.cur)
    } else {
        ok.tau = TRUE
        niter.tau = 0
        quad = sum(diff((phi.cur-phi.ref),diff=pen.order)^2)
    }
    ##
    ## -- (3) E-step --
    E = E.step(phi.cur) ; fg = E$fg ; rg = E$rg
    ##
    if (verbose){
        cat("Iteration",niter,"(",niter.phi,niter.tau,"): tau =",round(tau,2),
            "; Deviance =",round(g.cur$dev,2),
            "; EDF =",round(g.cur$edf,2),
            "; AIC =",round(g.cur$AIC,2),"\n")
        if (is.density) cat("  >>> int f(x)dx = :",g.cur$F0,"\n")
        if (nconstraints > 1) cat("  >>> Target moment value(s):",c(Mean0,Var0),
                                  "; Fitted value(s):",round(c(g.cur$F1,g.cur$F2),3),"\n")
        ## if (nconstraints > 0) cat("\n")
    }
    ## (4) Convergence criteria
    ## ------------------------
    ## converged = ((ok.phi & ok.tau) | (niter > itermax))
    ok.AIC = (abs(AIC.old - g.cur$AIC) < crit.tol)
    ok.dev = (abs(dev.old - g.cur$dev) < crit.tol)
    ok.dphi = ((L2(phi.old-phi.cur)/L2(phi.cur)) < .1*crit.tol)
    if (criterion=="AIC"){
        converged = (ok.AIC && ok.phi && ok.tau)  | (niter >= itermax)
    } else if (criterion=="dev"){
        converged = (ok.dev && ok.phi && ok.tau)  | (niter >= itermax)
    } else if (criterion=="L2"){
        converged = (ok.dphi && ok.phi && ok.tau)  | (niter >= itermax)
    }
    if ((niter >= itermax) & verbose){
        message("Convergence in density estimation did not occur after ",niter," iterations\n")
    }
  }
  elapsed.time <- (proc.time()-ptm)[1] ## <------------
  if (verbose) cat("\n-- Elapsed time:",elapsed.time,"--\n")
  ##
  ## Fitted hazard / cumulative hazard / density at the bin midpoints
  eta.grid = c(Bg %*% phi.cur)
  h.grid = exp(eta.grid)
  H.grid = cumsum(h.grid*du) -.5*du*h.grid ## <---- remove .5*du !!
  dens.grid = h.grid*exp(-H.grid)
  ## Fitted hazard / cumulative hazard / density at the bin limits
  h.bins = exp(c(Bb %*% phi.cur))
  H.bins = cumsum(c(0,.5*(h.bins[-nbins]+h.bins[-1])*du)) ## Trapeze integration
  dens.bins = h.bins*exp(-H.bins)
  ##
  expctd = rg*(h.grid*du)
  BWB = t(Bg)%*%(expctd*Bg) ## crossprod(sqrt(expctd)*Bg)
  ##
  ## Estimated functions in relation to the fitted density
  ## #####################################################
  ## (a) Estimated density function (as a function)
  ## ----------------------------------------------
  ddist = function(x,ddist.out=1e-12){
    out = (x<bins[1]) | (x>=tail(bins,n=1))
    ans = ddist.out + 0*x
    ##
    xx = x[!out]
    idx = (xx-bins[1])/du
    ilow = floor(idx)+1
    ##
    frac = idx%%1
    ## ans[!out] = pmax(0,splinefun(bins,dens.bins)(x[!out]))
    ans[!out] = dens.bins[ilow] + frac*(dens.bins[ilow+ceiling(frac)]-dens.bins[ilow])
    return(ans)
  }
  ## (b) Estimated cumulative hazard (as a function)
  ## -----------------------------------------------
  Hdist = function(x,Hdist.out=c(1e-12,Inf)){
    out = (x<bins[1]) | (x>tail(bins,n=1))
    ans = 0*x
    ans[x<bins[1]] = Hdist.out[1]
    ans[x>=tail(bins,n=1)] = Hdist.out[2]
    ##
    xx = x[!out]
    idx = (xx-bins[1])/du
    ilow = floor(idx)+1
    ##
    frac = idx%%1
    ans[!out] = H.bins[ilow] + frac*(H.bins[ilow+ceiling(frac)]-H.bins[ilow])
    ## ans[!out] = splinefun(bins,H.bins)(x[!out])
    return(ans)
  }
  ## (c) Estimated hazard (as a function)
  ## ------------------------------------
  hdist = function(x,hdist.out=c(1e-12,Inf)){
    out = (x<bins[1]) | (x>=tail(bins,n=1))
    ## ans = hdist.out + 0*x
    ans = 0*x
    ans[x<bins[1]] = hdist.out[1]
    ans[x>=tail(bins,n=1)] = hdist.out[2]
    ##
    xx = x[!out]
    idx = (xx-bins[1])/du
    ilow = floor(idx)+1
    ##
    frac = idx%%1
    ans[!out] = h.bins[ilow] + frac*(h.bins[ilow+ceiling(frac)]-h.bins[ilow])
    return(ans)
  }
  ##
  ## (d) Estimated cdf (as a function)
  ## ---------------------------------
  pdist = function(x) 1-exp(-Hdist(x))
  ##
  ## Estimated mean and variance
  mean.dist = sum(ugrid*du*dens.grid)
  var.dist = sum((ugrid-mean.dist)^2*du*dens.grid)
  ltau = ifelse(fixed.penalty, log(abs(tau0)), log(tau-tau.min))
  ##
  ## Prepare returned list <ans>
  ## ###########################
  ans = list()
  ans$converged = converged ## Convergence logical indicator
  ans$ddist = ddist ## Estimated density as a function
  ans$Hdist = Hdist ## Estimated cumulative hazard function
  ans$hdist = hdist ## Estimated hazard function
  ans$pdist = pdist ## Estimated cdf
  ans$phi = phi.cur ; ans$U.phi = U.phi ## Estimated spline parameters and Score function there
  ans$tau = tau ; ans$ltau = ltau
  ans$est = c(phi.cur,ltau)
  ans$fixed.phi = fixed.phi ## Indicates whether the spline parameters were fixed (i.e. not estimated)
  ans$phi.ref = phi.ref ## Reference value for <phi> with departure from it penalized
  ans$quad = quad
  ans$BWB = BWB ; ans$Prec = BWB+tau*Pd ; ans$Fisher = solve(BWB+tau*Pd)
  ans$ugrid = ugrid ; ans$bins = obj.data$bins ## Grid midpoints
  ans$du = ans$dbins = obj.data$dbins ## Width of small bins
  ans$h.grid = h.grid ## Hazard function at grid midpoints
  ans$H.grid = H.grid ## Cumulative hazard function at grid midpoints
  ans$dens.grid = dens.grid ## Density at grid midpoints
  ans$h.bins = h.bins ## Hazard function at bin limits
  ans$H.bins = H.bins ## Cumulative hazard function at bin limits
  ans$dens.bins = dens.bins ## Density at bin limits
  ans$expected = expctd ## Expected number of events within a small bin
  ans$is.density = is.density
  if (!is.null(Mean0)) ans$Mean0 = Mean0 ## Requested mean value
  if (!is.null(Var0)) ans$Var0 = Var0    ## Requested variance value
  ans$Finfty = Finfty.fun(phi.cur)  ## int f(x)dx
  ans$mean.dist = mean.dist ; ans$var.dist = var.dist ## Mean and variance of the fitted density
  ans$method = method ## LPS or Schall
  if (method=="LPS") ans$denom = denom
  if (method=="Schall") ans$ed = ed
  ans$ed = sum(diag(solve(BWB+tau*Pd)%*%BWB)) ## Effective number of parameters
  ans = c(ans,obj.data) ## Remind input object (including summary statistics, knots, density support, etc.)
  ans$iterations = niter
  ans$elapsed.time = elapsed.time
  ##
  return(structure(ans,class="densLPS"))
}
