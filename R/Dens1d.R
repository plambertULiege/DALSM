## Author: Philippe LAMBERT (ULiege, UCLouvain, Belgium), Sept 2017
###################################################################################
#' Object creation for density estimation from right- or interval-censored data
#' @description Object creation for density estimation from right- or interval-censored data using function \link{densityLPS}.
#'
#' @usage Dens1d(y, event=NULL, ymin=NULL, ymax=NULL,
#'        K=25, equid.knots=TRUE, pen.order=2, nbins=501)
#' @param y a n-vector (if no interval-censored data) or a nx2 matrix (left and right limits of the interval for IC data ; right limit set to Inf for right-censored data).
#' @param event a n-vector of observation indicators (0: right-censored ; 1: exactly observed or interval-censored).
#' @param ymin left limit of the variable support.
#' @param ymax right limit of the variable support.
#' @param K number of B-splines in the basis to approximate the log-hazard.
#' @param equid.knots logical indicating if equidistants knots are desired.
#' @param pen.order penalty order when equidistant knots (otherwise: penalty matrix computed to penalize the second derivative).
#' @param nbins number of small bins used for quadrature and approximations.
#'
#' @return A \link{Dens1d.object}, i.e. a list with summary measures and precomputed components required for density estimation using \code{\link{densityLPS}}.
#' @export
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references Lambert, P. (2021). Fast Bayesian inference using Laplace approximations
#' in nonparametric double additive location-scale models with right- and
#' interval-censored data.
#' \emph{Computational Statistics and Data Analysis}, 161: 107250.
#' <doi:10.1016/j.csda.2021.107250>
#'
#' @seealso \code{\link{Dens1d.object}}, \code{\link{densityLPS}}.
#'
#' @examples
#' library(DALSM)
#' data(DALSM_IncomeData)
#' resp = DALSM_IncomeData[,1:2]
#' head(resp,n=20)
#' temp = Dens1d(y=resp,ymin=0) ## Create Dens1d object from positive censored data
#' obj = densityLPS(temp) ## Density estimation from IC & RC data
#' plot(obj) ## Visualize the estimated density
#'
Dens1d = function(y, event=NULL, ymin=NULL, ymax=NULL,
                  K=25, equid.knots=TRUE, pen.order=2, nbins=501){
  if (is.data.frame(y)) y = as.matrix(y)
  ## n: number of units
  if (is.matrix(y)){
    if (any(y[,1] > y[,2])) stop("y[i,1] cannot be larger than y[i,2] !!")
    n = nrow(y)
  } else if (is.vector(y)) {
    n = length(y)
  } else {
    stop("'y' should be a matrix or a vector !")
  }
  ##
  if (is.matrix(y)){ ## IC data with (n x 2) matrix with left and right limit
    if (is.null(event)){
      event = ifelse(is.infinite(y[,2]),0,1) ## Set what 'event' indicators are supposed to be
    }
    ## For RC data with y[i,2]=Inf, set yup[i]=ylow[i]
    ylow = y[,1] ; yup = ifelse(is.infinite(y[,2]),y[,1],y[,2])
    ymid = apply(cbind(ylow,yup),1,mean) ## Impute center of the interval for IC data
    is.IC = (ylow!=yup) ## Interval censoring indicator
    n.IC = sum(is.IC) ## Number of interval-censored data
  } else {
    ylow = yup = ymid = y
    is.IC = rep(FALSE,length(ymid))
    n.IC = 0
    ## Unless indicated, for a vectorial response, set event = 1
    if (is.null(event)) event = rep(1,n)
  }
  ## n = length(ymid) ## Number of units
  ## if (is.null(event)) event = rep(1,n) ##Unless indicated, assumed that event occured for each unit
  is.uncensored = (ylow==yup) & (event==1)
  n.uncensored = sum(is.uncensored)
  is.RC = (event==0)
  n.RC = sum(is.RC)
  ##
  ## Knots & Penalty matrix
  obj.knots = qknots(cbind(ylow,yup), xmin=ymin,xmax=ymax,
                     equid.knots = equid.knots, pen.order=pen.order, K=K)
  knots = obj.knots$knots ## Knots
  ymin = min(knots) ; ymax = max(knots)
  Pd = obj.knots$Pd ## Penalty matrix
  pen.order = obj.knots$pen.order ## Penalty order can change if based on quantiles
  ##
  ## Grid on (ymin,ymax) for binning
  bins = seq(min(knots),max(knots),length=nbins+1)
  dbins = bins[2] - bins[1]
  ugrid = bins[-1] - .5*dbins ## Grid midpoints
  ngrid = length(ugrid)
  Bgrid = Bsplines(ugrid, knots) ## B-spline basis on a fine regular grid
  Bbins = Bsplines(bins, knots) ## B-spline basis at the bin limits
  ## Compute the number of events & number at risk within a given bin
  ##   with IC data assimilated to their midpoint
  fgrid = hist(ymid[event==1],breaks=bins,plot=FALSE)$counts ## Grid frequencies (= number of events within a given bin)
  temp = hist(ymid,breaks=bins,plot=FALSE)$counts ## Count the number of units with an event or right censoring within a bin
  rgrid = n-c(0,cumsum(temp)[-ngrid]) ## Number of units 'still at risk' in a given bin defined by the grid
  rgrid = rgrid - .5*fgrid ## Account for the fact that a person is at risk half the time when having the event
  ##
  ## (n x ngrid) Matrix 'C' of event or censoring indicators
  ##   For a unit with IC data, the bins with an non-empty intersection with the interval are indicated.
  ##   Unit with a precise event time or censoring time have the containing bin indicated.
  C = matrix(0,nrow=n,ncol=ngrid)
  C = 0 + (1-outer(ylow,bins[-1],">=")) * outer(yup,bins[-length(bins)],">")
  ## C = as(C,"sparseMatrix")
  ##
  ## Other data summaries
  Bs = Bsplines(sort(ymid), knots)
  BsB = t(Bs)%*%Bs
  ## ev = eigenvalues of  diag(1/sqrt(svd(Pd)$d)) %*% (1/n)*M %*% diag(1/sqrt(svd(Pd)$d))
  sP = svd(Pd)
  id1 = 1:(K-pen.order) ; id0 = (K-pen.order+1):K
  Bstilde = Bs %*% sP$u
  B1 = Bstilde[,id1] ; B0 = Bstilde[,id0]
  Mat = t(B1)%*%B1-t(B1)%*%B0%*%solve(t(B0)%*%B0+1e-6*diag(ncol(B0)))%*%t(B0)%*%B1
  ev=svd(diag(sqrt(1/sP$d[id1]))%*%((1/n)*Mat)%*%diag(sqrt(1/sP$d[id1])))$d
  ##
  ans = list(n=n, y=y, event=event, ymin=ymin,ymax=ymax,
             is.uncensored=is.uncensored, n.uncensored=n.uncensored,
             is.IC=is.IC, n.IC = n.IC,
             is.RC=is.RC, n.RC=n.RC,
             ylow=ylow, yup=yup, ymid=ymid,
             K=K, knots=knots, pen.order=pen.order, Pd=Pd,
             nbins=nbins, bins=bins, ugrid=ugrid, dbins=dbins, Bbins=Bbins,
             C=C,
             Bgrid=Bgrid, fgrid=fgrid, rgrid=rgrid,
             BsB=BsB,
             ev=ev
  )
  return(structure(ans,class="Dens1d"))
}
