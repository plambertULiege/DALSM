###################################################################################
## GOAL
##  Construct list with data, B-splines basis, Penalty matrix, etc.
##
## INPUT:
##  - y: data vector or matrix (if IC (=interval censored) data)
##  - ymin, ymax: support of the distribution
##  - event: vector of event indicators event[i] (0: right censored ; 1: event at y[i] or in (y[i,1],y[i,2]))
##  - K: number of B-splines in the basis
##  - pen.order: penalty order
##  - nbins: number of small bins
##
## OUTPUT:
##  Dens1d object = list with
##   - all inputs
##   - n.IC: number of interval censored data
##   - n.uncensored: number of uncensored data
##   - ylow: y if non-IC or y[,1] if IC data
##   - yup: y if non- IC or y[,2] if IC data
##   - ymid: y if non-IC or .5*(y[,1]+y[,2]) of IC data
##   - Dd, Pd: difference and penalty matrices of order <pen.order>
##   - bins: limits of the <nbins> small bins on the support
##   - ugrid: midpoints of the small bins
##   - dbins: width of a small bin
##   - Bbins: B-spline matrix at the bin limits
##   - fgrid: number of data per small bin
##   - rgrid: number of units 'still at risk' per small bin
##   - Bgrid: B-spline matrix at the small bin midpoints
##   - BsB: t(B)%*%B with B=B-spline matrix evaluated at <y>
##   - ev: eigenvalues of A'A with  A=(1/sqrt(n))*B*svd(Pd)$u*diag(1/svd(Pd)$d)
##
## Author: Philippe LAMBERT (ULg, UCL, Belgium), Sept 2017
###################################################################################
#' Creates an object for density estimation from right- or interval-censored data using \code{\link{densityGivenMeanVariance}}
#'
#' @param y a n-vector (if no interval-censored data) or a nx2 matrix (left and right limits of the interval-censored variable ; identical columns if right-censored)
#' @param event a n-vector of observation indicators (0: right-censored ; 1: exactly observed or interval-censored)
#' @param ymin left limit of the variable support
#' @param ymax right limit of the variable support
#' @param K number of B-splines in the basis to approximate the log-hazard
#' @param equid.knots logical indicating if equidistants knots are desired
#' @param pen.order penalty order when equidistant knots (otherwise: computed from the second derivative)
#' @param nbins number of small bins for quadrature and approximations
#'
#' @return A \code{Dens1d} object, i.e. a list with summary measures and precomputed components required for density estimation using \code{\link{densityGivenMeanVariance}}
#' @export
#'
#' @examples
Dens1d = function(y, event=NULL, ymin=NULL, ymax=NULL,
                  K=25, equid.knots=TRUE, pen.order=2, nbins=501){
  if (is.matrix(y)){ ## IC data with (n x 2) matrix with left and right limit
    ylow = y[,1] ; yup = y[,2]
    ymid = apply(y,1,mean) ## Impute center of the interval for IC data
    is.IC = (ylow!=yup) ## Interval censoring indicator
    n.IC = sum(is.IC) ## Number of interval-censored data
  } else {
    ylow = yup = ymid = y
    is.IC = rep(FALSE,length(ymid))
    n.IC = 0
  }
  n = length(ymid) ## Number of units
  if (is.null(event)) event = rep(1,n) ##Unless indicated, assumed that event occured for each unit
  is.uncensored = (ylow==yup) & (event==1)
  n.uncensored = sum(is.uncensored)
  is.RC = (event==0)
  n.RC = sum(is.RC)
  ##
  ## Knots & Penalty matrix
  obj.knots = qknots(y, xmin=ymin,xmax=ymax,
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
  ## (n x ngrid) Matrix <Cf> of event or censoring indicators
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
  ##
  ## return(ans)
  return(structure(ans,class="Dens1d"))
}
