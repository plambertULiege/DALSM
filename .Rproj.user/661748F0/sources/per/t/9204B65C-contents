### INPUT:
##  x: data that should be supported by the knots
##  xmin: minimum value for the knots
##  xmax: maximum value for the knots
##  equid.knots: TRUE if equidistant knots desired
##  pen.order: penalty order if  equid.knots = TRUE
##  K: number of B-splines in the basis
qknots = function(x, xmin=NULL,xmax=NULL, equid.knots = TRUE, pen.order=3, K=25){
  if (is.null(xmin)) xmin = min(x) - sd(x)
  if (is.null(xmax)) xmax = max(x) + sd(x)
  if (xmin > min(x)) xmin = min(x) - sd(x)
  if (xmax < max(x)) xmax = max(x) + sd(x)
  cnt = 0
  if (xmin < min(x)) cnt = cnt + 1
  if (xmax > max(x)) cnt = cnt + 1
  ##
  if (equid.knots){
    knots = seq(xmin,xmax,length=K-2)
    Dd = diff(diag(K),dif=pen.order) ## Penalty order
    Pd = t(Dd)%*%Dd #+ 1e-6*diag(KK)
  } else {
    knots = unique(c(xmin,quantile(x,p=seq(0,1,length=K-2-cnt)),xmax))
    ngrid = 501
    xgrid = seq(xmin,xmax,length=ngrid)
    d2B = D2Bsplines(xgrid,knots)
    Pd = 1/ngrid * (t(d2B)%*%d2B)
    pen.order = 2 ## Penalty based on 2nd derivatives of B-splines
  }
  ##
  ans = list(xmin=xmin,xmax=xmax,knots=knots,Pd=Pd,pen.order=pen.order)
  return(ans)
}
