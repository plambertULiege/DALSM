## GOAL:
##  Plot a credible region (in grey) for a curve together with its point estimate (lwd=2)
##
## INPUT:
##  x: vector of values where the curve is evaluated
##  mat: cbind(f.hat,f.low,f.up) is a matrix containing the point estimates
##       and the liming values for the credible region
##  add: if FALSE, creates graph from scratch, otherwise surperposes within the current plotting window
##
## OUTPUT:
##  Plot with credible region in grey and the point estimate for the curve in black
##
## -----------------------------------
## Philippe LAMBERT (ULiege, Oct 2018
## Email:  p.lambert@uliege.be
## Web: http://www.statsoc.ulg.ac.be
## -----------------------------------
plotRegion = function(x,mat,add=FALSE,xlim=range(x),ylim=range(mat),
                      lwd=2,xlab="",ylab="",main="",...){
  f = mat[,1] ; f.low = mat[,2] ; f.up = mat[,3]
  if (add==FALSE) plot(x,f,type="n",ylim=ylim,xlim=xlim,
                       lwd=lwd,xlab=xlab,ylab=ylab,main=main,...)
  polygon(c(x,rev(x)),c(f.low,rev(f.up)),col="grey",border=F)
  lines(x,f,lwd=lwd)
}
