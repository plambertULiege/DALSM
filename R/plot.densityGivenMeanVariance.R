plot.densHazardConstrainedFit = function(x,xlim=range(fit$bins),breaks=NULL,histRC=FALSE,xlab="",ylab="Density",main="",...){
  fit = x
  dens.grid = fit$ddist(fit$ugrid) ## Fitted density at bin midpoints
  ##
  if (is.null(breaks)) breaks="Sturges"
  temp = with(fit, hist(rep(ugrid,apply((C/apply(C,1,sum)),2,sum)),breaks=breaks,plot=FALSE))
  temp = with(fit, hist(rep(ugrid,apply((C/apply(C,1,sum)),2,sum)),breaks=breaks,freq=FALSE,
                        xlim=xlim,ylim=1.05*c(0,max(dens.grid,temp$density)),
                        xlab=xlab,main=main,col="grey"))
  n.RC = sum(fit$event==0)
  if (histRC & n.RC>0) temp = with(fit, hist(rep(ugrid,apply((C[event==0,]/apply(C[event==0,],1,sum)),2,sum)),breaks=breaks,freq=FALSE,col="gray87",add=T))
  with(fit, curve(ddist,lwd=2,add=T,n=501))
}
