###################################################################################
## Author: Philippe LAMBERT (ULiege, UCLouvain, Belgium), Oct 2018
###################################################################################
#' Print a summary of the information in a \code{densLPS.object}
#' @description Print summary information on the density estimate obtained by \code{densityLPS} from censored data with given mean and variance.
#' @usage \method{print}{densLPS}(x,...)
#' @param x a \code{\link{densLPS.object}}.
#' @param ... Optional additional print parameters.
#'
#' @return No returned value (just printed summary).
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references Lambert, P. (2021). Fast Bayesian inference using Laplace approximations
#' in nonparametric double additive location-scale models with right- and
#' interval-censored data.
#' \emph{Computational Statistics and Data Analysis}, 161: 107250.
#' <doi:10.1016/j.csda.2021.107250>
#'
#' @seealso \code{\link{densLPS.object}}, \code{\link{plot.densLPS}}, \code{\link{densityLPS}}.
#'
#' @export
#'
#' @examples
#' require(DALSM)
#' data(DALSM_IncomeData)
#' resp = DALSM_IncomeData[,1:2]
#' fit = DALSM(y=resp,
#'             formula1 = ~twoincomes+s(age)+s(eduyrs),
#'             formula2 = ~twoincomes+s(age)+s(eduyrs),
#'             data = DALSM_IncomeData)
#' plot(fit$derr)  ## Plot the estimated error density
#' print(fit$derr) ## ... and provide some descriptive elements on it
#'
print.densLPS = function(x,...){
  fit = x
  cat("-----------------------------------------------------------------------\n")
  cat("   Constrained Density/Hazard estimation from censored data using LPS  \n")
  cat("-----------------------------------------------------------------------\n")
  cat("INPUT:\n")
  cat("  Total sample size: ",fit$n,"\n")
  cat("  Uncensored data: ",fit$n.uncensored," (",round(100*fit$n.uncensored/fit$n,2)," percents)\n",sep="")
  cat("  Interval-censored (IC) data: ",fit$n.IC," (",round(100*sum(fit$n.IC)/fit$n,2)," percents)\n",sep="")
  cat("  Right-censored (RC) data: ",sum(1-fit$event)," (",round(100*sum(1-fit$event)/fit$n,2)," percents)\n",sep="")
  cat("  ---\n")
  if (fit$n.uncensored>0){
    a = with(fit,min(ylow[is.uncensored])) ; b = with(fit,max(ylow[is.uncensored]))
    cat("  Range of the uncensored data: (",a,",",b,")\n",sep="")
  }
  if (fit$n.IC>0){
    a = with(fit,min(y[is.IC,])) ; b = with(fit,max(y[is.IC,]))
    cat("  Range of the IC data: (",a,",",b,")\n",sep="")
  }
  if (fit$n.RC>0){
    a = with(fit,min(ylow[is.RC])) ; b = with(fit,max(ylow[is.RC]))
    cat("  Range of the RC data: (",a,",",b,")\n",sep="")
  }
  cat("  ---\n")
  cat("  Assumed support: (",fit$ymin,",",fit$ymax,")\n",sep="")
  cat("  Number of small bins on the support:",fit$nbins,"\n")
  cat("  Number of B-splines:",fit$K,"; Penalty order:",fit$pen.order,"\n")
  cat("\nOUTPUT:\n")
  cat("  Returned functions:  ddist, pdist, hdist, Hdist(x)\n")
  cat("  Parameter estimates:  phi, tau\n")
  cat("  Value of the estimated cdf at +infty:",round(fit$Finfty,3),"\n")
  if (!is.null(fit$Mean0)) cat("  Constraint on the Mean:",fit$Mean0," ; Fitted mean:",fit$mean.dist,"\n")
  if (!is.null(fit$Var0)) cat("  Constraint on the Variance:",fit$Var0," ; Fitted variance:",fit$var.dist,"\n")
  if ((!fit$is.density) & is.null(fit$Mean0) & is.null(fit$Var0)) cat("  Unconstrained estimation\n")
  if (is.null(fit$Mean0)) cat("  Mean of the fitted density:",fit$mean.dist,"\n")
  if (is.null(fit$Var0))  cat("  Variance of the fitted density:",fit$var.dist,"\n")
  cat("  Selected penalty parameter <tau>:",round(fit$tau,1),"\n")
  cat("  Effective number of parameters:",round(fit$ed,1),"\n")
  cat("-----------------------------------------------------------------------\n")
  cat("Elapsed time: ",round(fit$elapsed.time,1)," seconds  (",fit$iter," iterations)\n",sep="")
  cat("-----------------------------------------------------------------------\n")
}

