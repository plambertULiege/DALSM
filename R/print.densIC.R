###################################################################################
## Author: Philippe LAMBERT (ULiege, UCLouvain, Belgium), Oct 2018
###################################################################################
#' Print a summary of the information in a \code{densIC.object}
#' @description Print summary information on the density estimate obtained by \code{densityIC} from censored data with given mean and variance.
#' @usage \method{print}{densIC}(x,...)
#' @param x a code{densIC.object}.
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
#' @seealso \code{\link{densIC.object}}, \code{\link{plot.densIC}}, \code{\link{densityIC}}.
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
print.densIC = function(x,...){
  fit = x
  cat("** Constrained Density/Hazard estimation from right- and Interval-censored data **\n")
  cat("INPUT:\n")
  message(sprintf("  Total sample size: %s",fit$n))
  message(sprintf("  Uncensored data: %s (%s percents)",fit$n.uncensored,round(100*fit$n.uncensored/fit$n,2)))
  message(sprintf("  Interval Censored data: %s (%s percents)",fit$n.IC,round(100*sum(fit$n.IC)/fit$n,2)))
  message(sprintf("  Right censored data: %s (%s percents)",sum(1-fit$event),round(100*sum(1-fit$event)/fit$n,2)))
  cat("  ---\n")
  if (fit$n.uncensored>0) message(sprintf("  Range of the uncensored data: (%.1e,%.1e)",with(fit,min(ylow[is.uncensored])),with(fit,max(ylow[is.uncensored]))))
  if (fit$n.IC>0) message(sprintf("  Range of the Interval Censored data: (%.1e,%.1e)",with(fit,min(y[is.IC,])),with(fit,max(y[is.IC,]))))
  if (fit$n.RC>0) message(sprintf("  Range of the right censored data: (%.1e,%.1e)",with(fit,min(ylow[is.RC])),with(fit,max(ylow[is.RC]))))
  cat("  ---\n")
  message(sprintf("  Assumed support: (%.1e,%.1e)",fit$ymin,fit$ymax))
  message(sprintf("  Number of small bins on the support: %s",fit$nbins))
  message(sprintf("  Number of B-splines: %s ; Penalty order: %s",fit$K,fit$pen.order))
  cat("\nOUTPUT:\n")
  message(sprintf("  Value of the estimated cdf at +infty: %.1e",fit$Finfty))
  if (!is.null(fit$Mean0)) message(sprintf("  Constraint on the Mean: %s ; Fitted Mean: %.1e",fit$Mean0,fit$mean.dist))
  if (!is.null(fit$Var0)) message(sprintf("  Constraint on the Variance: %s ; Fitted Variance: %.1e",fit$Var0,fit$var.dist))
  if ((!fit$is.density) & is.null(fit$Mean0) & is.null(fit$Var0)) message(sprintf("  Unconstrained estimation"))
  if (is.null(fit$Mean0)) message(sprintf("     Mean of the fitted density: %.1e",fit$mean.dist))
  if (is.null(fit$Var0))  message(sprintf("     Variance of the fitted density: %.1e",fit$var.dist))
  message(sprintf("  Selected penalty parameter <tau>: %.1f",fit$tau))
  message(sprintf("  Effective number of parameters: %.1f",fit$ed))
  cat("  Parameter estimates:  phi, tau\n")
  cat("  Returned functions:  ddist, pdist, hdist, Hdist(x)\n")
  message(sprintf("\n** %s iterations (%.1f seconds) **",fit$iter,fit$elapsed.time))
}

