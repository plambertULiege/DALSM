###################################################################################
## GOAL
##  print, plot fitted density objects obtained for right censored data
##
## INPUT:
##  - list with class name "densHazardFit" correponding to output of densHazardFit
##
## OUTPUT:
##  - print: summary of data input and output content
##  - plot: histogram with fitted density
##
## Author: Philippe LAMBERT (ULg, UCL, Belgium), Oct 2018
###################################################################################
print.densHazardConstrainedFit = function(x,...){
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
  cat("  Estimated Hazard, cumulative hazard and density at grid midpoints <ugrid>\n")
  cat("  Returned functions:  ddist, pdist, hdist, Hdist(x)\n")
  message(sprintf("\n** %s iterations (%.1f seconds) **",fit$iter,fit$elapsed.time))
}
