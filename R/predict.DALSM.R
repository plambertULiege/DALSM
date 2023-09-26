#' Prediction based on a DALSM model
#' @description Estimated conditional mean and standard deviation of the response based on a DALSM object
#'  for given covariate values in a data frame 'newdata'. Conditional quantiles can also
#'  be computed.
#'
#' @usage
#' \method{predict}{DALSM}(object, newdata, probs, ...)
#'
#' @param object a \code{\link{DALSM.object}}.
#' @param newdata an optional data frame in which to look for variables with which to predict. If omitted, the covariate values in the original data frame used to fit the DALSM model are considered.
#' @param probs probability levels of the requested conditional quantiles.
#' @param ... further arguments passed to or from other methods.
#'
#' @return Returns a list containing:
#' \itemize{
#' \item{\code{mu} : \verb{ }}{estimated conditional mean.}
#' \item{\code{sd} : \verb{ }}{estimated conditional standard deviation.}
#' \item{\code{quant} : \verb{ }}{estimated quantiles (at probability level \code{probs}) of the fitted conditional response in the DALSM model.}
#' \item{\code{qerr} : \verb{ }}{quantiles (at probability level \code{probs}) of the fitted error distribution in the DALSM model.}
#' \item{\code{probs} : \verb{ }}{a reminder of the requested probability levels for the fitted quantiles.}
#' }
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
#' data2 = data.frame(age=c(40,60),eduyrs=c(18,12))
#' predict(fit, data = DALSM_IncomeData, newdata=data2, probs=c(.2,.5,.8))
predict.DALSM = function(object, newdata=NULL, probs, ...){
    data = object$data ## Data frame used when calling the DALSM function
    if (is.null(newdata)){
        ## message("The data frame <newdata> for which 'predictions' are desired is missing !")
        ## return(NA)
      newdata = NULL
      npred = nrow(data) ## Just provide fitted values to original data
    } else {
      npred = nrow(newdata) ## Give predictions for new covariate values
    }
    ##
    ## Values for <mu> for the <newdata>
    temp = DesignFormula(object$formula1, plyr::rbind.fill(data,newdata))
    if (!identical(object$regr1$knots,temp$knots)){
        message("Predictions requested for covariate values outside their range in the fitted model for location !")
        return(NA)
    }
    Xcal.1 = tail(temp$Xcal, npred)
    mu.pred = Xcal.1 %*% object$psi1
    ## Values for <sd> for the <newdata>
    temp = DesignFormula(object$formula2, plyr::rbind.fill(data,newdata))
    if (!identical(object$regr2$knots,temp$knots)){
        message("Predictions requested for covariate values outside their range in the fitted model for dispersion !")
        return(NA)
    }
    Xcal.2 = tail(temp$Xcal, npred)
    sd.pred = exp(Xcal.2 %*% object$psi2)
    ## Quantiles
    if (!missing(probs)){
        if (!(all(probs>=0) & all(probs<=1))){
            message("The requested quantiles should be in (0,1) !")
            return(NA)
        }
        nquant = length(probs)
        qerr = numeric(nquant) ; names(qerr) = probs
        quant.pred = matrix(nrow=npred,ncol=nquant)
        colnames(quant.pred) = probs
        low = min(object$derr$knots) ; up = max(object$derr$knots) ; mid = .5*(low+up)
        for (k in 1:nquant){
            qerr[k] = optim(mid, fn = function(y) abs(object$derr$pdist(y)-probs[k]),
                            lower=low, upper=up,
                            method="L-BFGS-B")$par
            quant.pred[,k] = mu.pred + qerr[k]*sd.pred
        }
    }
    ans = list(mu=c(mu.pred), sd=c(sd.pred))
    if (!missing(probs)) ans = c(ans, list(quant=quant.pred, qerr=qerr, probs=probs))
    return(ans)
}



