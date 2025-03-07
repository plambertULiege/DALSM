#' Akaike Information Criterion (AIC) of a DALSM object.
#'
#' @description
#' Akaike Information Criterion (AIC) for the fitted DALSM model in a \code{DALSM.object}.
#'
#' @usage \method{AIC}{DALSM}(object, ..., k=2)
#'
#' @param object A \code{\link{DALSM.object}}.
#' @param k The penalty per parameter to be used. (Default: k=2 for the classical AIC).
#' @param ... Other optional DALSM objects.
#'
#' @details Akaike information criterion for the fitted model in a DALSM object, with a penalty calculated using the total effective degrees of freedom, -2log(L) + 2*ED.tot, smaller values being preferred during model selection.
#'
#' @return The AIC of the fitted DALSM model in \code{x}.
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references Time-varying exogenous covariates with frequently changing values in double additive cure survival model: an application to fertility.
#' \emph{Journal of the Royal Statistical Society, Series A}. <doi:10.1093/jrsssa/qnaf035>
#'
#' @examples
#' require(DALSM)
#' data(DALSM_IncomeData)
#' resp = DALSM_IncomeData[,1:2]
#' fit = DALSM(y=resp,
#'             formula1 = ~twoincomes+s(age)+s(eduyrs),
#'             formula2 = ~twoincomes+s(age)+s(eduyrs),
#'             data = DALSM_IncomeData)
#' AIC(fit)
#'
#' @seealso \code{\link{DALSM}}, \code{\link{DALSM.object}}, \code{\link{BIC.DALSM}}, \code{\link{logEvid}}
#'
#' @export
#'
AIC.DALSM <- function(object, ..., k=2){
    obj = object
    lls = function(obj) return(ans = c(dev=obj$dev, edf=obj$ED.tot, nobs=obj$sw))
    if (!missing(...)) {
        vals = sapply(list(obj,...), lls)
        val <- data.frame(edf = round(vals[2L, ],2), AIC = vals[1L, ] + k * vals[2L, ])
        nos <- na.omit(vals[3L, ])
        if (length(nos) && any(nos != nos[1L])) warning("models are not all fitted to the same number of observations")
        Call <- match.call()
        Call$k <- NULL
        row.names(val) <- as.character(Call[-1L])
        val
    } else {
        vals = unname(lls(obj))
        vals[1L] + k * vals[2L]
    }
}
