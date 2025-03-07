#' Bayesian Information Criterion (BIC) of a DALSM object.
#'
#' @description
#' Bayesian Information Criterion (BIC) for the fitted DALSM model in a \code{DALSM.object}.
#'
#' @usage \method{BIC}{DALSM}(object, ...)
#'
#' @param object A \code{\link{DALSM.object}}.
#' @param ... Other optional DALSM objects.
#'
#' @details Bayesian (Schwarz) information for the fitted model in a DALSM object, with a penalty calculated using the total effective degrees of freedom and the total number of observed or interval-censored data, -2log(L) + log(n.uncensored)*ED.tot, smaller values being preferred during model selection.
#'
#' @return The BIC of the fitted DALSM model in \code{x}.
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
#' BIC(fit)
#'
#' @seealso \code{\link{DALSM}}, \code{\link{DALSM.object}}, \code{\link{BIC.DALSM}}, \code{\link{logEvid}}
#'
#' @export
#'
BIC.DALSM <- function(object, ...){
    obj = object
    lls = function(obj) return(ans = c(dev=obj$dev, edf=obj$ED.tot, d=obj$n.uncensored))
    if (!missing(...)) {
        vals = sapply(list(obj,...), lls)
        val <- data.frame(edf = round(vals[2L, ],2), BIC = vals[1L, ] + log(vals[3L, ]) * vals[2L, ])
        nos <- na.omit(vals[3L, ])
        if (length(nos) && any(nos != nos[1L])) warning("models are not all fitted to the same number of observations")
        Call <- match.call()
        Call$k <- NULL
        row.names(val) <- as.character(Call[-1L])
        val
    } else {
        vals = unname(lls(obj))
        vals[1L] + log(vals[3L]) * vals[2L]
    }
}
