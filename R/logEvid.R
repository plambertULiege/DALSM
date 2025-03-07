#' Generic function for computing log-evidence
#'
#' This is a generic function returning the log-evidence for the class of the input object.
#'
#' @param object An object for which log-evidence is to be computed.
#' @param ... Additional arguments (not used for now).
#' @export
logEvid <- function(object, ...) {
  UseMethod("logEvid")
}


#' Log-evidence of a DALSM object.
#'
#' @description
#' The log-evidence of the fitted DALSM model in a \code{DALSM.object}.
#'
#' @param object A \code{\link{DALSM.object}}.
#' @param ... Optionally more DALSM objects.
#'
#' @details Provides the log-evidence (or log-marginal likelihood) of the fitted DALSM model in a given \code{\link{DALSM.object}}, where the evidence is the marginal posterior of the penalty parameters at their selected values.
#' @return The log-evidence of the DALSM model in \code{object}.
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
#' logEvid(fit)
#'
#' @seealso \code{\link{DALSM}}, \code{\link{DALSM.object}}, \code{\link{AIC.DALSM}}, \code{\link{BIC.DALSM}}
#'
#' @export
#'
logEvid.DALSM <- function(object, ...){
    obj = object
  lls = function(obj) return(ans = c(logEvid=obj$logEvid, edf=obj$ED.tot, nobs=obj$sw))
  if (!missing(...)) {
    vals = sapply(list(obj,...), lls)
    val <- data.frame(edf = round(vals[2L, ],2), logEvid = vals[1L, ])
    nos <- na.omit(vals[3L, ])
    if (length(nos) && any(nos != nos[1L])) warning("models are not all fitted to the same number of observations")
    Call <- match.call()
    Call$k <- NULL
    row.names(val) <- as.character(Call[-1L])
    val
  } else {
    vals = unname(lls(obj))
    vals[1L]
  }
}


