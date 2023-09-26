#' Income data
#'
#' @docType data
#'
#' @description Income (in euros) available per person in 756 Belgian households (ESS 2016 ; Lambert 2021, Section 4)
#' for respondents aged 25-55 when the main source of income
#' comes from wages or salaries. The data were extracted from the European Social Survey (2016)
#' and reworked to quantify the resources available at the individual level
#' in the selected Belgian households.
#'
#' @usage data(DALSM_IncomeData)
#'
#' @format A data frame with 756 rows and 5 columns:
#' \itemize{
#' \item{\code{inc.low} : \verb{ }}{lower bound for the individual income data.}
#' \item{\code{inc.up} : \verb{ }}{upper bound for the individual income data (\code{Inf} if right-censored).}
#' \item{\code{twoincomes} : \verb{ }}{\code{0}: one-income household ; \code{1}: two-income household.}
#' \item{\code{age} : \verb{ }}{age of the respondent.}
#' \item{\code{eduyrs} : \verb{ }}{number of years of education completed.}
#' }
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references Lambert, P. (2021). Fast Bayesian inference using Laplace approximations
#' in nonparametric double additive location-scale models with right- and
#' interval-censored data.
#' \emph{Computational Statistics and Data Analysis}, 161: 107250.
#' <doi:10.1016/j.csda.2021.107250>
#'
#' @references European Social Survey Round 8 Data (2016)
#'  Data file edition 2.1. NSD - Norwegian Centre for Research Data, Norway.
#'
"DALSM_IncomeData"
