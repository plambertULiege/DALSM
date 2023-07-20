#' Income data, see Lambert (2021, Section 4).
#'
#' @docType data
#'
#' @description Money available per person in 756 Belgian households (Lambert 2021, Section 4)
#' for respondents aged 25-55 when the main source of income
#' comes from wages or salaries. The data were extracted from the European Social Survey (2016)
#' and reworked to quantify the resources available at the individual level
#' in the selected Belgian households.
#'
#' @usage data(DALSM_IncomeData)
#'
#' @format A data frame with 756 rows and 5 columns.
#' \describe{
#'  \item{\code{inc.low}}{Lower bound for the individual income data.}
#'  \item{\code{inc.up}}{Upper bound for the individual income data (Inf if right-censored).}
#'  \item{\code{twoincomes}}{\code{1} Two-income household, \code{0} One-income household.}
#'  \item{\code{age}}{Age of the respondent}
#'  \item{\code{eduyrs}}{Number of years of education completed.}
#' }
#'
#' @references  Lambert, P. (2021). Fast Bayesian inference using Laplace approximations
#' in nonparametric double additive location-scale models
#' with right- and interval-censored data. Computational Statistics and Data Analysis.
#' https://doi.org/10.1016/j.csda.2021.107250
#'
#' @references European Social Survey Round 8 Data (2016)
#'  Data file edition 2.1. NSD - Norwegian Centre for Research Data, Norway.
#'
"DALSM_IncomeData"
