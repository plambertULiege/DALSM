#' Object created by \code{Dens1d} to prepare for density estimation from censored data using \code{densityLPS}
#'
#' @description An object returned by function \code{\link{Dens1d}} to prepare for density estimation with given mean and variance from censored data using function \code{\link{densityLPS}}.
#'
#' @return A \code{Dens1d.object} is a list with the following elements:
#'
#' \itemize{
#' \item{\code{n} : \verb{ }}{total sample size.}
#' \item{\code{y} : \verb{ }}{a n-vector if no interval-censored data, a nx2 matrix otherwise (for interval-censored data: left and right limits of the interval ; for right-censored data: the right limit is set to Inf).}
#' \item{\code{event} : \verb{ }}{a n-vector of observation indicators (0: right-censored ; 1: exactly observed or interval-censored).}
#' \item{\code{ymin} : \verb{ }}{left limit of the variable support.}
#' \item{\code{ymax} : \verb{ }}{right limit of the variable support.}
#' \item{\code{is.uncensored} : \verb{ }}{boolean n-vector indicating if the corresponding <y> value is not right- or interval-censored.}
#' \item{\code{n.uncensored} : \verb{ }}{total number of uncensored data.}
#' \item{\code{is.IC} : \verb{ }}{boolean n-vector indicating if the corresponding \code{y} value is interval-censored.}
#' \item{\code{n.IC} : \verb{ }}{total number of interval-censored data.}
#' \item{\code{is.RC} : \verb{ }}{n-vector of booleans indicating if the corresponding \code{y} value is right-censored.}
#' \item{\code{n.RC} : \verb{ }}{total number of right-censored data.}
#' \item{\code{ylow, yup} : \verb{ }}{n-vector with the lower and upper limits of the interval data (when interval-censored). When \code{y} is exactly observed or right-censored, the two values are identical.}
#' \item{\code{ymid} : \verb{ }}{n-vector containing the mean of \code{y.low} and \code{y.up}.}
#' \item{\code{K} : \verb{ }}{number of B-splines in the basis used to model the log-hazard.}
#' \item{\code{knots} : \verb{ }}{vector of knots for the B-spline basis.}
#' \item{\code{pen.order} : \verb{ }}{penalty order.}
#' \item{\code{Pd} : \verb{ }}{penalty matrix in the P-spline model.}
#' \item{\code{nbins} : \verb{ }}{number of small bins used to partition the variable support.}
#' \item{\code{bins} : \verb{ }}{small bins of equal width used to partition the variable support (cf. binning).}
#' \item{\code{ugrid} : \verb{ }}{midpoints of the small bins.}
#' \item{\code{dbins} : \verb{ }}{width of a small bin.}
#' \item{\code{Bbins} : \verb{ }}{((nbins +1) x K)-matrix containing the B-spline basis evaluated at the bin limits.}
#' \item{\code{C} : \verb{ }}{(n x nbins) matrix of event or censoring indicators \eqn{C_{ij}} for unit 'i' and bin 'j'. For a unit with IC data, the bins with an non-empty intersection with the interval are indicated. When the unit is associated to a precise event time or censoring time in bin 'j', then \eqn{C_{ij}=1} and 0 for other bins.}
#' \item{\code{Bgrid} : \verb{ }}{(nbins x K)-matrix containing the B-spline basis evaluated at the bin midpoints.}
#' \item{\code{fgrid} : \verb{ }}{nbins-vector with the estimated density values at the bin midpoints.}
#' \item{\code{rgrid} : \verb{ }}{nbins-vector with the number of units 'still at risk' of an "event" in a given bin.}
#' \item{\code{BsB, ev} : \verb{ }}{technical quantities used in the estimation procedure, see the code in Dens1d.R for more details.}
#' }
#'
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references Lambert, P. (2021). Fast Bayesian inference using Laplace approximations
#' in nonparametric double additive location-scale models with right- and
#' interval-censored data.
#' \emph{Computational Statistics and Data Analysis}, 161: 107250.
#' <doi:10.1016/j.csda.2021.107250>
#'
#' @seealso \code{\link{Dens1d}}, \code{\link{densityLPS}}, \code{\link{densLPS.object}}
#'
#' @name Dens1d.object
NULL

