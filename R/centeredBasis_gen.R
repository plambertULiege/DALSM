#' Generation of a recentered cubic B-spline basis matrix in additive models
#' @description Generation of a cubic B-spline basis matrix with recentered columns
#'  to handle the identifiability constraint in additive models. See Wood (CRC Press 2017, pp. 175-176) for more details.
#' @param x vector of values where to compute the "recentered" B-spline basis.
#' @param knots vector of knots (that should cover the values in <x>).
#' @param cm (Optional) values subtracted from each column of the original B-spline matrix.
#' @param pen.order penalty order for the B-spline parameters (Default: 2).
#'
#' @return List containing
#' \itemize{
#' \item{\code{B} : \verb{ }}{centered cubic B-spline matrix (with columns recentered to have mean 0 over equi-spaced x values on the range of the knots).}
#' \item{\code{Dd} : \verb{ }}{difference matrix (of order <pen.order>) for the associated centered B-spline matrix.}
#' \item{\code{Pd} : \verb{ }}{penalty matrix (of order <pen.order>) for the associated centered B-spline matrix.}
#' \item{\code{K} : \verb{ }}{number of centered B-splines in the basis.}
#' \item{\code{cm} : \verb{ }}{values subtracted from each column of the original B-spline matrix. By default, this is a vector containing the mean of each column in the original B-spline matrix.}
#'}
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references Lambert, P. (2021). Fast Bayesian inference using Laplace approximations
#' in nonparametric double additive location-scale models with right- and
#' interval-censored data.
#' \emph{Computational Statistics and Data Analysis}, 161: 107250.
#' <doi:10.1016/j.csda.2021.107250>
#'
#' @export
#'
#' @examples
#' x = seq(0,1,by=.01)
#' knots = seq(0,1,length=5)
#' obj = centeredBasis.gen(x,knots)
#' matplot(x,obj$B,type="l",ylab="Centered B-splines")
#' colMeans(obj$B)
centeredBasis.gen = function(x,knots,cm=NULL,pen.order=2){
  if ((max(x)>max(knots))|(min(x)<min(knots))){
    cat("The knots do no cover the values of x !!\n")
    return(NULL)
  }
  ##
  knots.x = knots ; pen.order.x = pen.order ## Only corrected on 1st May 2020 !!!
  ##
  temp = cubicBsplines::Bsplines(x, knots.x)
  idx = 1:(ncol(temp)-1)
  B = temp[,idx,drop=FALSE] ## Reduce the number of B-splines by 1 !!
  if (is.null(cm)){
    ## Consider a B-spline matrix over a dense equi-spaced set of values on range(knots)
    ##  to compute a recentering value for each column of the matrix
    B.ref = Bsplines(seq(min(knots),max(knots),length=100), knots.x)[,idx]
    cm = colMeans(B.ref)
  }
  B = sweep(B,2,cm) ## Subtract column mean from each column
  Dd = diff(diag(ncol(temp)),dif=pen.order.x)[,idx] ## Penalty order <---- Important to remove last column !!
  Pd = t(Dd) %*% Dd
  return(list(B=B,Dd=Dd,Pd=Pd,K=ncol(B),cm=cm))
}
