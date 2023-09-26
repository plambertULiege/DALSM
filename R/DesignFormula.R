## -----------------------------------
## Philippe LAMBERT (ULiege, Oct 2018)
## Email:  p.lambert@uliege.be
## Web: http://www.statsoc.ulg.ac.be
## -----------------------------------
#' Internal function extracting design matrices from formulas in the DALSM function and computing penalty related matrices
#' @description Internal function extracting design matrices from formulas in the DALSM function and computing penalty related matrices.
#'
#' @param formula a formula describing the fixed effects and the additive terms in a regression model.
#' @param data a dataframe containing the data.
#' @param K number of B-splines to describe an additive term.
#' @param pen.order desired penalty order for the spline parameters in the additive terms.
#' @param n number of units (Default: number of rows in the design matrix constructed from the formula and the data frame).
#'
#' @return a list with
#' \itemize{
#' \item{\code{Z} : \verb{ }}{(n x nfixed) design matrix with fixed effects (including a first column of 1).}
#' \item{\code{X} : \verb{ }}{(n x J) design matrix with the covariates involved in the additive terms.}
#' \item{\code{nfixed} : \verb{ }}{number of fixed effect regression parameters.}
#' \item{\code{J} : \verb{ }}{number of additive terms.}
#' \item{\code{K} : \verb{ }}{number of B-splines in a basis used to estimate an additive term.}
#' \item{\code{Bx} : \verb{ }}{list with J objects (one per additive term) including (B,Dd,Pd,K,cm).}
#' \item{\code{Pd.x, Dd.x} : \verb{ }}{penalty and difference penalty matrices applied on the spline parameters of an additive term.}
#' \item{\code{knots.x} : \verb{ }}{list of length J with the knots associated to each of the J additive terms.}
#' \item{\code{Bcal} : \verb{ }}{column-stacked matrix with the J centered B-spline bases.}
#' \item{\code{Xcal} : \verb{ }}{Z column-stacked with the J centered B-spline bases to yield the full design matrix (with column labels).}
#' \item{\code{pen.order} : \verb{ }}{penalty order for the spline parameters in the additive terms.}
#' \item{\code{additive.lab} : \verb{ }}{labels for the columns in <Bcal> associated to the additive terms.}
#' \item{\code{lambda.lab} : \verb{ }}{labels for the penalty parameters.}
#' }
#' @author Philippe Lambert \email{p.lambert@uliege.be}
#' @references Lambert, P. (2021). Fast Bayesian inference using Laplace approximations
#' in nonparametric double additive location-scale models with right- and
#' interval-censored data.
#' \emph{Computational Statistics and Data Analysis}, 161: 107250.
#' <doi:10.1016/j.csda.2021.107250>
#'
#' @keywords internal
#'
DesignFormula = function(formula, data, K=10, pen.order=2, n=NULL){
  if (!inherits(formula, "formula"))
    stop("Incorrect model formula")
  if ((formula=="~1")&(missing(data))){
    if (is.null(n)){
      cat("Model with only the intercept: the sample size <n> or a data frame should be provided !\n")
      return(NULL)
    }
    XX = model.matrix(~ 1, data = data.frame(rep(1,n)))
  } else {
    ## Extract design matrix
    if (missing(data)) {
      mf <- stats::model.frame(formula)
      XX <- stats::model.matrix(mf)
    }
    else {
      mf <- stats::model.frame(formula, data = data)
      XX <- stats::model.matrix(mf, data = data)
    }
  }
  ## Identify additive terms
  smterms <- grepl("s(", colnames(XX), fixed = TRUE)
  X = subset(XX, select=smterms) ## Covariates with additive effects in the design matrix
  ## Reorder the design matrix to start with 'fixed effect' covariates
  idx = c(which(!smterms),which(smterms))
  XX <- subset(XX,select=idx) ## Reordering
  if (any(is.infinite(XX)))
    stop("Covariates contain Inf, NA or NaN values")
  J <- sum(smterms) ## Nbr of additive terms
  nfixed <- ncol(XX) - J
  n <- nrow(XX) ## Nbr of units
  ##
  Z <- subset(XX, select=1:nfixed) ## 'Fixed effects' part of the design matrix
  if (ncol(Z) == 1)
    colnames(Z) <- "(Intercept)"
  ## Some labels
  fixed.lab = colnames(Z) ## Labels of the fixed effects
  additive.lab = lambda.lab = Bcal = NULL
  ## Additives terms
  Bx = NULL ## B-spline for the additive terms
  if (J > 0){
    K = floor(K) ## Make sure that this is an integer
    ## Number of knots
    nknots = K-1 ## Number of knots for the cubic B-splines in the basis
    if (!is.vector(nknots, mode = "numeric") || length(nknots) > 1 || is.na(nknots))
      stop("nknots must be numeric of length 1")
    if (K < 5 || K > 60)
      stop("K must be between 5 and 60")
    ## Penalty order
    pen.order <- floor(pen.order)
    if (pen.order < 1 || pen.order > 4)
      stop("Penalty order must be between 1 and 4")
    knots.x = NULL
    for (j in 1:J){  ## Loop over functional components in location
      xj = XX[,nfixed+j]
      knots.x[[j]] = seq(min(xj),max(xj),length=nknots)
      ## B-spline matrix for the jth additive term
      Bx[[j]] = centeredBasis.gen(xj,knots=knots.x[[j]],cm=NULL,pen.order)
      colnames(Bx[[j]]$B) = paste(colnames(XX)[nfixed+j],".",1:K,sep="")
    }
    ## Global design matrix
    Bcal = NULL ## Global design matrix for the B-splines of the additive terms
    for (j in 1:J) Bcal = cbind(Bcal,Bx[[j]]$B)
    ## Extra labels
    # additive.lab = colnames(XX)[-(1:nfixed)]
    additive.lab = unname(sapply(colnames(XX)[-(1:nfixed)], function(x) substring(x,3,nchar(x)-1)))
    lambda.lab = paste("lambda.",additive.lab,sep="")
  }
  if (J > 0) {
    Xcal = cbind(Z,Bcal) ## Global design matrix
  } else Xcal = Z
  ##
  Pd.x = Dd.x = NULL
  if (J > 0){
    Pd.x = Bx[[1]]$Pd ## Same penalty matrix for all functional components
    Dd.x = Bx[[1]]$Dd ## Same difference matrix for all functional components
  }
  ##
  ans = list(Z=Z,X=X,Xcal=Xcal,
             nfixed=nfixed,J=J,K=K)
  if (J > 0){
    ans = c(ans,
            list(
                 Bcal=Bcal,
                 Bx=Bx,Pd.x=Pd.x,Dd.x=Dd.x,knots.x=knots.x,
                 pen.order=pen.order,additive.lab=additive.lab,lambda.lab=lambda.lab))
  }
  return(ans)
}
