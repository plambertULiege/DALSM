## Weighted histogram function (adapted from stats::hist.default and weights::wtd.hist)
## -------------------------------------------------------------------------------------
#' Histogram from weighted data
#' @description Computes a histogram from weighted data.
#'   If \code{plot=TRUE}, the resulting object of class \code{histogram} is
#'   plotted before it is returned.
#' @usage weighted.hist(x, weight = NULL, breaks = "Sturges", freq = NULL,
#' probability = !freq, include.lowest = TRUE, right = TRUE, fuzz = 1e-07,
#' density = NULL, angle = 45, col = "lightgray", border = NULL,
#' main = paste("Histogram of", xname), xlim = range(breaks), ylim = NULL,
#' xlab = xname, ylab, axes = TRUE, plot = TRUE, labels = FALSE, nclass = NULL,
#' warn.unused = TRUE, ...)
#'
#' @param x a vector of values for which the histogram is desired.
#' @param weight optional vector of the same length as \code{x} defining the weights associated each each element in \code{x}. By default, it is \code{rep(1,length(x))}.
#' @param breaks one of:
#' \itemize{
#'    \item a vector giving the breakpoints between histogram cells,
#'    \item a single number giving the number of cells for the histogram,
#'    \item a character string naming an algorithm to compute the number of cells (see \sQuote{Details}),
#'    \item a function to compute the number of cells.
#' }
#' @param freq logical; if \code{TRUE}, the histogram graphic is a representation of frequencies, the \code{counts} component of the result; if \code{FALSE}, probability densities, component density, are plotted (so that the histogram has a total area of one). Defaults to \code{TRUE} \emph{if and only if} \code{breaks} are equidistant (and \code{probability} is not specified).
#' @param probability an \emph{alias} for \code{!freq}, for S compatibility.
#' @param include.lowest logical; if \code{TRUE}, an \code{x[i]} equal to the \code{breaks} value will be included in the first (or last, for \code{right = FALSE}) bar. This will be ignored (with a warning) unless \code{breaks} is a vector.
#' @param right logical; if \code{TRUE}, the histogram cells are right-closed (left open) intervals.
#' @param fuzz non-negative number, for the case when the data is \emph{pretty} and some observations \code{x[.]} are close but not exactly on a \code{break}. For counting fuzzy breaks proportional to \code{fuzz} are used. The default is occasionally suboptimal.
#' @param density the density of shading lines, in lines per inch. The default value of \code{NULL} means that no shading lines are drawn. Non-positive values of \code{density} also inhibit the drawing of shading lines.
#' @param angle the slope of shading lines, given as an angle in degrees (counter-clockwise).
#' @param col a colour to be used to fill the bars.
#' @param border the color of the border around the bars.  The default is to use the standard foreground color.
#' @param main main title.
#' @param xlim range of \code{x} for plotting.
#' @param ylim range of \code{y} values.
#' @param xlab label for the horizontal axis.
#' @param ylab label for the vertical axis.
#' @param axes logical.  If \code{TRUE} (default), axes are draw if the plot is drawn.
#' @param plot logical.  If \code{TRUE} (default), a histogram is plotted, otherwise a list of breaks and counts is returned.  In the latter case, a warning is used if (typically graphical) arguments are specified that only apply to the \code{plot=TRUE} case.
#' @param labels logical or character.  Additionally draw labels on top of bars, if not \code{FALSE}; see \code{plot.histogram} in the \code{graphics} package.
#' @param nclass numeric (integer).  For S(-PLUS) compatibility only, \code{nclass} is equivalent to \code{breaks} for a scalar or character argument.
#' @param warn.unused logical.  If \code{plot=FALSE} and \code{warn.unused=TRUE}, a warning will be issued when graphical parameters are passed to \code{hist.default()}.
#' @param ... further arguments and graphical parameters passed to \code{plot.histogram} and thence to title and axis (if \code{plot=TRUE}).
#'
#' @return a list comparable to what \code{hist.default} produces and a plot the weighted histogram if \code{plot=TRUE}.
#' @author Philippe Lambert \email{p.lambert@uliege.be}. The code for this function is fully based on \code{hist.default} from the graphics package in R, with minor modifications to account for weights.
#' @export
#'
#' @examples
#' ## Example 1
#' set.seed(123)
#' mround <- function(x,base) base*round(x/base) ## Rounding function
#' y = mround(rgamma(500,10,2), base=.5) ## Rounded data
#' tab = table(y) ; print(tab) ## Empirical distribution of the rounded data
#' ## Histogram from the empirical distribution
#' res = weighted.hist(as.numeric(names(tab)), weight=c(tab), xlab="",main="")
#'
#' ## Example 2
#' ## Generate categorized data
#' set.seed(123)
#' brks = c(0,3,5,7,9,15) ; ycat = cut(rgamma(500,10,2),breaks=brks)
#' tab = table(ycat) ## Empirical distribution of the categorized data
#' print(tab)
#' ymid = .5 * (brks[-6]+brks[-1]) ; w = c(tab) ## Midpoints and their weights
#' ## Histogram from weighted data
#' res = weighted.hist(ymid,weight=w,breaks=brks,xlab="",main="")
weighted.hist <- function (x, weight = NULL, breaks = "Sturges", freq = NULL, probability = !freq,
                           include.lowest = TRUE, right = TRUE, fuzz = 1e-07, density = NULL,
                           angle = 45, col = "lightgray", border = NULL,
                           main = paste("Histogram of", xname),
                           xlim = range(breaks), ylim = NULL, xlab = xname,
                           ylab, axes = TRUE, plot = TRUE, labels = FALSE, nclass = NULL,
                           warn.unused = TRUE, ...)
{
    if (!is.numeric(x))
        stop("'x' must be numeric")
    if (is.null(weight))
        weight <- rep(1, length(x))
    xname <- deparse1(substitute(x), collapse = "\n")
    ## n <- length(x <- x[is.finite(x)])
    ## n <- as.integer(n)
    weight <- weight[is.finite(x)]
    n <- sum(weight)
    x <- x[is.finite(x)]
    ##
    if (is.na(n))
        stop("invalid length(x)")
    use.br <- !missing(breaks)
    if (use.br) {
        if (!missing(nclass))
            warning("'nclass' not used when 'breaks' is specified")
    }
    else if (!is.null(nclass) && length(nclass) == 1L)
        breaks <- nclass
    use.br <- use.br && (nB <- length(breaks)) > 1L
    if (use.br)
        breaks <- sort(breaks)
    else {
        if (!include.lowest) {
            include.lowest <- TRUE
            warning("'include.lowest' ignored as 'breaks' is not a vector")
        }
        if (is.character(breaks)) {
            breaks <- match.arg(tolower(breaks), c("sturges",
                                                   "fd", "freedman-diaconis", "scott"))
            breaks <- switch(breaks, sturges = nclass.Sturges(x),
                             `freedman-diaconis` = , fd = nclass.FD(x), scott = nclass.scott(x),
                             stop("unknown 'breaks' algorithm"))
        }
        else if (is.function(breaks)) {
            breaks <- breaks(x)
        }
        if (length(breaks) == 1) {
            if (!is.numeric(breaks) || !is.finite(breaks) ||
                breaks < 1L)
                stop("invalid number of 'breaks'")
            if (breaks > 1e+06) {
                warning(gettextf("'breaks = %g' is too large and set to 1e6",
                                 breaks), domain = NA)
                breaks <- 1000000L
            }
            breaks <- pretty(range(x), n = breaks, min.n = 1)
            nB <- length(breaks)
            if (nB <= 1)
                stop(gettextf("hist.default: pretty() error, breaks=%s",
                              format(breaks)), domain = NA)
        }
        else {
            if (!is.numeric(breaks) || length(breaks) <= 1)
                stop(gettextf("Invalid breakpoints produced by 'breaks(x)': %s",
                              format(breaks)), domain = NA)
            breaks <- sort(breaks)
            nB <- length(breaks)
            use.br <- TRUE
        }
    }
    nB <- as.integer(nB)
    if (is.na(nB))
        stop("invalid length(breaks)")
    h <- as.double(diff(breaks))
    equidist <- !use.br || diff(range(h)) < 1e-07 * mean(h)
    if (!use.br && any(h <= 0))
        stop("'breaks' are not strictly increasing")
    freq1 <- freq
    if (is.null(freq)) {
        freq1 <- if (!missing(probability))
                     !as.logical(probability)
                 else equidist
    }
    else if (!missing(probability) && any(probability == freq))
        stop("'probability' is an alias for '!freq', however they differ.")
    stopifnot(`fuzz must be non-negative` = fuzz >= 0)
    diddle <- fuzz * if (nB > 5)
                         stats::median(h)
                     else if (nB <= 3)
                         diff(range(x))
                     else min(h[h > 0])
    fuzz <- if (right)
                c(if (include.lowest) -diddle else diddle, rep.int(diddle,
                                                                   nB - 1L))
            else c(rep.int(-diddle, nB - 1L), if (include.lowest) diddle else -diddle)
    fuzzybreaks <- breaks + fuzz
    ##counts <- .Call(C_BinCount, x, fuzzybreaks, right, include.lowest)
    counts <- as.numeric(xtabs(weight ~ cut(x, fuzzybreaks)))
    ##
    if (any(counts < 0L))
        stop("negative 'counts'. Internal Error.", domain = NA)
    ## if (sum(counts) < n)
    if (sum(counts) < n - .01)
        stop("some 'x' not counted; maybe 'breaks' do not span range of 'x'")
    dens <- counts/(n * h)
    mids <- 0.5 * (breaks[-1L] + breaks[-nB])
    r <- structure(list(breaks = breaks, counts = counts, density = dens,
                        mids = mids, xname = xname, equidist = equidist), class = "histogram")
    if (plot) {
        plot(r, freq = freq1, col = col, border = border, angle = angle,
             density = density, main = main, xlim = xlim, ylim = ylim,
             xlab = xlab, ylab = ylab, axes = axes, labels = labels,
             ...)
        invisible(r)
    }
    else {
        if (warn.unused) {
            nf <- names(formals())
            nf <- nf[is.na(match(nf, c("x", "breaks", "nclass", "plot",
                                       "include.lowest", "weight", "right", "fuzz")))]
            ## "plot", "include.lowest", "right", "fuzz")))]
            missE <- lapply(nf, function(n) substitute(missing(.),
                                                       list(. = as.name(n))))
            not.miss <- !vapply(missE, eval, NA, envir = environment())
            if (any(not.miss))
                warning(sprintf(ngettext(sum(not.miss), "argument %s is not made use of",
                                         "arguments %s are not made use of"), paste(sQuote(nf[not.miss]),
                                                                                    collapse = ", ")), domain = NA)
        }
        r
    }
}
