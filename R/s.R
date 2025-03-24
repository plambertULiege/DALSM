#' Internal function to specify smooth terms in DALSM formulas
#' @keywords internal
#' @export
#' @param x Name of the variable for which an additive term is requested.
#' @param ref (Optional) reference value for \code{x} where the additive term is zero.
#' @return The submitted variable for which an additive term is required.
s <- function(x,ref=NULL) x
