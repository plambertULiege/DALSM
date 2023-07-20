## Log-determinant of a positive-definite matrix
##
## Input:
##  x: positive definite matrix
##
## Output:
##  log(det(x))
logDet.fun = function(x){
  if(any(is.nan(x))) return(NA)
  if(any(!is.finite(x))) return(NA)
  eival = svd(x)$d
  idx = which(eival > 0)
  return(sum(log(eival[idx])))
}
