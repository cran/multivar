#' Shrinakge operator for fista sparse.
#'
#'@keywords internal
shrinkage <- function( y, tau, w ){
  if(!is.null(w)){
    tau <- w * tau
  }
  out <- sign(y)*pmax(abs(y) - tau, rep(0, length(0)))
  return(out)
}
