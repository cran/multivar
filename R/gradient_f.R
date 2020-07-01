#' Shrinakge operator for fista sparse.
#'
#'@keywords internal
gradient_f <- function( x, AtA, Atb ){
 out <- AtA %*% x - Atb
 return(out)
}
