#' Proximal operator.
#'
#'@keywords internal
function_proximal <- function( y, A, b, L, lambda, AtA, Atb, w){
  Y <- y - 1/L * gradient_f(y, AtA, Atb);
  out <- shrinkage(Y, 2*lambda/L, w);
  return(out)
}
