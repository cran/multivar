#' Q function.
#'
#'@keywords internal
function_Q <- function( x, y, A, b, L, AtA, Atb){
  N <- prod(dim(x));
  out <- function_f(y, A, b) + sum(matlab::reshape((x - y)*gradient_f(y, AtA, Atb),c(N,1))) + 0.5 * L * norm(x - y,type = "F")^2;
  return(out)
}
