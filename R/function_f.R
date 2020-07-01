#' Objective function for fista sparse.
#'
#'@keywords internal
function_f <- function( x, A, b ){
  out <- 0.5 * norm(A %*% x - b, type = "F")^2
  return(out)
}
