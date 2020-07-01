#' Objective function g.
#'
#'@keywords internal
function_g <- function( x, lambda, w){
  if(is.null(w)){
    out <- lambda * sum(abs(x))
  } else {
    out <- lambda * sum(w*abs(x))
  }
  return(out)
}
