#' Estimate a Sparse Multiple-Subject Vector Autoregression (VAR) Model
#'
#' Function for estimating multiple-subjbect Vector Autoregression (VAR) models
#' using Fast Iterative Shrinkage-Thresholding Algorithm (FISTA; Beck and Teboulle, 2009)
#'
#' @param A An N x P design matrix.
#' @param b An N x P outcome matrix.
#' @param lambda Regularization parameter.
#' @param x_true Numeric matrix containing the true transition matrix (if available).
#' @param niter Integer giving the maximum number of iterations.
#' @param backtrack Logical. If backtracking should be used in the FISTA algorithm.
#' @param w Numeric matrix containing the weights (if available).
#' @param conv Convergance criterion.
#' @details
#'
#' \strong{Function Under Development}
#'
#' This is a prototype function and is currently under development.
#'
#' @references
#'
#' Fisher, Z.F., Kim, Y., and Pipiras, V. (Under Review) Penalized
#' Estimation and Forecasting of Multiple Subject Intensive Longitudinal Data.
#'
#' Beck A. and Teboulle, M. (2009). A Fast Iterative Shrinkage-Thresholding
#' Algorithm for Linear Inverse Problems. SIAM J. Img. Sci. 2, 1, 183–202.
#'
#' @keywords var lasso
#'
#' @examples
#'
#' theta    <- matrix(rnorm(9),3,3)
#' data     <- var_sim(20, theta, diag(.1,3))
#' datalag  <- embed(data, 2)
#' b        <- datalag[,1:3]
#' A        <- datalag[,4:6]
#' fista_sparse(A, b, 1, theta, niter = 1, backtrack = TRUE)
#'
#'
#' @export
fista_sparse <- function(A, b, lambda, x_true, niter, backtrack, w = NULL, conv=1e-10){

  if(is.null(x_true)){
    x_true <- matrix(0, ncol(A), ncol(A))
  }

  out <- out.cpu <- out.obj <- out.relerr <- c()
  t   <- tnew <- 1
  x   <- xnew <- y <- array(0, dim(x_true))

  if (backtrack == 1){
    L <- norm(A, "2")^2/5 # L is stepsize or t in V&B
    eta <- 2
  } else {
    L <- norm(A, "2")^2
  }

  AtA <- t(A) %*% A
  Atb <- t(A) %*% b
  t0  <- Sys.time()

  for(i in 1:niter){

    if (backtrack == 1){

      L_bar <- L;
      found <- FALSE;

        while (found == FALSE){

            prox <- function_proximal(y, A, b, L_bar, lambda, AtA, Atb, w);

            if (function_f(prox, A, b) <= function_Q(prox, y, A, b, L_bar, AtA, Atb)){
                found <- TRUE;
            } else {
                L_bar <- L_bar * eta;
            }

        }

        L = L_bar;

    }


    x <- xnew;
    xnew <- prox;
    t <- tnew;
    tnew <- (1 + sqrt(1 + 4*t^2))/2;
    y    <- xnew + (t - 1) / tnew * (xnew - x);


    out.cpu   <- c(out.cpu, Sys.time()-t0)
    out.obj   <- c(out.obj, function_f(xnew, A, b) + function_g(xnew, lambda, w))
    out.relerr <- c(out.relerr, norm(xnew - x_true,'F')/norm(x_true,'F'))

    if(i > 1){
      if(abs((out.obj[i]-out.obj[i-1])/out.obj[i]) < conv){
        break
      }
    }

  }

  out.x <- xnew

  return(list(
    out.cpu = out.cpu,
    out.obj = out.obj,
    out.relerr = out.relerr,
    out.x = out.x
  ))
}

