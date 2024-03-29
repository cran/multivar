#include <RcppArmadillo.h>
using namespace arma; 
using namespace Rcpp;
using std::exp;
using std::sqrt;
using std::log;
using arma::pow;
using arma::mat;
using arma::size;
using arma::sign;
using arma::max;
using arma::accu;
using arma::sum;
using arma::field;
using arma::norm;

//[[Rcpp::depends(RcppArmadillo)]]

void showMatrix(mat X) {
    Rcout << "Armadillo matrix is" << std::endl << X << std::endl;
}

mat mvrnormArma(int n, vec mu, mat sigma) {
   int ncols = sigma.n_cols;
   mat Y = arma::randn(n, ncols);
   return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}

mat shrinkage2(mat u, double tau, Rcpp::Nullable<Rcpp::NumericVector> w = R_NilValue) {
  mat z; z.zeros(size(u));
  mat taum(size(u));
  taum.fill(tau);
  
  if (w.isNotNull()) {
    //colvec wa = Rcpp::as<arma::colvec>(w);
    //mat wmat = arma::repmat(wa, 1,u.n_cols);
    mat wmat = Rcpp::as<arma::mat>(w);
    taum = taum % wmat;
  } 

  return sign(u) % max(abs(u)-taum, z);
}


double function_f2(mat x, mat A, mat b) {
  double out = 0.5 * pow(norm(A * x - b, "fro"),2);
  return(out);
}


mat gradient_f2(mat x, mat AtA, mat Atb) {
  mat out = AtA * x - Atb;
  return(out);
}


double function_g2(mat x, double lambda, Rcpp::Nullable<Rcpp::NumericVector> w = R_NilValue){
    mat xabs; xabs.zeros(size(x));
    if (w.isNotNull()) {
      //colvec wa = Rcpp::as<arma::colvec>(w);
      //mat wmat = arma::repmat(wa, 1,x.n_cols);
      mat wmat = Rcpp::as<arma::mat>(w);
      xabs = abs(x) % wmat;
    } else {
      xabs = abs(x);
    }
  return(lambda * accu(xabs));
}


double function_Q2(mat x, mat y, mat A, mat b, double L, mat AtA, mat Atb){
  double out = function_f2(y,A,b) + accu((x - y)%gradient_f2(y, AtA, Atb)) + 0.5 * L * pow(norm(x - y,"fro"),2);
  return(out);
}

mat function_proximal2(mat y, mat A, mat b, double L, double lambda, mat AtA, mat Atb,Rcpp::Nullable<Rcpp::NumericVector> w = R_NilValue){
  mat out = shrinkage2(y - 1/L * gradient_f2(y, AtA, Atb), 2*lambda/L,w);
  return(out);
}

mat var_sim(int n, mat coef, mat sigma) {
   // under development
   vec mu;  mu.zeros(coef.n_cols);
   mat errors = mvrnormArma(n, mu, sigma);
   
   int m = errors.n_cols; 

   mat simdata(n,m);
   simdata.row(0) = arma::zeros<arma::mat>(1,m);
   for (int row=1; row<m; row++) {
      simdata.row(row) = simdata.row(row-1)*trans(coef)+errors.row(row);
   }
   return simdata;
}

//[[Rcpp::export]]
List fista_sparse(const mat& A, const mat& b, double lambda, const mat& x_true, int niter, bool backtrack, Rcpp::Nullable<Rcpp::NumericVector> w = R_NilValue, double conv=1e-10){

  int index = 0;
  double  t = 1, tnew = 1, L = 1, eta = 0, L_bar = 0;
  mat x; x.zeros(size(x_true));
  mat xnew; xnew.zeros(size(x_true));
  mat y; y.zeros(size(x_true));
  mat prox; prox.zeros(size(x_true));
  
  vec obj;  obj.zeros(niter);
  vec err;  err.zeros(niter);

  if (backtrack){
    L = pow(norm(A, 2),2)/5;
    eta = 2.0;
  } else {
    L = pow(norm(A, 2),2);
  }

  const mat& AtA = A.t() * A;
  const mat& Atb = A.t() * b;
  const mat& xtrue = x_true;

  for (int i = 0; i < niter; i++) {

    if(backtrack){

      L_bar = L;
      bool found = false;

      while(!found){
      
      prox = function_proximal2(y, A, b, L_bar, lambda, AtA, Atb, w);
        
      if (function_f2(prox, A, b) <= function_Q2(prox, y, A, b, L_bar, AtA, Atb)){
        found = true;
      } else {
        L_bar = L_bar * eta;
      }

      }

      L = L_bar;

     }

     x = xnew;
     xnew = prox;
     t = tnew;
     tnew = (1 + sqrt(1 + 4*pow(t,2)))/2;
     y = xnew + (t - 1) / tnew * (xnew - x);

     obj(i) = function_f2(xnew, A, b) + function_g2(xnew, lambda, w);
     err(i) = norm(xnew - xtrue,"fro")/norm(xtrue,"fro");

    
    if(i > 0){
      if(fabs((obj(i)-obj(i-1))/obj(i)) < conv){
       index = i;
       break;
      }
    }

  }

  return(List::create(Named("out.x") = xnew.t() , Named("out.obj") = obj, Named("out.relerr") = err, Named("index") = index));
  
}




