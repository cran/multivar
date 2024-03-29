// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "multivar_types.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// fista_sparse
List fista_sparse(const mat& A, const mat& b, double lambda, const mat& x_true, int niter, bool backtrack, Rcpp::Nullable<Rcpp::NumericVector> w, double conv);
RcppExport SEXP _multivar_fista_sparse(SEXP ASEXP, SEXP bSEXP, SEXP lambdaSEXP, SEXP x_trueSEXP, SEXP niterSEXP, SEXP backtrackSEXP, SEXP wSEXP, SEXP convSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< const mat& >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const mat& >::type x_true(x_trueSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< bool >::type backtrack(backtrackSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type w(wSEXP);
    Rcpp::traits::input_parameter< double >::type conv(convSEXP);
    rcpp_result_gen = Rcpp::wrap(fista_sparse(A, b, lambda, x_true, niter, backtrack, w, conv));
    return rcpp_result_gen;
END_RCPP
}
// showValue
void showValue(int x);
RcppExport SEXP _multivar_showValue(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type x(xSEXP);
    showValue(x);
    return R_NilValue;
END_RCPP
}
// norm2
double norm2(NumericVector x);
RcppExport SEXP _multivar_norm2(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(norm2(x));
    return rcpp_result_gen;
END_RCPP
}
// ST1a
double ST1a(double z, double gam);
RcppExport SEXP _multivar_ST1a(SEXP zSEXP, SEXP gamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type z(zSEXP);
    Rcpp::traits::input_parameter< double >::type gam(gamSEXP);
    rcpp_result_gen = Rcpp::wrap(ST1a(z, gam));
    return rcpp_result_gen;
END_RCPP
}
// ST3a
colvec ST3a(colvec z, colvec gam);
RcppExport SEXP _multivar_ST3a(SEXP zSEXP, SEXP gamSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< colvec >::type z(zSEXP);
    Rcpp::traits::input_parameter< colvec >::type gam(gamSEXP);
    rcpp_result_gen = Rcpp::wrap(ST3a(z, gam));
    return rcpp_result_gen;
END_RCPP
}
// ind
uvec ind(int n2, int m);
RcppExport SEXP _multivar_ind(SEXP n2SEXP, SEXP mSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n2(n2SEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    rcpp_result_gen = Rcpp::wrap(ind(n2, m));
    return rcpp_result_gen;
END_RCPP
}
// FISTA
mat FISTA(const mat& Y, const mat& Z, mat& B, mat& W, const rowvec lambda1, const double eps, double step);
RcppExport SEXP _multivar_FISTA(SEXP YSEXP, SEXP ZSEXP, SEXP BSEXP, SEXP WSEXP, SEXP lambda1SEXP, SEXP epsSEXP, SEXP stepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< mat& >::type B(BSEXP);
    Rcpp::traits::input_parameter< mat& >::type W(WSEXP);
    Rcpp::traits::input_parameter< const rowvec >::type lambda1(lambda1SEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< double >::type step(stepSEXP);
    rcpp_result_gen = Rcpp::wrap(FISTA(Y, Z, B, W, lambda1, eps, step));
    return rcpp_result_gen;
END_RCPP
}
// lamloopFISTA
cube lamloopFISTA(NumericVector beta_, const mat& Y, const mat& Z, NumericVector W_, const mat& lambda1, const double eps, const colvec& YMean2, const colvec& ZMean2, mat& B1, double step);
RcppExport SEXP _multivar_lamloopFISTA(SEXP beta_SEXP, SEXP YSEXP, SEXP ZSEXP, SEXP W_SEXP, SEXP lambda1SEXP, SEXP epsSEXP, SEXP YMean2SEXP, SEXP ZMean2SEXP, SEXP B1SEXP, SEXP stepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type beta_(beta_SEXP);
    Rcpp::traits::input_parameter< const mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const mat& >::type Z(ZSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type W_(W_SEXP);
    Rcpp::traits::input_parameter< const mat& >::type lambda1(lambda1SEXP);
    Rcpp::traits::input_parameter< const double >::type eps(epsSEXP);
    Rcpp::traits::input_parameter< const colvec& >::type YMean2(YMean2SEXP);
    Rcpp::traits::input_parameter< const colvec& >::type ZMean2(ZMean2SEXP);
    Rcpp::traits::input_parameter< mat& >::type B1(B1SEXP);
    Rcpp::traits::input_parameter< double >::type step(stepSEXP);
    rcpp_result_gen = Rcpp::wrap(lamloopFISTA(beta_, Y, Z, W_, lambda1, eps, YMean2, ZMean2, B1, step));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_multivar_fista_sparse", (DL_FUNC) &_multivar_fista_sparse, 8},
    {"_multivar_showValue", (DL_FUNC) &_multivar_showValue, 1},
    {"_multivar_norm2", (DL_FUNC) &_multivar_norm2, 1},
    {"_multivar_ST1a", (DL_FUNC) &_multivar_ST1a, 2},
    {"_multivar_ST3a", (DL_FUNC) &_multivar_ST3a, 2},
    {"_multivar_ind", (DL_FUNC) &_multivar_ind, 2},
    {"_multivar_FISTA", (DL_FUNC) &_multivar_FISTA, 7},
    {"_multivar_lamloopFISTA", (DL_FUNC) &_multivar_lamloopFISTA, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_multivar(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
