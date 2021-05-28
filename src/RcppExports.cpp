// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// blb_mix
List blb_mix(NumericVector y, int b, int s, int r, int d);
RcppExport SEXP _bumblb_blb_mix(SEXP ySEXP, SEXP bSEXP, SEXP sSEXP, SEXP rSEXP, SEXP dSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    rcpp_result_gen = Rcpp::wrap(blb_mix(y, b, s, r, d));
    return rcpp_result_gen;
END_RCPP
}
// gmm
List gmm(NumericVector X, int d, Nullable<NumericVector> pi_, Nullable<NumericVector> mu_, Nullable<NumericVector> sd_, int max_iter, double tol);
RcppExport SEXP _bumblb_gmm(SEXP XSEXP, SEXP dSEXP, SEXP pi_SEXP, SEXP mu_SEXP, SEXP sd_SEXP, SEXP max_iterSEXP, SEXP tolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericVector >::type X(XSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type pi_(pi_SEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type mu_(mu_SEXP);
    Rcpp::traits::input_parameter< Nullable<NumericVector> >::type sd_(sd_SEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    rcpp_result_gen = Rcpp::wrap(gmm(X, d, pi_, mu_, sd_, max_iter, tol));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bumblb_blb_mix", (DL_FUNC) &_bumblb_blb_mix, 5},
    {"_bumblb_gmm", (DL_FUNC) &_bumblb_gmm, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_bumblb(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
