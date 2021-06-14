// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// blb_mix
Rcpp::List blb_mix(Rcpp::NumericVector y, int b, int s, int r, int d, Rcpp::Nullable<Rcpp::NumericVector> pr_, Rcpp::Nullable<Rcpp::NumericVector> pi_, Rcpp::Nullable<Rcpp::NumericVector> mu_, Rcpp::Nullable<Rcpp::NumericVector> sd_, int max_iter, double tol, double beta, double c, Rcpp::Nullable<Rcpp::NumericVector> schedule_);
RcppExport SEXP _bumblb_blb_mix(SEXP ySEXP, SEXP bSEXP, SEXP sSEXP, SEXP rSEXP, SEXP dSEXP, SEXP pr_SEXP, SEXP pi_SEXP, SEXP mu_SEXP, SEXP sd_SEXP, SEXP max_iterSEXP, SEXP tolSEXP, SEXP betaSEXP, SEXP cSEXP, SEXP schedule_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type y(ySEXP);
    Rcpp::traits::input_parameter< int >::type b(bSEXP);
    Rcpp::traits::input_parameter< int >::type s(sSEXP);
    Rcpp::traits::input_parameter< int >::type r(rSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type pr_(pr_SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type pi_(pi_SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type mu_(mu_SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type sd_(sd_SEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type schedule_(schedule_SEXP);
    rcpp_result_gen = Rcpp::wrap(blb_mix(y, b, s, r, d, pr_, pi_, mu_, sd_, max_iter, tol, beta, c, schedule_));
    return rcpp_result_gen;
END_RCPP
}
// gmm
Rcpp::List gmm(Rcpp::NumericVector Y, int d, Rcpp::Nullable<Rcpp::NumericVector> pi_, Rcpp::Nullable<Rcpp::NumericVector> mu_, Rcpp::Nullable<Rcpp::NumericVector> sd_, int max_iter, double tol, double beta, double c, Rcpp::Nullable<Rcpp::NumericVector> schedule_);
RcppExport SEXP _bumblb_gmm(SEXP YSEXP, SEXP dSEXP, SEXP pi_SEXP, SEXP mu_SEXP, SEXP sd_SEXP, SEXP max_iterSEXP, SEXP tolSEXP, SEXP betaSEXP, SEXP cSEXP, SEXP schedule_SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< int >::type d(dSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type pi_(pi_SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type mu_(mu_SEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type sd_(sd_SEXP);
    Rcpp::traits::input_parameter< int >::type max_iter(max_iterSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< double >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double >::type c(cSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericVector> >::type schedule_(schedule_SEXP);
    rcpp_result_gen = Rcpp::wrap(gmm(Y, d, pi_, mu_, sd_, max_iter, tol, beta, c, schedule_));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_bumblb_blb_mix", (DL_FUNC) &_bumblb_blb_mix, 14},
    {"_bumblb_gmm", (DL_FUNC) &_bumblb_gmm, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_bumblb(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
