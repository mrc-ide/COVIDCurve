// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// NatCubic_SplineGrowth_loglike_cubicspline
Rcpp::List NatCubic_SplineGrowth_loglike_cubicspline(Rcpp::NumericVector params, int param_i, Rcpp::List data, Rcpp::List misc);
RcppExport SEXP _COVIDCurve_NatCubic_SplineGrowth_loglike_cubicspline(SEXP paramsSEXP, SEXP param_iSEXP, SEXP dataSEXP, SEXP miscSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type params(paramsSEXP);
    Rcpp::traits::input_parameter< int >::type param_i(param_iSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type data(dataSEXP);
    Rcpp::traits::input_parameter< Rcpp::List >::type misc(miscSEXP);
    rcpp_result_gen = Rcpp::wrap(NatCubic_SplineGrowth_loglike_cubicspline(params, param_i, data, misc));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_COVIDCurve_NatCubic_SplineGrowth_loglike_cubicspline", (DL_FUNC) &_COVIDCurve_NatCubic_SplineGrowth_loglike_cubicspline, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_COVIDCurve(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
