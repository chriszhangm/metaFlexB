// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// main_draw
List main_draw(int M, arma::mat xct, arma::mat nct);
RcppExport SEXP _metaFlexB_main_draw(SEXP MSEXP, SEXP xctSEXP, SEXP nctSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type M(MSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xct(xctSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type nct(nctSEXP);
    rcpp_result_gen = Rcpp::wrap(main_draw(M, xct, nct));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_metaFlexB_main_draw", (DL_FUNC) &_metaFlexB_main_draw, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_metaFlexB(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
