// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// calc_M2_c
arma::mat calc_M2_c(int nSize, int nSpecies, arma::mat N, arma::vec scLinf, arma::mat ration, arma::mat wgt, arma::cube suitability, double phiMin, double otherFood);
RcppExport SEXP _LeMANS_calc_M2_c(SEXP nSizeSEXP, SEXP nSpeciesSEXP, SEXP NSEXP, SEXP scLinfSEXP, SEXP rationSEXP, SEXP wgtSEXP, SEXP suitabilitySEXP, SEXP phiMinSEXP, SEXP otherFoodSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nSize(nSizeSEXP);
    Rcpp::traits::input_parameter< int >::type nSpecies(nSpeciesSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type N(NSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type scLinf(scLinfSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type ration(rationSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type wgt(wgtSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type suitability(suitabilitySEXP);
    Rcpp::traits::input_parameter< double >::type phiMin(phiMinSEXP);
    Rcpp::traits::input_parameter< double >::type otherFood(otherFoodSEXP);
    rcpp_result_gen = Rcpp::wrap(calc_M2_c(nSize, nSpecies, N, scLinf, ration, wgt, suitability, phiMin, otherFood));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_LeMANS_calc_M2_c", (DL_FUNC) &_LeMANS_calc_M2_c, 9},
    {NULL, NULL, 0}
};

RcppExport void R_init_LeMANS(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
