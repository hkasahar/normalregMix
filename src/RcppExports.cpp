// This file was generated by Rcpp::compileAttributes
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// cppNormalmixPMLE
List cppNormalmixPMLE(NumericMatrix bs, NumericVector ys, NumericMatrix zs, NumericVector mu0s, NumericVector sigma0s, int m, int p, double an, int maxit, int ninits, double tol, double tau, int h, int k);
RcppExport SEXP normalregMix_cppNormalmixPMLE(SEXP bsSEXP, SEXP ysSEXP, SEXP zsSEXP, SEXP mu0sSEXP, SEXP sigma0sSEXP, SEXP mSEXP, SEXP pSEXP, SEXP anSEXP, SEXP maxitSEXP, SEXP ninitsSEXP, SEXP tolSEXP, SEXP tauSEXP, SEXP hSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type bs(bsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ys(ysSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type zs(zsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu0s(mu0sSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma0s(sigma0sSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type an(anSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< int >::type ninits(ninitsSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    __result = Rcpp::wrap(cppNormalmixPMLE(bs, ys, zs, mu0s, sigma0s, m, p, an, maxit, ninits, tol, tau, h, k));
    return __result;
END_RCPP
}
// cppRegmixPMLE
List cppRegmixPMLE(NumericMatrix bs, NumericVector ys, NumericMatrix xs, NumericMatrix zs, NumericVector mu0s, NumericVector sigma0s, int m, int p, double an, int maxit, int ninits, double tol, double tau, int h, int k);
RcppExport SEXP normalregMix_cppRegmixPMLE(SEXP bsSEXP, SEXP ysSEXP, SEXP xsSEXP, SEXP zsSEXP, SEXP mu0sSEXP, SEXP sigma0sSEXP, SEXP mSEXP, SEXP pSEXP, SEXP anSEXP, SEXP maxitSEXP, SEXP ninitsSEXP, SEXP tolSEXP, SEXP tauSEXP, SEXP hSEXP, SEXP kSEXP) {
BEGIN_RCPP
    Rcpp::RObject __result;
    Rcpp::RNGScope __rngScope;
    Rcpp::traits::input_parameter< NumericMatrix >::type bs(bsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type ys(ysSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type xs(xsSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type zs(zsSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type mu0s(mu0sSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type sigma0s(sigma0sSEXP);
    Rcpp::traits::input_parameter< int >::type m(mSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< double >::type an(anSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< int >::type ninits(ninitsSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< int >::type h(hSEXP);
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    __result = Rcpp::wrap(cppRegmixPMLE(bs, ys, xs, zs, mu0s, sigma0s, m, p, an, maxit, ninits, tol, tau, h, k));
    return __result;
END_RCPP
}
