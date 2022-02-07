// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// irls_gaussian_identity
Rcpp::List irls_gaussian_identity(arma::mat A, arma::vec b, arma::vec eta, arma::vec offset);
RcppExport SEXP _glmhdfe_irls_gaussian_identity(SEXP ASEXP, SEXP bSEXP, SEXP etaSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type offset(offsetSEXP);
    rcpp_result_gen = Rcpp::wrap(irls_gaussian_identity(A, b, eta, offset));
    return rcpp_result_gen;
END_RCPP
}
// fe_gaussian_identity_foc
double fe_gaussian_identity_foc(arma::vec lhs, arma::vec theta);
RcppExport SEXP _glmhdfe_fe_gaussian_identity_foc(SEXP lhsSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type lhs(lhsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(fe_gaussian_identity_foc(lhs, theta));
    return rcpp_result_gen;
END_RCPP
}
// irls_gaussian_log
Rcpp::List irls_gaussian_log(arma::mat A, arma::vec b, arma::vec eta, arma::vec offset);
RcppExport SEXP _glmhdfe_irls_gaussian_log(SEXP ASEXP, SEXP bSEXP, SEXP etaSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type offset(offsetSEXP);
    rcpp_result_gen = Rcpp::wrap(irls_gaussian_log(A, b, eta, offset));
    return rcpp_result_gen;
END_RCPP
}
// fe_gaussian_log_foc
double fe_gaussian_log_foc(arma::vec lhs, arma::vec theta);
RcppExport SEXP _glmhdfe_fe_gaussian_log_foc(SEXP lhsSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type lhs(lhsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(fe_gaussian_log_foc(lhs, theta));
    return rcpp_result_gen;
END_RCPP
}
// irls_poisson_log
Rcpp::List irls_poisson_log(arma::mat A, arma::vec b, arma::vec eta, arma::vec offset);
RcppExport SEXP _glmhdfe_irls_poisson_log(SEXP ASEXP, SEXP bSEXP, SEXP etaSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type offset(offsetSEXP);
    rcpp_result_gen = Rcpp::wrap(irls_poisson_log(A, b, eta, offset));
    return rcpp_result_gen;
END_RCPP
}
// fe_poisson_log_foc
double fe_poisson_log_foc(arma::vec lhs, arma::vec theta);
RcppExport SEXP _glmhdfe_fe_poisson_log_foc(SEXP lhsSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type lhs(lhsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(fe_poisson_log_foc(lhs, theta));
    return rcpp_result_gen;
END_RCPP
}
// irls_gamma_log
Rcpp::List irls_gamma_log(arma::mat A, arma::vec b, arma::vec eta, arma::vec offset);
RcppExport SEXP _glmhdfe_irls_gamma_log(SEXP ASEXP, SEXP bSEXP, SEXP etaSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type offset(offsetSEXP);
    rcpp_result_gen = Rcpp::wrap(irls_gamma_log(A, b, eta, offset));
    return rcpp_result_gen;
END_RCPP
}
// fe_gamma_log_foc
double fe_gamma_log_foc(arma::vec lhs, arma::vec theta);
RcppExport SEXP _glmhdfe_fe_gamma_log_foc(SEXP lhsSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type lhs(lhsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(fe_gamma_log_foc(lhs, theta));
    return rcpp_result_gen;
END_RCPP
}
// irls_inverse_gaussian_log
Rcpp::List irls_inverse_gaussian_log(arma::mat A, arma::vec b, arma::vec eta, arma::vec offset);
RcppExport SEXP _glmhdfe_irls_inverse_gaussian_log(SEXP ASEXP, SEXP bSEXP, SEXP etaSEXP, SEXP offsetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type A(ASEXP);
    Rcpp::traits::input_parameter< arma::vec >::type b(bSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type offset(offsetSEXP);
    rcpp_result_gen = Rcpp::wrap(irls_inverse_gaussian_log(A, b, eta, offset));
    return rcpp_result_gen;
END_RCPP
}
// fe_inverse_gaussian_log_foc
double fe_inverse_gaussian_log_foc(arma::vec lhs, arma::vec theta);
RcppExport SEXP _glmhdfe_fe_inverse_gaussian_log_foc(SEXP lhsSEXP, SEXP thetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type lhs(lhsSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type theta(thetaSEXP);
    rcpp_result_gen = Rcpp::wrap(fe_inverse_gaussian_log_foc(lhs, theta));
    return rcpp_result_gen;
END_RCPP
}
// ols
Rcpp::List ols(arma::mat X, arma::vec y, bool full);
RcppExport SEXP _glmhdfe_ols(SEXP XSEXP, SEXP ySEXP, SEXP fullSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< bool >::type full(fullSEXP);
    rcpp_result_gen = Rcpp::wrap(ols(X, y, full));
    return rcpp_result_gen;
END_RCPP
}
// wls
Rcpp::List wls(arma::mat X, arma::vec y, arma::vec w, bool full);
RcppExport SEXP _glmhdfe_wls(SEXP XSEXP, SEXP ySEXP, SEXP wSEXP, SEXP fullSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    Rcpp::traits::input_parameter< bool >::type full(fullSEXP);
    rcpp_result_gen = Rcpp::wrap(wls(X, y, w, full));
    return rcpp_result_gen;
END_RCPP
}
// wdemean
arma::vec wdemean(arma::vec var, arma::vec weight);
RcppExport SEXP _glmhdfe_wdemean(SEXP varSEXP, SEXP weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type var(varSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type weight(weightSEXP);
    rcpp_result_gen = Rcpp::wrap(wdemean(var, weight));
    return rcpp_result_gen;
END_RCPP
}
// demean
arma::vec demean(arma::vec var);
RcppExport SEXP _glmhdfe_demean(SEXP varSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type var(varSEXP);
    rcpp_result_gen = Rcpp::wrap(demean(var));
    return rcpp_result_gen;
END_RCPP
}
// wdemean_list
Rcpp::List wdemean_list(Rcpp::List var, Rcpp::NumericVector weight);
RcppExport SEXP _glmhdfe_wdemean_list(SEXP varSEXP, SEXP weightSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type var(varSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type weight(weightSEXP);
    rcpp_result_gen = Rcpp::wrap(wdemean_list(var, weight));
    return rcpp_result_gen;
END_RCPP
}
// demean_list
Rcpp::List demean_list(Rcpp::List var);
RcppExport SEXP _glmhdfe_demean_list(SEXP varSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type var(varSEXP);
    rcpp_result_gen = Rcpp::wrap(demean_list(var));
    return rcpp_result_gen;
END_RCPP
}
// gradient
arma::mat gradient(arma::mat X, arma::vec score);
RcppExport SEXP _glmhdfe_gradient(SEXP XSEXP, SEXP scoreSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type score(scoreSEXP);
    rcpp_result_gen = Rcpp::wrap(gradient(X, score));
    return rcpp_result_gen;
END_RCPP
}
// vcov_estimate
Rcpp::List vcov_estimate(arma::mat X, arma::mat grad, arma::mat cluster_meat);
RcppExport SEXP _glmhdfe_vcov_estimate(SEXP XSEXP, SEXP gradSEXP, SEXP cluster_meatSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type grad(gradSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type cluster_meat(cluster_meatSEXP);
    rcpp_result_gen = Rcpp::wrap(vcov_estimate(X, grad, cluster_meat));
    return rcpp_result_gen;
END_RCPP
}
// demean_criterion
arma::mat demean_criterion(arma::mat X, arma::vec w);
RcppExport SEXP _glmhdfe_demean_criterion(SEXP XSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(demean_criterion(X, w));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_glmhdfe_irls_gaussian_identity", (DL_FUNC) &_glmhdfe_irls_gaussian_identity, 4},
    {"_glmhdfe_fe_gaussian_identity_foc", (DL_FUNC) &_glmhdfe_fe_gaussian_identity_foc, 2},
    {"_glmhdfe_irls_gaussian_log", (DL_FUNC) &_glmhdfe_irls_gaussian_log, 4},
    {"_glmhdfe_fe_gaussian_log_foc", (DL_FUNC) &_glmhdfe_fe_gaussian_log_foc, 2},
    {"_glmhdfe_irls_poisson_log", (DL_FUNC) &_glmhdfe_irls_poisson_log, 4},
    {"_glmhdfe_fe_poisson_log_foc", (DL_FUNC) &_glmhdfe_fe_poisson_log_foc, 2},
    {"_glmhdfe_irls_gamma_log", (DL_FUNC) &_glmhdfe_irls_gamma_log, 4},
    {"_glmhdfe_fe_gamma_log_foc", (DL_FUNC) &_glmhdfe_fe_gamma_log_foc, 2},
    {"_glmhdfe_irls_inverse_gaussian_log", (DL_FUNC) &_glmhdfe_irls_inverse_gaussian_log, 4},
    {"_glmhdfe_fe_inverse_gaussian_log_foc", (DL_FUNC) &_glmhdfe_fe_inverse_gaussian_log_foc, 2},
    {"_glmhdfe_ols", (DL_FUNC) &_glmhdfe_ols, 3},
    {"_glmhdfe_wls", (DL_FUNC) &_glmhdfe_wls, 4},
    {"_glmhdfe_wdemean", (DL_FUNC) &_glmhdfe_wdemean, 2},
    {"_glmhdfe_demean", (DL_FUNC) &_glmhdfe_demean, 1},
    {"_glmhdfe_wdemean_list", (DL_FUNC) &_glmhdfe_wdemean_list, 2},
    {"_glmhdfe_demean_list", (DL_FUNC) &_glmhdfe_demean_list, 1},
    {"_glmhdfe_gradient", (DL_FUNC) &_glmhdfe_gradient, 2},
    {"_glmhdfe_vcov_estimate", (DL_FUNC) &_glmhdfe_vcov_estimate, 3},
    {"_glmhdfe_demean_criterion", (DL_FUNC) &_glmhdfe_demean_criterion, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_glmhdfe(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
