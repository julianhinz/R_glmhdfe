#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace R;

/////
// weighted demeaning
/////
// [[Rcpp::export]]
arma::vec wdemean(arma::vec var, arma::vec weight)
{
  return (var - sum(weight % var) / sum(weight));
}

// [[Rcpp::export]]
arma::vec demean(arma::vec var)
{
  return (var - sum(var) / var.n_elem);
}

// [[Rcpp::export]]
Rcpp::List wdemean_list(Rcpp::List var,
                        Rcpp::NumericVector weight)
{
  int n = var.size();
  NumericVector v;
  for (int i = 0; i < n; i++)
  {
    v = var[i];
    var[i] = v - sum(weight * v) / sum(weight);
  }
  return var;
}

// [[Rcpp::export]]
Rcpp::List demean_list(Rcpp::List var)
{
  int n = var.size();
  NumericVector v;
  for (int i = 0; i < n; i++)
  {
    v = var[i];
    var[i] = v - mean(v);
  }
  return var;
}

/////
// return gradient for score and rhs_pv
/////
// [[Rcpp::export]]
arma::mat gradient(arma::mat X, arma::vec score)
{
  return(X.each_col() % score);
}

/////
// estimate different vcov matrices
/////
// [[Rcpp::export]]
Rcpp::List vcov_estimate(arma::mat X, arma::mat grad, arma::mat cluster_meat)
{
  arma::mat inv_hessian = arma::pinv(arma::trans(X) * X);
  arma::mat op_grad = arma::trans(grad) * grad;

  return Rcpp::List::create(
    Rcpp::Named("vcov_hessian") = inv_hessian,
    Rcpp::Named("vcov_opg") = arma::pinv(op_grad),
    Rcpp::Named("vcov_robust") = inv_hessian * op_grad * inv_hessian,
    Rcpp::Named("vcov_clustered") = inv_hessian * cluster_meat * inv_hessian
  );
}

/////
// demean criterion
/////
// [[Rcpp::export]]
arma::mat demean_criterion(arma::mat X, arma::vec w)
{
  X = X.each_col() % sqrt(w);
  X = arma::trans(X) * X;
  return X.diag();
}
