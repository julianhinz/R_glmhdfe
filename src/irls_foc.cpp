#include <RcppArmadillo.h>
#include <omp.h>
//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace R;

/////
// Gaussian identity
/////
// [[Rcpp::export]]
Rcpp::List irls_gaussian_identity(arma::mat A, arma::vec b, arma::vec eta, arma::vec offset)
{
  int nobs = A.n_rows;
  int nvars = A.n_cols;
  arma::vec W(nobs);
  arma::vec z(nobs);
  arma::vec beta;
  arma::vec mu(nobs);
  arma::vec e(nobs);
  arma::vec r(nobs);
  arma::vec sse(nobs);
  double dev;

  // mu
  mu = offset+eta;

  // z
  z = eta+(b-mu);

  //beta
  beta = arma::solve(A.t()*A, A.t()*z);

  // deviance residuals
  r = pow(b - mu, 2);
  dev = sum(r);

  //sse
  e = (b - A*beta);
  sse = e.t()*e;

  return Rcpp::List::create(
    Rcpp::Named("beta") = beta.t(),
    Rcpp::Named("eta") = A*beta,
    Rcpp::Named("mu") = offset+eta,
    Rcpp::Named("nobs") = nobs,
    Rcpp::Named("nvars") = nvars,
    Rcpp::Named("dev") = dev
  );
}

// [[Rcpp::export]]
double fe_gaussian_identity_foc(arma::vec lhs, arma::vec theta)
{
  return mean(lhs - theta);
}


/////
// Gaussian log
/////
// [[Rcpp::export]]
Rcpp::List irls_gaussian_log(arma::mat A, arma::vec b, arma::vec eta, arma::vec offset)
{
  int nobs = A.n_rows;
  int nvars = A.n_cols;
  arma::vec W(nobs);
  arma::vec z(nobs);
  arma::vec beta;
  arma::vec mu(nobs);
  arma::vec e(nobs);
  arma::vec r(nobs);
  arma::vec sse(nobs);
  double dev;

  // mu
  mu = exp(offset+eta);

  // weights
  W = mu % mu;

  // z
  z = eta+(b-mu)/mu;

  //beta
  beta = arma::solve(A.t()*(W % A.each_col()), A.t()*(W % z));

  // deviance residuals
  r = pow(b - mu, 2);
  dev = sum(r);

  // sse
  e = (b - A*beta);
  sse = e.t()*e;

  return Rcpp::List::create(
    Rcpp::Named("beta") = beta.t(),
    Rcpp::Named("eta") = A*beta,
    Rcpp::Named("mu") = exp(offset+eta),
    Rcpp::Named("nobs") = nobs,
    Rcpp::Named("nvars") = nvars,
    Rcpp::Named("dev") = dev
  );
}

// [[Rcpp::export]]
double fe_gaussian_log_foc(arma::vec lhs, arma::vec theta)
{
  return log(sum(lhs % exp(theta)) / sum(pow(exp(theta),2)));
}


/////
// Poisson log
/////
// [[Rcpp::export]]
Rcpp::List irls_poisson_log(arma::mat A, arma::vec b, arma::vec eta, arma::vec offset)
{
  int nobs = A.n_rows;
  int nvars = A.n_cols;
  arma::vec W(nobs);
  arma::vec z(nobs);
  arma::vec beta;
  arma::vec mu(nobs);
  arma::vec e(nobs);
  arma::vec r(nobs);
  arma::vec sse(nobs);
  double dev;

  // mu
  mu = exp(offset+eta);

  // weights
  W = mu;

  // z
  z = eta+(b-mu)/mu;

  //beta
  beta = arma::solve(A.t()*(W % A.each_col()), A.t()*(W % z));

  // deviance
  r = 2 * (b % log(b/mu) - (b - mu));
  arma::uvec zeros = find(b == 0); // find zeros
  arma::uvec not_zeros = find(b != 0); // find zeros
  arma::vec replace = 2*mu.elem(zeros); // replace them
  dev = sum(r.elem(not_zeros)) + sum(replace);

  // sse
  e = (b - A*beta);
  sse = e.t()*e;

  return Rcpp::List::create(
    Rcpp::Named("beta") = beta.t(),
    Rcpp::Named("eta") = A*beta,
    Rcpp::Named("mu") = exp(offset+eta),
    Rcpp::Named("nobs") = nobs,
    Rcpp::Named("nvars") = nvars,
    Rcpp::Named("dev") = dev
  );
}

// [[Rcpp::export]]
double fe_poisson_log_foc(arma::vec lhs, arma::vec theta)
{
  return log(sum(lhs) / sum(exp(theta)));
}


/////
// Gamma log
/////
// [[Rcpp::export]]
Rcpp::List irls_gamma_log(arma::mat A, arma::vec b, arma::vec eta, arma::vec offset)
{
  int nobs = A.n_rows;
  int nvars = A.n_cols;
  arma::vec W(nobs);
  arma::vec z(nobs);
  arma::vec beta;
  arma::vec mu(nobs);
  arma::vec e(nobs);
  arma::vec r(nobs);
  arma::vec sse(nobs);
  double dev;

  // mu
  mu = exp(offset+eta);

  // z
  z = eta+(b-mu)/mu;

  //beta
  beta = arma::solve(A.t()*A, A.t()*z);

  // deviance residuals
  r = 2 * ((b - mu) / mu - log(b/mu));
  arma::uvec zeros = find(b == 0); // find zeros
  arma::uvec not_zeros = find(b != 0); // find zeros
  arma::vec replace = 2*mu.elem(zeros); // replace them
  dev = sum(r.elem(not_zeros)) + sum(replace);

  //sse
  e = (b - A*beta);
  sse = e.t()*e;

  return Rcpp::List::create(
    Rcpp::Named("beta") = beta.t(),
    Rcpp::Named("eta") = A*beta,
    Rcpp::Named("mu") = exp(offset+eta),
    Rcpp::Named("nobs") = nobs,
    Rcpp::Named("nvars") = nvars,
    Rcpp::Named("dev") = dev
  );
}

// [[Rcpp::export]]
double fe_gamma_log_foc(arma::vec lhs, arma::vec theta)
{
  return log(mean(lhs / exp(theta)));
}

/////
// Inverse Gaussian log
/////
// [[Rcpp::export]]
Rcpp::List irls_inverse_gaussian_log(arma::mat A, arma::vec b, arma::vec eta, arma::vec offset)
{
  int nobs = A.n_rows;
  int nvars = A.n_cols;
  arma::vec W(nobs);
  arma::vec z(nobs);
  arma::vec beta;
  arma::vec mu(nobs);
  arma::vec e(nobs);
  arma::vec r(nobs);
  arma::vec sse(nobs);
  double dev;

  // mu
  mu = exp(offset+eta);

  // weights
  W = pow(mu,-1);

  // z
  z = eta+(b-mu)/mu;

  //beta
  beta = arma::solve(A.t()*(W % A.each_col()), A.t()*(W % z));

  // deviance residuals
  r = pow(b - mu, 2);
  dev = sum(r);

  // sse
  e = (b - A*beta);
  sse = e.t()*e;

  return Rcpp::List::create(
    Rcpp::Named("beta") = beta.t(),
    Rcpp::Named("eta") = A*beta,
    Rcpp::Named("mu") = exp(offset+eta),
    Rcpp::Named("nobs") = nobs,
    Rcpp::Named("nvars") = nvars,
    Rcpp::Named("dev") = dev
  );
}

// [[Rcpp::export]]
double fe_inverse_gaussian_log_foc(arma::vec lhs, arma::vec theta)
{
  return log(sum(lhs % pow(exp(theta),-2)) / sum(pow(exp(theta),-1)));
}
