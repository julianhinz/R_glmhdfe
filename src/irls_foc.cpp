#include <RcppArmadillo.h>
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
  arma::vec score;
  double dev;

  // mu
  mu = offset+eta;

  // z
  z = eta+(b-mu);

  //beta
  beta = arma::solve(A.t()*A, A.t()*z);

  // update mu with new beta
  mu = offset + A*beta;

  // deviance residuals (with updated mu)
  r = pow(b - mu, 2);
  dev = sum(r);

  //sse
  e = (b - A*beta);
  sse = e.t()*e;

  // score
  score = A.t()*(b - mu);

  return Rcpp::List::create(
    Rcpp::Named("beta") = beta.t(),
    Rcpp::Named("eta") = A*beta,
    Rcpp::Named("mu") = mu,
    Rcpp::Named("nobs") = nobs,
    Rcpp::Named("nvars") = nvars,
    Rcpp::Named("dev") = dev,
    Rcpp::Named("sse") = sse,
    Rcpp::Named("score") = score
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
  arma::vec score;
  double dev;

  // mu
  mu = exp(offset+eta);

  // weights
  W = mu % mu;

  // z
  z = eta+(b-mu)/mu;

  //beta
  beta = arma::solve(A.t()*(W % A.each_col()), A.t()*(W % z));

  // update mu with new beta
  mu = exp(offset + A*beta);

  // deviance residuals (with updated mu)
  r = pow(b - mu, 2);
  dev = sum(r);

  // sse
  e = (b - A*beta);
  sse = e.t()*e;

  // score
  score = A.t()*(W % ((b - mu) / mu));

  return Rcpp::List::create(
    Rcpp::Named("beta") = beta.t(),
    Rcpp::Named("eta") = A*beta,
    Rcpp::Named("mu") = mu,
    Rcpp::Named("nobs") = nobs,
    Rcpp::Named("nvars") = nvars,
    Rcpp::Named("dev") = dev,
    Rcpp::Named("sse") = sse,
    Rcpp::Named("score") = score
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
  arma::vec score;
  double dev;

  // mu
  mu = exp(offset+eta);

  // weights
  W = mu;

  // z
  z = eta+(b-mu)/mu;

  //beta
  beta = arma::solve(A.t()*(W % A.each_col()), A.t()*(W % z));

  // update mu with new beta
  mu = exp(offset + A*beta);

  // deviance (with updated mu)
  r = 2 * (b % log(b/mu) - (b - mu));
  arma::uvec zeros = find(b == 0); // find zeros
  arma::uvec not_zeros = find(b != 0); // find non-zeros
  arma::vec replace = 2*mu.elem(zeros); // replace them
  dev = sum(r.elem(not_zeros)) + sum(replace);

  // sse
  e = (b - A*beta);
  sse = e.t()*e;

  // score
  score = A.t()*(b - mu);

  return Rcpp::List::create(
    Rcpp::Named("beta") = beta.t(),
    Rcpp::Named("eta") = A*beta,
    Rcpp::Named("mu") = mu,
    Rcpp::Named("nobs") = nobs,
    Rcpp::Named("nvars") = nvars,
    Rcpp::Named("dev") = dev,
    Rcpp::Named("sse") = sse,
    Rcpp::Named("score") = score
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
  arma::vec score;
  double dev;

  // mu
  mu = exp(offset+eta);

  // z
  z = eta+(b-mu)/mu;

  //beta
  beta = arma::solve(A.t()*A, A.t()*z);

  // update mu with new beta
  mu = exp(offset + A*beta);

  // deviance residuals (with updated mu)
  r = 2 * ((b - mu) / mu - log(b/mu));
  arma::uvec zeros = find(b == 0); // find zeros
  arma::uvec not_zeros = find(b != 0); // find non-zeros
  arma::vec replace = 2*mu.elem(zeros); // replace them
  dev = sum(r.elem(not_zeros)) + sum(replace);

  //sse
  e = (b - A*beta);
  sse = e.t()*e;

  // score
  score = A.t()*((b - mu) / mu);

  return Rcpp::List::create(
    Rcpp::Named("beta") = beta.t(),
    Rcpp::Named("eta") = A*beta,
    Rcpp::Named("mu") = mu,
    Rcpp::Named("nobs") = nobs,
    Rcpp::Named("nvars") = nvars,
    Rcpp::Named("dev") = dev,
    Rcpp::Named("sse") = sse,
    Rcpp::Named("score") = score
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
  arma::vec score;
  double dev;

  // mu
  mu = exp(offset+eta);

  // weights
  W = pow(mu,-1);

  // z
  z = eta+(b-mu)/mu;

  //beta
  beta = arma::solve(A.t()*(W % A.each_col()), A.t()*(W % z));

  // update mu with new beta
  mu = exp(offset + A*beta);

  // deviance residuals: inverse Gaussian deviance = (b - mu)^2 / (b * mu^2)
  r = pow(b - mu, 2) / (b % pow(mu, 2));
  dev = sum(r);

  // sse
  e = (b - A*beta);
  sse = e.t()*e;

  // score
  score = A.t()*(pow(mu, -1) % ((b - mu) / mu));

  return Rcpp::List::create(
    Rcpp::Named("beta") = beta.t(),
    Rcpp::Named("eta") = A*beta,
    Rcpp::Named("mu") = mu,
    Rcpp::Named("nobs") = nobs,
    Rcpp::Named("nvars") = nvars,
    Rcpp::Named("dev") = dev,
    Rcpp::Named("sse") = sse,
    Rcpp::Named("score") = score
  );
}

// [[Rcpp::export]]
double fe_inverse_gaussian_log_foc(arma::vec lhs, arma::vec theta)
{
  return log(sum(lhs % pow(exp(theta),-2)) / sum(pow(exp(theta),-1)));
}


/////
// ordinary least squares
/////

// [[Rcpp::export]]
Rcpp::List ols(arma::mat X, arma::vec y, bool full)
{
  int nobs = X.n_rows;
  int nvars = X.n_cols;
  arma::vec e(nobs);
  arma::vec mse;
  arma::mat hessian;
  arma::vec score;
  arma::vec beta;
  Rcpp::List output;

  // hessian
  hessian = X.t()*X;

  // score
  score = X.t()*y;

  // beta
  beta = arma::solve(hessian, score);

  if (!full) {
    output = output.create(
      Rcpp::Named("beta") = beta
    );
  }

  if (full) {
    // residuals and mse
    e = (y - X*beta);
    mse = e.t()*e / nobs;

    output = output.create(
      Rcpp::Named("beta") = beta.t(),
      Rcpp::Named("mu") = X*beta,
      Rcpp::Named("residuals") = e,
      Rcpp::Named("hessian") = hessian,
      Rcpp::Named("score") = score,
      Rcpp::Named("nobs") = nobs,
      Rcpp::Named("nvars") = nvars,
      Rcpp::Named("mse") = mse
    );
  }

  return output;
}

/////
// weighted least squares
/////

// [[Rcpp::export]]
Rcpp::List wls(arma::mat X, arma::vec y, arma::vec w, bool full)
{
  int nobs = X.n_rows;
  int nvars = X.n_cols;
  arma::vec e(nobs);
  arma::vec mse;
  arma::mat hessian;
  arma::vec score;
  arma::vec beta;
  Rcpp::List output;

  // hessian
  hessian = X.t()*(w % X.each_col());

  // score
  score = X.t()*(w % y);

  // beta
  beta = arma::solve(hessian, score);

  if (!full) {
    output = output.create(
      Rcpp::Named("beta") = beta
    );
  }

  if (full) {
    // residuals and mse
    e = (y - X*beta);
    mse = e.t()*e / nobs;

    output = output.create(
      Rcpp::Named("beta") = beta.t(),
      Rcpp::Named("mu") = X*beta,
      Rcpp::Named("residuals") = e,
      Rcpp::Named("hessian") = hessian,
      Rcpp::Named("score") = score,
      Rcpp::Named("nobs") = nobs,
      Rcpp::Named("nvars") = nvars,
      Rcpp::Named("mse") = mse
    );
  }

  return output;
}
