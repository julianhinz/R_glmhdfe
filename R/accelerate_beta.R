#' Acceleration algorithm for betas
#'
#' @param beta_list List of previous estimated coefficients
#' @param beta Vector of current estimated coefficients
#' @param h Acceleration iteration
#' @param rhs_var Righthandside variables

accelerate_beta = function(beta_list, beta, h, rhs_var) {
  step1 = beta_list[(h - length(beta)):(h - 1), ] - beta_list[(h - length(beta) - 1):(h - 2), ] %>% as.matrix()
  step2 = beta_list[(h - length(beta) + 1):(h), ] - beta_list[(h - length(beta)):(h - 1), ] %>% as.matrix()
  lambda = t(step2) %*% solve(t(step1))
  beta_old = beta
  beta = beta + 0.5*solve(diag(nrow = length(beta)) - lambda) %*% step1[1, ]
  beta = c(t(beta))
  names(beta) = rhs_var
  return(beta)
}
