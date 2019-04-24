#' Generalized linear model with high-dimensional fixed effects
#'
#' This function allows for arbitrary sets of fixed effects in a generalized linear model.
#'
#' @param formula Formula that describes lefthandside variable, righthandside variables of interest, sets of fixed effects and error clustering, e.g. `y ~ x | fe1 + fe2 | cluster1 + cluster2`
#' @param data `data.table` or `data.frame` with data used in the regression
#' @param family Estimator used, currently limited to gaussian(link = "identity"), gaussian(link = "log"), poisson(link = "log"), Gamma(link = "log") and inverse.gaussian(link = "log")
#' @param beta Initial beta vector, defaults to 0
#' @param tolerance Iteration stops when tolerance is reached (numeric)
#' @param max_iterations Iteration stops when maximum number of iteration is reached (integer)
#' @param fe_iterations Iterations for fixed effects to "catch up" when initial beta vector is provided
#' @param accelerate Acceleration algorithm for finding betas
#' @param accelerate_iterations Number of iterations before starting acceleration algorithm
#' @param accelerate_aux_vector Include fixed effects vectors in IRLS? Usually increases convergence speed
#' @param compute_vcov Compute variance-covariance matrix? Can be computed ex-post when data from estimation provided.
#' @param demean_variables Demean variables to be used in estimation of variance-covariance matrix?
#' @param demean_iterations Iterations for demeaning algorithm
#' @param demean_tolerance Demeaning tolerance
#' @param include_fe Return fixed effects?
#' @param include_data Return data used in estimation?
#' @param include_data_vcov Return data used in variance-covariance matrix estimation?
#' @param skip_checks Skip checks before starting procedure? Character vector can include "separation", "complete_cases", "multicollinearity".
#' @param trace Show some information during estimation
#' @param verbose Show more information during estimation for the impatient
#'
#' @import data.table
#' @import Formula
#' @import stringr
#' @import crayon
#' @importFrom lfe felm
#' @importFrom stats gaussian coef complete.cases dnorm model.frame model.matrix na.omit pnorm terms
#'
#' @export glmhdfe

glmhdfe <- function(formula,
                    data,
                    family = gaussian(),
                    beta = NULL,
                    tolerance = 1e-8,
                    max_iterations = 10000,
                    fe_iterations = 5,
                    accelerate = T,
                    accelerate_iterations = 10,
                    accelerate_aux_vector = T,
                    compute_vcov = T,
                    demean_variables = F,
                    demean_iterations = 100,
                    demean_tolerance = tolerance,
                    include_fe = T,
                    include_data = F,
                    include_data_vcov = F,
                    skip_checks = NULL,
                    trace = F,
                    verbose = F){

  pretty_message(verbose | trace, "Preparations", time = T, linebreak = T)

  # define local variables
  . <- NULL
  keep <- NULL
  lhs <- NULL
  eta <- NULL
  offset <- NULL
  theta <- NULL
  mu <- NULL
  fe <- NULL
  score <- NULL
  w <- NULL

  # sanity checks
  formula <- Formula(formula)
  if (length(formula)[[1]] == 0) stop(call. = F, "No dependent variable specified")
  if (length(formula)[[1]] > 1) stop(call. = F, "More than one dependent variable specified")
  if (nrow(data) == 0 | length(data) == 0) stop(call. = F, "No data specified")

  # hard-coded family-link combination
  family_link_hc <- c("gaussian_identity", "gaussian_log", "poisson_log", "Gamma_log", "inverse.gaussian_log")
  family_link <- str_c(family[["family"]], "_", family[["link"]])
  if (!family_link %in% family_link_hc) stop(call. = F, "This distribution not yet supported.")

  # get lhs, rhs variables and fixed effects
  all_vars <- attr(terms(formula, lhs = 0), "term.labels")
  lhs_var <- deparse(formula[[2]])
  rhs_var <- attr(terms(formula, lhs = 0, rhs = 1), "term.labels")
  rhs_var_length <- length(rhs_var)
  if (rhs_var_length == 0) rhs_var <- NULL
  if (length(formula)[[2]] == 1) rhs_fe_length <- 0
  if (length(formula)[[2]] > 1) {
    rhs_fe <- attr(terms(formula, lhs = 0, rhs = 2), "term.labels")
    rhs_fe_length <- length(rhs_fe)
    if (rhs_fe_length > 0) {
      rhs_fe_est <- str_c(rhs_fe, "_est")
    } else {
      rhs_fe <- NULL
      rhs_fe_est <- NULL
    }
  }
  if (is.null(rhs_var) & is.null(rhs_fe)) stop(call. = F, "No dependent variable specified")

  # clustered standard errors
  rhs_cluster_se <- NULL
  if (length(formula)[[2]] == 3 && (attr(formula, "rhs")[[3]] != "0")) {
    rhs_cluster_se <- attr(terms(formula, lhs = 0, rhs = 3), "term.labels")
  }

  # save call in call_object
  call_object <- structure(list(), class = "glmhdfe_call")
  call_object[["family"]] <- family[["family"]]
  call_object[["link"]] <- family[["link"]]
  call_object[["family_link"]] <- family_link
  call_object[["formula"]] <- formula
  call_object[["lhs_var"]] <- lhs_var
  call_object[["rhs_var"]] <- rhs_var
  call_object[["rhs_fe"]] <- rhs_fe
  call_object[["rhs_cluster_se"]] <- rhs_cluster_se

  # skip checks
  if (is.null(skip_checks)) skip_checks <- ""
  if ("all" %in% skip_checks) skip_checks <- c("separation", "complete_cases", "multicollinearity")

  # set up data table with lhs, rhs variables, fixed effects and clusters
  pretty_message(verbose, "set up model frame", task = T)
  if (!is.data.table(data)) data = as.data.table(data) # not by reference, i.e. setDT, to keep original data untouched
  if (length(colnames(data)[all_vars %in% colnames(data)]) < length(all_vars)) {
    cat("\n")
    stop(call. = F, str_c("Variables ", str_c(all_vars[!all_vars %in% colnames(data)], collapse = ", "), " not in data"))
  }
  model_frame_cols <- lhs_var
  if (!is.null(rhs_var)) model_frame_cols <- append(model_frame_cols, rhs_var)
  if (!is.null(rhs_fe)) model_frame_cols <- append(model_frame_cols, rhs_fe)
  if (!is.null(rhs_cluster_se)) model_frame_cols <- append(model_frame_cols, rhs_cluster_se)
  data <- data[, lapply(model_frame_cols, function(x) eval(parse_text(x)))] # allows transformations of the data in the formula
  setnames(data, model_frame_cols)
  if (length(rhs_var[rhs_var == "(Intercept)"]) == 1) data[, "(Intercept)" := 1]
  setnames(data, lhs_var, "lhs") # big speed increase because not calling get(lhs_var)
  pretty_message(verbose, checkmark = T)

  # only complete cases
  pretty_message(verbose, "subset to complete cases", task = T)
  if(!"complete_cases" %in% skip_checks) data <- data[complete.cases(data)]
  pretty_message(verbose, checkmark = T)

  # keep only those linear combinations of groups where there is variation
  nobs_data <- nrow(data)
  if (!"separation" %in% skip_checks & rhs_fe_length > 0) {
    pretty_message(verbose, "remove combinations of groups without variation", task = T)
    data[, keep := T]
    repeat {
      nobs <- nrow(data)
      # for (i in list_combinations(rhs_fe)) data[, keep := keep & not_min_max(lhs), by = mget(i)]
      for (i in rhs_fe) data[, keep := keep & not_min_max(lhs), by = get(i)]
      data <- data[keep == T]
      if (nrow(data) == nobs) break
    }
    data[, keep := NULL]
    pretty_message(verbose, str_c(": dropped ", nobs_data - nobs, " and kept ", nobs," observations"), linebreak = T)
    if (nobs == 0) {
      cat("\n")
      stop(call. = F, "All observations dropped.")
    }
  } else {
    nobs <- nobs_data
  }

  # detect and remove multicollinear, get starting values for beta
  if (!is.null(rhs_var) && !"multicollinearity" %in% skip_checks) {
    pretty_message(verbose, "remove multicollinear variables", task = T)
    if (!is.null(rhs_fe)) {
      help_formula <- as.Formula(str_c("lhs ~ ", str_c("`", rhs_var, "`", collapse = " + "), " | ", str_c(rhs_fe, collapse = " + ")))
    } else {
      help_formula <- as.Formula(str_c("lhs ~ 0 + ", str_c("`", rhs_var, "`", collapse = " + ")))
    }
    suppressWarnings(help_reg <- felm(help_formula, data))
    if ((verbose | trace) &&
        (sum(is.nan(coef(help_reg))) > 0)) cat(str_c(": ", str_c(names(coef(help_reg))[is.nan(coef(help_reg))], collapse =", "), " dropped"))
    rhs_var <- rhs_var[!is.nan(coef(help_reg))]
    pretty_message(verbose, checkmark = T)
  }

  # save data information in info_object
  info_object <- structure(list(), class = "glmhdfe_info")
  info_object[["nobs"]] <- nobs # number of observations used
  info_object[["nobs_data"]] <- nobs_data # number of observations of data
  info_object[["nvars"]] <- rhs_var_length # number of variables of interest
  if (!is.null(rhs_fe)) {
    info_object[["nfe"]] <- as.list(data[,lapply(.SD, function(x) length(unique(x))), .SDcols = rhs_fe]) # number of levels of fixed effects
  } else {
    info_object[["nfe"]] <- 0
  }
  info_object[["p"]] <- info_object[["nvars"]] + sum(unlist(info_object[["nfe"]])) # number of all variables
  info_object[["df"]] <- info_object[["nobs"]] - info_object[["p"]] - 1 # degrees of freedom
  rm(nobs_data, nobs)

  # set key for efficient allocation in memory
  if (!is.null(rhs_fe)) setkeyv(data, rhs_fe)

  # prepare for irls
  pretty_message(verbose, "prepare for IRLS", task = T)
  if (!is.null(rhs_var)) {
    vars <- as.matrix(data[, rhs_var, with = F])
  } else {
    vars <- as.matrix(rep_len(0, nrow(data)))
  }
  pretty_message(verbose, checkmark = T)

  # initialize estimated fixed effects
  if (rhs_fe_length > 0) {
    pretty_message(verbose, "initialize estimated fixed effects", task = T)
    data[, str_c(rhs_fe, "_est") := 0.01]
    pretty_message(verbose, checkmark = T)
  }

  # initialize offset, eta, beta and deviance
  pretty_message(verbose, "initialize offset, beta and deviance", task = T)
  if (is.null(beta)) {
    beta <- t(rep(0.01, rhs_var_length))
  } else {
    beta <- t(beta)
  }
  if (!is.null(rhs_var)) names(beta) <- rhs_var
  if (!is.null(rhs_fe)) offset_eval <- parse_text(rhs_fe_est)
  if (!is.null(rhs_var)) {
    data[, eta := tcrossprod(vars,beta)]
  } else {
    data[, eta := 0.01]
  }
  if (!is.null(rhs_fe)) {
    data[, offset := eval(offset_eval)]
  } else {
    data[, offset := 0]
  }

  dev <- sum(family[["dev.resids"]](data[["lhs"]], family[["linkinv"]](data[["offset"]]), 1))
  pretty_message(verbose, checkmark = T)

  # iteration
  pretty_message(verbose | trace, "Begin iterations", time = T, linebreak = T)
  change <- 1
  change_fe <- 1
  accelerate_i <- 0
  beta_list <- as.matrix(beta)
  change_list <- as.matrix(t(c(Sys.time(),0)))
  h <- 1

  if ((family_link %in% family_link_hc)) {

    # routine for hard-coded family-link combinations with nice fixed effects expression
    repeat {

      pretty_message(verbose | trace, str_c("Iteration ", h), time = T, linebreak = T)

      # compare criterion
      dev_old <- dev

      # update fixed effects
      if (!is.null(rhs_fe)) {

        pretty_message(verbose, "update fixed effects", task = T)

        # if beta provided, update fixed effects more than once to "catch up"
        if (h == 1) { #&& sum(beta) > 0){
          Z <- fe_iterations
        } else {
          Z <- 1
        }

        for (z in seq(Z)) for (i in seq(rhs_fe_length)){

          # set theta
          data[, theta := eval(parse_text(c("eta", rhs_fe_est[-i])))]

          if (family_link == "gaussian_identity") data[, rhs_fe_est[i] := fe_gaussian_identity_foc(lhs,theta), by = eval(rhs_fe[i])]
          if (family_link == "gaussian_log") data[, rhs_fe_est[i] := fe_gaussian_log_foc(lhs,theta), by = eval(rhs_fe[i])]
          if (family_link == "poisson_log") data[, rhs_fe_est[i] := fe_poisson_log_foc(lhs,theta), by = eval(rhs_fe[i])]
          if (family_link == "Gamma_log") data[, rhs_fe_est[i] := fe_gamma_log_foc(lhs,theta), by = eval(rhs_fe[i])]
          if (family_link == "inverse.gaussian_log") data[, rhs_fe_est[i] := fe_inverse_gaussian_log_foc(lhs,theta), by = eval(rhs_fe[i])]

        }
        pretty_message(verbose, checkmark = T)
      }

      # recompute offset
      if (!is.null(rhs_fe)) data[, offset := eval(offset_eval)]

      # compute new beta
      if (!is.null(rhs_var)) {

        if (!is.null(rhs_fe) & accelerate_aux_vector) vars <- as.matrix(data[, c(rhs_var, rhs_fe_est), with = F])

        if (family_link == "gaussian_identity") beta <- irls_gaussian_identity(A = vars, b = data[["lhs"]], eta = data[["eta"]], offset = data[["offset"]])
        if (family_link == "gaussian_log") beta <- irls_gaussian_log(A = vars, b = data[["lhs"]], eta = data[["eta"]], offset = data[["offset"]])
        if (family_link == "poisson_log") beta <- irls_poisson_log(A = vars, b = data[["lhs"]], eta = data[["eta"]], offset = data[["offset"]])
        if (family_link == "Gamma_log") beta <- irls_gamma_log(A = vars, b = data[["lhs"]], eta = data[["eta"]], offset = data[["offset"]])
        if (family_link == "inverse.gaussian_log") beta <- irls_inverse_gaussian_log(A = vars, b = data[["lhs"]], eta = data[["eta"]], offset = data[["offset"]])

        if (accelerate_aux_vector) beta[["beta"]] <- t(beta[["beta"]][,seq(rhs_var_length)])
        pretty_message(verbose,
                       str_c("- updated beta: ", str_c(round(beta[["beta"]], 4), collapse = ", "),
                             ", change: ", formatC(sqrt(sum((beta_list[h - 1, ] - beta_list[h, ])^2)), format = "e", digits = 2)),
                       linebreak = T)

        data[, eta := beta[["eta"]]]
        beta_list <- rbind(beta_list, as.matrix(beta[["beta"]]))
        dev <- beta[["dev"]]

        # accelerate
        if (accelerate) {
          if (h > rhs_var_length) accelerate_i <- accelerate_i + 1
          pretty_message(verbose, str_c("acceleration earliest in ", accelerate_iterations - accelerate_i + 1, " iterations"), task = T, linebreak = T)
          if (accelerate_i > accelerate_iterations) {
            if (sum(abs(beta_list[h - 2, ] - beta_list[h - 1, ]) > abs(beta_list[h - 1, ] - beta_list[h, ])) == length(beta[["beta"]])) {
              beta[["beta"]] <- accelerate_beta(beta_list, t(beta[["beta"]]), h, rhs_var)
              data[, eta := tcrossprod(vars[,seq(rhs_var_length)],t(beta[["beta"]]))]
              pretty_message(verbose, str_c("accelerated to: ", str_c(round(beta[["beta"]], 4), collapse = ", ")), task = T, linebreak = T)
              accelerate_i <- 0
            }
          }
        }

      } else {
        dev <- sum(family[["dev.resids"]](data[["lhs"]], family[["linkinv"]](data[["offset"]]), 1))
      }

      # continue?
      change <- (dev_old - dev) / dev
      change_list <- rbind(change_list, as.matrix(t(c(Sys.time(), change))))
      if (h < 3) {
        pretty_message(verbose, str_c("current deviance: ", round(dev, digits = 2),
                                      ", change: ", formatC(change, format = "e", digits = 2)),
                       task = T, linebreak = T)
      } else {
        pretty_message(verbose, str_c("current deviance: ", round(dev, digits = 2),
                                      ", change: ", formatC(change, format = "e", digits = 2),
                                      ", ETA of convergence: ", as.POSIXct(predict_convergence_time(change_list[,1],change_list[,2], tolerance), origin="1970-01-01")),
                       task = T, linebreak = T)
      }

      if (h == max_iterations) break
      if (abs(change) <= tolerance) break

      h <- h + 1
    }

    # extract mu
    if (!is.null(rhs_var)) {
      data[, mu := beta[["mu"]]]
    } else {
      data[, mu := offset]
    }

    # demean variables
    if (!is.null(rhs_var) && (compute_vcov | demean_variables)) {

      # score
      data[, score := lhs - mu]
      if (call_object[["link"]] == "log") data[, score := score / mu]

      # weights
      if (call_object[["family_link"]] %in% c("gaussian_identity", "Gamma_log")) data[, w := 1]
      if (call_object[["family_link"]] == "gaussian_log") data[, w := mu^2]
      if (call_object[["family_link"]] == "poisson_log") data[, w := mu]

      data <- cbind(data, demean_variable(variable = c("lhs", rhs_var, "score"),
                                          fixed_effects = rhs_fe,
                                          weights = "w",
                                          data = data,
                                          demean_iterations = demean_iterations,
                                          demean_tolerance = demean_tolerance,
                                          verbose = verbose))
    }

    # standard errors
    if (!is.null(rhs_var) && compute_vcov) {
      vcov_object <- compute_vcov(data = data,
                                  call = call_object,
                                  info = info_object,
                                  custom_cluster = NULL,
                                  demean_iterations = demean_iterations,
                                  demean_tolerance = demean_tolerance,
                                  include_data_vcov = include_data_vcov,
                                  verbose = verbose,
                                  trace = trace)
    }

    # retrieve fixed effects
    if (!is.null(rhs_fe) && include_fe) {
      fe_object <- list()
      for (i in seq(rhs_fe_length)){
        fe_object[[rhs_fe[i]]] <- unique(data[, list(variable = get(rhs_fe[i]), value = get(rhs_fe_est[i]))])
      }
      fe_object <- structure(fe_object, class = "glmhdfe_fe")
    }

  } else {

    # routine for generic family-link combination with IRLS and lmhdfe
    stop(call. = F, "Routine for generic family-link combination not yet implemented")
  }

  # format coefficients
  if (!is.null(rhs_var)) {
    beta[["beta"]] <- c(beta[["beta"]])
    names(beta[["beta"]]) <- rhs_var
  } else {
    beta <- list()
    beta[["beta"]] <- NULL
  }

  # rename lhs
  setnames(data, "lhs", lhs_var)

  # prepare output
  object <- list()
  object[["call"]] <- call_object # call object includes information on formula, family and link function
  object[["info"]] <- info_object # info_object includes information on number of variables, degrees of freedom, etc.
  object[["coefficients"]] <- beta[["beta"]] # estimated coefficients
  object[["response"]] <- data[[lhs_var]] # response vector
  object[["fitted.values"]] <- data[["mu"]] # fitted values
  object[["residuals"]] <- object[["response"]] - object[["fitted.values"]] # residuals
  if (include_data == T) object[["data"]] <- data # data matrix
  if (!is.null(rhs_fe) && include_fe == T) object[["fe"]] <- fe_object # estimated fixed effects
  if (!is.null(rhs_var) && compute_vcov == T) object[["vcov"]] <- vcov_object # variance covariance matrix object, may be robust or clustered
  if (!is.null(rhs_var)) object[["score"]] <- beta[["score"]] # score
  if (!is.null(rhs_var)) object[["sse"]] <- beta[["sse"]] # sum of squared errors
  if (!is.null(rhs_var)) object[["deviance"]] <- beta[["dev"]] # deviance
  if (verbose == T) object[["beta_list"]] <- beta_list
  if (verbose == T) object[["change_list"]] <- change_list

  return(structure(object, class = "glmhdfe"))
}
