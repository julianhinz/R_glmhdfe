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
#' @import crayon
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
                    accelerate = TRUE,
                    accelerate_iterations = 10,
                    accelerate_aux_vector = TRUE,
                    compute_vcov = TRUE,
                    demean_variables = FALSE,
                    demean_iterations = 100,
                    demean_tolerance = tolerance,
                    include_fe = TRUE,
                    include_data = FALSE,
                    include_data_vcov = FALSE,
                    skip_checks = NULL,
                    trace = FALSE,
                    verbose = FALSE){

  pretty_message(verbose | trace, "Preparations", time = TRUE, linebreak = TRUE)

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
  if (length(formula)[[1]] == 0) stop(call. = FALSE, "No dependent variable specified")
  if (length(formula)[[1]] > 1) stop(call. = FALSE, "More than one dependent variable specified")
  if (nrow(data) == 0 | length(data) == 0) stop(call. = FALSE, "No data specified")

  # hard-coded family-link combination
  family_link_hc <- c("gaussian_identity", "gaussian_log", "poisson_log", "Gamma_log", "inverse.gaussian_log")
  family_link <- paste0(family[["family"]], "_", family[["link"]])
  if (!family_link %in% family_link_hc) stop(call. = FALSE, "This distribution not yet supported.")

  # get lhs, rhs variables and fixed effects
  all_vars <- attr(terms(formula, lhs = 0), "term.labels")
  lhs_var <- deparse(formula[[2]])
  rhs_var <- attr(terms(formula, lhs = 0, rhs = 1), "term.labels")
  rhs_var_length <- length(rhs_var)
  if (rhs_var_length == 0) rhs_var <- NULL
  rhs_fe <- NULL
  rhs_fe_est <- NULL
  if (length(formula)[[2]] == 1) rhs_fe_length <- 0
  if (length(formula)[[2]] > 1) {
    rhs_fe <- attr(terms(formula, lhs = 0, rhs = 2), "term.labels")
    rhs_fe_length <- length(rhs_fe)
    if (rhs_fe_length > 0) {
      rhs_fe_est <- paste0(rhs_fe, "_est")
    } else {
      rhs_fe <- NULL
      rhs_fe_est <- NULL
    }
  }
  if (is.null(rhs_var) & is.null(rhs_fe)) stop(call. = FALSE, "No dependent variable specified")

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
  pretty_message(verbose, "set up model frame", task = TRUE)
  if (!is.data.table(data)) data = as.data.table(data) # not by reference, i.e. setDT, to keep original data untouched
  vars_needed <- unique(all.vars(formula))
  missing_vars <- setdiff(vars_needed, colnames(data))
  if (length(missing_vars) > 0) {
    cat("\n")
    stop(call. = FALSE, paste0("Variables ", paste0(missing_vars, collapse = ", "), " not in data"))
  }
  model_frame_cols <- lhs_var
  if (!is.null(rhs_var)) model_frame_cols <- append(model_frame_cols, rhs_var)
  if (!is.null(rhs_fe)) model_frame_cols <- append(model_frame_cols, rhs_fe)
  if (!is.null(rhs_cluster_se)) model_frame_cols <- append(model_frame_cols, rhs_cluster_se)
  data <- data[, lapply(model_frame_cols, function(expr_text) eval(parse_text(expr_text)))] # allows transformations of the data in the formula
  setnames(data, model_frame_cols)
  if (length(rhs_var[rhs_var == "(Intercept)"]) == 1) data[, "(Intercept)" := 1]
  setnames(data, lhs_var, "lhs") # big speed increase because not calling get(lhs_var)
  pretty_message(verbose, checkmark = TRUE)

  # only complete cases
  pretty_message(verbose, "subset to complete cases", task = TRUE)
  if(!"complete_cases" %in% skip_checks) data <- data[complete.cases(data)]
  pretty_message(verbose, checkmark = TRUE)

  # keep only those linear combinations of groups where there is variation
  nobs_data <- nrow(data)
  if (!"separation" %in% skip_checks & rhs_fe_length > 0) {
    pretty_message(verbose, "remove combinations of groups without variation", task = TRUE)
    data[, keep := TRUE]
    repeat {
      nobs <- nrow(data)
      # for (i in list_combinations(rhs_fe)) data[, keep := keep & not_min_max(lhs), by = mget(i)]
      for (i in rhs_fe) data[, keep := keep & not_min_max(lhs), by = get(i)]
      data <- data[keep == TRUE]
      if (nrow(data) == nobs) break
    }
    data[, keep := NULL]
    pretty_message(verbose, paste0(": dropped ", nobs_data - nobs, " and kept ", nobs," observations"), linebreak = TRUE)
    if (nobs == 0) {
      cat("\n")
      stop(call. = FALSE, "All observations dropped.")
    }
  } else {
    nobs <- nobs_data
  }

  # detect and remove multicollinear variables using QR decomposition
  if (!is.null(rhs_var) && !"multicollinearity" %in% skip_checks) {
    pretty_message(verbose, "remove multicollinear variables", task = TRUE)

    # Build design matrix, demeaning by FEs if present
    if (!is.null(rhs_fe)) {
      dm_result <- demean_variable(variable = rhs_var,
                                   fixed_effects = rhs_fe,
                                   data = data,
                                   demean_iterations = demean_iterations,
                                   demean_tolerance = demean_tolerance,
                                   verbose = FALSE)
      design_mat <- as.matrix(dm_result)
    } else {
      design_mat <- as.matrix(data[, rhs_var, with = FALSE])
    }

    # QR decomposition to detect multicollinearity
    qr_result <- qr(design_mat)
    qr_rank <- qr_result$rank
    if (qr_rank < length(rhs_var)) {
      keep_cols <- qr_result$pivot[seq_len(qr_rank)]
      dropped_vars <- rhs_var[-keep_cols]
      if (verbose | trace) cat(paste0(": ", paste0(dropped_vars, collapse = ", "), " dropped"))
      rhs_var <- rhs_var[keep_cols]
      rhs_var_length <- length(rhs_var)
    }
    pretty_message(verbose, checkmark = TRUE)
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
  pretty_message(verbose, "prepare for IRLS", task = TRUE)
  if (!is.null(rhs_var)) {
    vars <- as.matrix(data[, rhs_var, with = FALSE])
  } else {
    vars <- as.matrix(rep_len(0, nrow(data)))
  }
  pretty_message(verbose, checkmark = TRUE)

  # initialize estimated fixed effects
  if (rhs_fe_length > 0) {
    pretty_message(verbose, "initialize estimated fixed effects", task = TRUE)
    data[, paste0(rhs_fe, "_est") := 0.01]
    pretty_message(verbose, checkmark = TRUE)
  }

  # initialize offset, eta, beta and deviance
  pretty_message(verbose, "initialize offset, beta and deviance", task = TRUE)
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
  pretty_message(verbose, checkmark = TRUE)

  # function dispatch tables (avoid repeated string comparisons in loop)
  fe_foc_fn <- switch(family_link,
    gaussian_identity = fe_gaussian_identity_foc,
    gaussian_log = fe_gaussian_log_foc,
    poisson_log = fe_poisson_log_foc,
    Gamma_log = fe_gamma_log_foc,
    inverse.gaussian_log = fe_inverse_gaussian_log_foc
  )
  irls_fn <- switch(family_link,
    gaussian_identity = irls_gaussian_identity,
    gaussian_log = irls_gaussian_log,
    poisson_log = irls_poisson_log,
    Gamma_log = irls_gamma_log,
    inverse.gaussian_log = irls_inverse_gaussian_log
  )

  # iteration
  pretty_message(verbose | trace, "Begin iterations", time = TRUE, linebreak = TRUE)
  change <- 1
  change_fe <- 1
  accelerate_i <- 0

  # pre-allocate beta_list
  beta_list <- matrix(NA_real_, nrow = max_iterations + 1, ncol = rhs_var_length)
  beta_list[1, ] <- beta
  beta_list_row <- 1

  change_list <- as.matrix(t(c(as.integer(Sys.time()),0)))
  h <- 1

  if ((family_link %in% family_link_hc)) {

    # routine for hard-coded family-link combinations with nice fixed effects expression
    repeat {

      pretty_message(verbose | trace, paste0("Iteration ", h), time = TRUE, linebreak = TRUE)

      # compare criterion
      dev_old <- dev

      # update fixed effects
      if (!is.null(rhs_fe)) {

        pretty_message(verbose, "update fixed effects", task = TRUE)

        # if beta provided, update fixed effects more than once to "catch up"
        if (h == 1) { #&& sum(beta) > 0){
          Z <- fe_iterations
        } else {
          Z <- 1
        }

        for (z in seq(Z)) for (i in seq(rhs_fe_length)){

          # set theta
          data[, theta := eval(parse_text(c("eta", rhs_fe_est[-i])))]

          data[, rhs_fe_est[i] := fe_foc_fn(lhs, theta), by = eval(rhs_fe[i])]

        }
        pretty_message(verbose, checkmark = TRUE)
      }

      # recompute offset
      if (!is.null(rhs_fe)) data[, offset := eval(offset_eval)]

      # compute new beta
      if (!is.null(rhs_var)) {

        if (!is.null(rhs_fe) & accelerate_aux_vector) vars <- as.matrix(data[, c(rhs_var, rhs_fe_est), with = FALSE])

        beta <- irls_fn(A = vars, b = data[["lhs"]], eta = data[["eta"]], offset = data[["offset"]])

        if (accelerate_aux_vector) beta[["beta"]] <- t(beta[["beta"]][,seq(rhs_var_length)])

        beta_list_row <- beta_list_row + 1
        beta_list[beta_list_row, ] <- beta[["beta"]]

        pretty_message(verbose,
                       paste0("- updated beta: ", paste0(round(beta[["beta"]], 4), collapse = ", "),
                             ", change: ", formatC(sqrt(sum((beta_list[beta_list_row - 1, ] - beta_list[beta_list_row, ])^2)), format = "e", digits = 2)),
                       linebreak = TRUE)

        data[, eta := beta[["eta"]]]
        dev <- beta[["dev"]]

        # accelerate
        if (accelerate) {
          if (h > rhs_var_length) accelerate_i <- accelerate_i + 1
          pretty_message(verbose, paste0("acceleration earliest in ", max(accelerate_iterations - accelerate_i + 1, 1), " iterations"), task = TRUE, linebreak = TRUE)
          if (accelerate_i > accelerate_iterations) {
            if (sum(abs(beta_list[beta_list_row - 2, ] - beta_list[beta_list_row - 1, ]) > abs(beta_list[beta_list_row - 1, ] - beta_list[beta_list_row, ])) == length(beta[["beta"]])) {
              beta[["beta"]] <- accelerate_beta(beta_list[seq_len(beta_list_row), , drop = FALSE], t(beta[["beta"]]), beta_list_row, rhs_var)
              eta_accel <- as.vector(vars[, seq(rhs_var_length), drop = FALSE] %*% as.numeric(beta[["beta"]]))
              data[, eta := eta_accel]
              pretty_message(verbose, paste0("accelerated to: ", paste0(round(beta[["beta"]], 4), collapse = ", ")), task = TRUE, linebreak = TRUE)
              accelerate_i <- 0
            }
          }
        }

      } else {
        dev <- sum(family[["dev.resids"]](data[["lhs"]], family[["linkinv"]](data[["offset"]]), 1))
      }

      # continue?
      change <- (dev_old - dev) / dev
      change_list <- rbind(change_list, as.matrix(t(c(as.integer(Sys.time()), change))))
      if (h < 3) {
        pretty_message(verbose, paste0("current deviance: ", round(dev, digits = 2),
                                      ", change: ", formatC(change, format = "e", digits = 2)),
                       task = TRUE, linebreak = TRUE)
      } else {
        pretty_message(verbose, paste0("current deviance: ", round(dev, digits = 2),
                                      ", change: ", formatC(change, format = "e", digits = 2)),
                       task = TRUE, linebreak = TRUE)
      }

      if (h == max_iterations) break
      if (abs(change) <= tolerance) break

      h <- h + 1
    }

    # convergence warning
    if (h == max_iterations) {
      warning("glmhdfe did not converge after ", max_iterations, " iterations")
    }

    # trim beta_list to actual size
    beta_list <- beta_list[seq_len(beta_list_row), , drop = FALSE]

    # recompute mu from final beta and offset
    if (!is.null(rhs_var)) {
      eta_final <- as.vector(as.matrix(data[, rhs_var, with = FALSE]) %*% as.numeric(beta[["beta"]]))
      data[, mu := family[["linkinv"]](eta_final + offset)]
    } else {
      data[, mu := family[["linkinv"]](offset)]
    }

    # demean variables
    if (!is.null(rhs_var) && (compute_vcov | demean_variables)) {

      # score
      data[, score := lhs - mu]
      if (call_object[["link"]] == "log") data[, score := score / mu]

      # weights
      if (call_object[["family_link"]] == "gaussian_identity") {
        data[, w := 1]
        w = "w"
      }
      if (call_object[["family_link"]] == "gaussian_log") {
        data[, w := mu^2]
        w = "w"
      }
      if (call_object[["family_link"]] == "poisson_log") {
        data[, w := mu]
        w = "w"
      }
      if (call_object[["family_link"]] == "Gamma_log") {
        data[, w := 1]
        w = "w"
      }
      if (call_object[["family_link"]] == "inverse.gaussian_log") {
        data[, w := 1/mu]
        w = "w"
      }

      data <- cbind(data, demean_variable(variable = c("lhs", rhs_var, "score"),
                                          fixed_effects = rhs_fe,
                                          weights = w,
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
    stop(call. = FALSE, "Routine for generic family-link combination not yet implemented")
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
  if (include_data == TRUE) object[["data"]] <- data # data matrix
  if (!is.null(rhs_fe) && include_fe == TRUE) object[["fe"]] <- fe_object # estimated fixed effects
  if (!is.null(rhs_var) && compute_vcov == TRUE) object[["vcov"]] <- vcov_object # variance covariance matrix object, may be robust or clustered
  if (!is.null(rhs_var)) object[["score"]] <- beta[["score"]] # score
  if (!is.null(rhs_var)) object[["sse"]] <- beta[["sse"]] # sum of squared errors
  if (!is.null(rhs_var)) object[["deviance"]] <- beta[["dev"]] # deviance
  if (verbose == TRUE) object[["beta_list"]] <- beta_list
  if (verbose == TRUE) object[["change_list"]] <- change_list
  object[["family"]] <- family # store family for S3 methods

  return(structure(object, class = "glmhdfe"))
}
