#' Extract coefficient estimates
#'
#' @param object Object of class \code{glmhdfe}
#' @param ...	further arguments
#'
#' @export
coef.glmhdfe = function (object, ...) object[["coefficients"]]

#' Extract fitted values
#'
#' @param object Object of class \code{glmhdfe}
#' @param ...	further arguments
#'
#' @export
fitted.glmhdfe = function (object, ...) object[["fitted.values"]]

#' Extract response residuals
#'
#' @param object Object of class \code{glmhdfe}
#' @param ...	further arguments
#'
#' @export
residuals.glmhdfe = function (object, ...) object[["residuals"]]

#' Predict method for glmhdfe objects
#'
#' @param object Object of class \code{glmhdfe}
#' @param newdata Optional new data frame for prediction. If NULL, returns fitted values.
#' @param type Type of prediction: "link" for linear predictor, "response" for response scale
#' @param ...	further arguments
#' @importFrom stats predict
#'
#' @export
predict.glmhdfe = function(object, newdata = NULL, type = c("response", "link"), ...) {
  type <- match.arg(type)

  if (is.null(newdata)) {
    if (type == "response") {
      return(object[["fitted.values"]])
    } else {
      family <- object[["family"]]
      if (is.null(family)) family <- get_family(object)
      return(family[["linkfun"]](object[["fitted.values"]]))
    }
  }

  # prediction on new data (linear predictor from X*beta only, no FEs)
  if (!is.null(object[["call"]][["rhs_fe"]])) {
    warning("predict.glmhdfe uses only X*beta and ignores fixed effects", call. = FALSE)
  }

  rhs_var <- object[["call"]][["rhs_var"]]
  if (is.null(rhs_var)) stop("No coefficients to predict with", call. = FALSE)

  if (!is.data.table(newdata)) newdata <- as.data.table(newdata)
  X_dt <- newdata[, lapply(rhs_var, function(expr_text) eval(parse_text(expr_text)))]
  setnames(X_dt, rhs_var)
  if (length(rhs_var[rhs_var == "(Intercept)"]) == 1) X_dt[, "(Intercept)" := 1]
  X <- as.matrix(X_dt)
  eta <- c(X %*% object[["coefficients"]])

  if (type == "link") {
    return(eta)
  } else {
    family <- object[["family"]]
    if (is.null(family)) family <- get_family(object)
    return(family[["linkinv"]](eta))
  }
}

#' Extract variance-covariance matrix
#'
#' @param object Object of class \code{glmhdfe}
#' @param type Type of vcov: "hessian" (default), "robust", "outer_product_gradient", or "clustered".
#'   When clusters were specified and type is not given, defaults to "clustered".
#' @param ...	further arguments
#' @importFrom stats vcov
#'
#' @export
vcov.glmhdfe = function(object, type = NULL, ...) {
  if (is.null(object[["vcov"]])) stop("No variance-covariance matrix available", call. = FALSE)

  if (is.null(type)) {
    if (!is.null(object[["vcov"]][["vcov_clustered"]])) {
      type <- "clustered"
    } else {
      type <- "hessian"
    }
  }

  vcov_name <- paste0("vcov_", type)
  if (type == "outer_product_gradient") vcov_name <- "vcov_opg"
  result <- object[["vcov"]][[vcov_name]]
  if (is.null(result)) stop(paste0("vcov type '", type, "' not available"), call. = FALSE)
  return(result)
}

#' Number of observations
#'
#' @param object Object of class \code{glmhdfe}
#' @param ...	further arguments
#' @importFrom stats nobs
#'
#' @export
nobs.glmhdfe = function(object, ...) {
  object[["info"]][["nobs"]]
}

#' Log-likelihood for glmhdfe objects
#'
#' @param object Object of class \code{glmhdfe}
#' @param ...	further arguments
#'
#' @importFrom stats dpois dgamma dnorm logLik
#'
#' @export
logLik.glmhdfe = function(object, ...) {
  y <- object[["response"]]
  mu <- object[["fitted.values"]]
  family_link <- object[["call"]][["family_link"]]
  n <- object[["info"]][["nobs"]]
  p <- length(object[["coefficients"]]) + sum(unlist(object[["info"]][["nfe"]]))

  ll <- switch(family_link,
    gaussian_identity = ,
    gaussian_log = {
      sigma2 <- sum((y - mu)^2) / n
      sum(dnorm(y, mean = mu, sd = sqrt(sigma2), log = TRUE))
    },
    poisson_log = {
      sum(dpois(y, lambda = mu, log = TRUE))
    },
    Gamma_log = {
      phi <- sum((y/mu - 1)^2) / n
      shape <- 1 / phi
      sum(dgamma(y, shape = shape, rate = shape / mu, log = TRUE))
    },
    inverse.gaussian_log = {
      lambda <- n / sum((y - mu)^2 / (y * mu^2))
      sum(0.5 * (log(lambda) - log(2 * pi * y^3)) - lambda * (y - mu)^2 / (2 * y * mu^2))
    }
  )

  structure(ll, df = p, nobs = n, class = "logLik")
}

#' Confidence intervals for glmhdfe coefficients
#'
#' @param object Object of class \code{glmhdfe}
#' @param parm Specification of which parameters to give confidence intervals for
#' @param level Confidence level (default 0.95)
#' @param type Type of vcov to use: "hessian", "robust", "outer_product_gradient", or "clustered"
#' @param ...	further arguments
#'
#' @importFrom stats qnorm confint
#'
#' @export
confint.glmhdfe = function(object, parm = NULL, level = 0.95, type = NULL, ...) {
  cf <- coef(object)
  if (is.null(cf)) stop("No coefficients available", call. = FALSE)

  V <- vcov.glmhdfe(object, type = type)
  se <- sqrt(diag(V))

  a <- (1 - level) / 2
  z <- qnorm(1 - a)

  ci <- cbind(cf - z * se, cf + z * se)
  colnames(ci) <- paste0(formatC(c(a, 1 - a) * 100, format = "f", digits = 1), " %")

  if (!is.null(parm)) {
    if (is.numeric(parm)) {
      ci <- ci[parm, , drop = FALSE]
    } else {
      ci <- ci[parm, , drop = FALSE]
    }
  }

  return(ci)
}

#' Print \code{glmhdfe} object
#'
#' @param x Object of class \code{glmhdfe}
#' @param digits Integer indicating the number of decimal places (defaults to 4)
#' @param ...	further arguments
#'
#' @export
print.glmhdfe = function (x, ..., digits = 4) {
  print(coef(x), digits = digits)
}

#' Print \code{summary.glmhdfe}
#'
#' @description \code{print.summary.glmhdfe} is a generic function which displays summary statistics from objects
#' returned by \code{summary.glmhdfe}.
#' @param x Object of class \code{summary.glmhdfe}.
#' @param digits Integer indicating the number of decimal places (defaults to 4)
#' @param ...	further arguments
#'
#' @importFrom stats printCoefmat
#'
#' @export
print.summary.glmhdfe <- function(x, ..., digits = 4) {
  cat("Model information:\n")
  cat("- family: ", x[["family"]], " with ", x[["link"]], " link\n", sep = "")
  cat("- estimated equation: ")
  print(x[["formula"]])
  cat("\n")
  cat("Model fit:\n")
  printCoefmat(x[["coefficients"]], P.values = TRUE, has.Pvalue = TRUE, digits = digits)
  if (!is.null(x[["clustering"]])) cat("Note: Standard errors clustered on ", paste0(x[["clustering"]], collapse = " + "), ".", sep = "")
  cat("\n\nn = ", x[["nobs"]],", deviance = ", round(x[["deviance"]], digits = digits), "\n", sep = "")
}

#' Summarizing models of class \code{glmhdfe}
#'
#' @description Summary statistics for objects of class \code{glmhdfe}.
#' @param object Object of class \code{glmhdfe}.
#' @param se_type Standard errors. Select from "hessian" (default), "robust", "outer_product_gradient", or "clustered".
#'   When clusters were specified in the formula, defaults to "clustered".
#' @param ...	further arguments
#'
#' @return
#' Returns an object of class \code{summary.glmhdfe} which is a list of summary statistics of
#' \code{object}.
#' @export
summary.glmhdfe = function(object, se_type = NULL, ...) {
  if(is.null(object[["coefficients"]])) stop(call. = FALSE, "No point coefficients estimated, only fixed effects.")
  if (is.null(object[["vcov"]])) stop(call. = FALSE, "No variance-covariance matrix available. Refit with compute_vcov = TRUE.")

  # default se_type: "clustered" when clusters were specified, otherwise "hessian"
  if (is.null(se_type)) {
    if (!is.null(object[["call"]][["rhs_cluster_se"]])) {
      se_type <- "clustered"
    } else {
      se_type <- "hessian"
    }
  }

  if(se_type == "clustered" && is.null(object[["vcov"]][["vcov_clustered"]])) stop(call. = FALSE, "No cluster defined.")
  beta <- object[["coefficients"]]
  if (se_type == "hessian") se <- sqrt(diag(object[["vcov"]][["vcov_hessian"]]))
  if (se_type == "robust") se <- sqrt(diag(object[["vcov"]][["vcov_robust"]]))
  if (se_type == "outer_product_gradient") se <- sqrt(diag(object[["vcov"]][["vcov_opg"]]))
  if (se_type == "clustered") se <- sqrt(diag(object[["vcov"]][["vcov_clustered"]]))
  z <- beta / se
  p <- 2.0 * pnorm(-abs(z))
  coefficients <- cbind(beta, se, z, p)
  rownames(coefficients) <- names(beta)
  colnames(coefficients) <- c("Estimate", "Std. error", "z value", "Pr(>|z|)")

  output <- list(coefficients = coefficients,
                 deviance = object[["deviance"]],
                 nobs = object[["info"]][["nobs_data"]],
                 formula = object[["call"]][["formula"]],
                 family = object[["call"]][["family"]],
                 link = object[["call"]][["link"]])
  if(se_type == "clustered") output[["clustering"]] = object[["call"]][["rhs_cluster_se"]]

  structure(output, class = "summary.glmhdfe")
}

#' Helper to reconstruct family object from call
#' @param object glmhdfe object
#' @keywords internal
get_family <- function(object) {
  family_name <- object[["call"]][["family"]]
  link_name <- object[["call"]][["link"]]
  do.call(family_name, list(link = link_name))
}
