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
  if (!is.null(x[["clustering"]])) cat("Note: Standard errors clustered on ", str_c(x[["clustering"]], collapse = " + "), ".", sep = "")
  cat("\n\nn = ", x[["nobs"]],", deviance = ", round(x[["deviance"]], digits = digits), "\n", sep = "")
}

#' Summarizing models of class \code{glmhdfe}
#'
#' @description Summary statistics for objects of class \code{glmhdfe}.
#' @param object Object of class \code{glmhdfe}.
#' @param se_type Standard errors. Select from "hessian" (default), "robust", "outer_product_gradient", or "clustered"
#' @param ...	further arguments
#'
#' @return
#' Returns an object of class \code{summary.glmhdfe} which is a list of summary statistics of
#' \code{object}.
#' @export
summary.glmhdfe = function(object, se_type = "hessian", ...) {
  if(is.null(object[["coefficients"]])) stop(call. = F, "No point coefficients estimated, only fixed effects.")
  if(se_type == "clustered" && is.null(object[["vcov"]][["vcov_clustered"]])) stop(call. = F, "No cluster defined.")
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
