#' Compute variance covariance matrix for standard errors from given glmhdfe_data, glmhdfe_call and glmhdfe_info objects
#'
#' @param data glmhdfe_data object
#' @param call glmhdfe_call object
#' @param info glmhdfe_info object
#' @param custom_cluster Character vector with variable names, overwrites clustering from glmhdfe_call object
#' @param demean_iterations Iterations for demeaning algorithm
#' @param demean_tolerance Demeaning tolerance
#' @param include_data_vcov Return data used in variance-covariance matrix estimation?
#' @param trace Show some information during estimation
#' @param verbose Show more information during estimation for the impatient
#'
#' @import stringr
#' @import Formula
#' @importFrom purrr map reduce
#'
#' @export compute_vcov
compute_vcov <- function (data, call, info,
                          custom_cluster = NULL,
                          demean_iterations = 100,
                          demean_tolerance = 1e-8,
                          include_data_vcov = F,
                          trace = F,
                          verbose = F) {

  # define local variables
  score  <- NULL
  mu <- NULL
  w <- NULL

  pretty_message(verbose | trace, "Compute variance covariance matrix", time = T)

  if (!is.null(custom_cluster)) call[["rhs_cluster_se"]] <- custom_cluster

  pretty_message(verbose = (verbose | trace && (!is.null(call[["rhs_cluster_se"]]))),
                 str_c(" clustered on ", str_c(call[["rhs_cluster_se"]], collapse = " + ")))
  pretty_message(verbose | trace, linebreak = T)

  # prepare variables
  if (!"score_dm" %in% colnames(data)) {

    pretty_message(verbose, "prepare variables", task = T)

    # score
    data[, score := get(call[["lhs_var"]]) - mu]
    if (call[["link"]] == "log") data[, score := score / mu]

    # weights
    if (call[["family_link"]] %in% c("gaussian_identity", "Gamma_log")) data[, w := 1]
    if (call[["family_link"]] == "gaussian_log") data[, w := mu^2]
    if (call[["family_link"]] == "poisson_log") data[, w := mu]
    pretty_message(verbose, checkmark = T)

    data <- cbind(data, demean_variable(variable = c(call[["rhs_var"]], "score"),
                                        fixed_effects = call[["rhs_fe"]],
                                        weights = "w",
                                        data = data,
                                        demean_iterations = demean_iterations,
                                        demean_tolerance = demean_tolerance,
                                        trace = trace,
                                        verbose = verbose))
  }

  # construct gradient and meat for clustered standard errors
  if (!is.null(call[["rhs_cluster_se"]])) {
    vcov_grad <- cbind(data[, call[["rhs_cluster_se"]], with = F],
                       gradient(as.matrix(data[,str_c(call[["rhs_var"]], "_dm"), with = F]), data[["score_dm"]]))
    setnames(vcov_grad, c(call[["rhs_cluster_se"]], str_c(call[["rhs_var"]], "_dm")))

    pretty_message(verbose, "compute meat for clustering", task = T)

    vcov_cluster_meat <- list_combinations(call[["rhs_cluster_se"]]) %>%
      map(~vcov_grad[, lapply(.SD, sum), by = eval(.x), .SDcols = str_c(call[["rhs_var"]], "_dm")]) %>%
      map(~((-1)^(is.even(length(colnames(.x))-length(call[["rhs_var"]])))*
              crossprod(as.matrix(.x[,str_c(call[["rhs_var"]], "_dm"), with = F])))) %>%
      reduce(`+`)
    pretty_message(verbose, checkmark = T)

  } else {
    vcov_grad <- data.table(gradient(as.matrix(data[,str_c(call[["rhs_var"]], "_dm"), with = F]), data[["score_dm"]]))
    setnames(vcov_grad, str_c(call[["rhs_var"]], "_dm"))
    vcov_cluster_meat <- diag(1, nrow = length(call[["rhs_var"]]), ncol = length(call[["rhs_var"]]))
  }

  # construct formula and estimate
  pretty_message(verbose, "estimate and return object", task = T)
  vcov_object <- vcov_estimate(X = as.matrix(data[,str_c(call[["rhs_var"]], "_dm"), with = F]),
                               grad = as.matrix(vcov_grad[,str_c(call[["rhs_var"]], "_dm"), with = F]),
                               cluster_meat = vcov_cluster_meat)
  if (!is.null(call[["rhs_cluster_se"]])) {
    attr(vcov_object[["vcov_clustered"]], "clustering") <- call[["rhs_cluster_se"]]
  } else{
    vcov_object[["vcov_clustered"]] <- NULL
  }
  if (!is.null(custom_cluster)) attr(vcov_object[["vcov_clustered"]], "clustering") <- custom_cluster

  pretty_message(verbose, checkmark = T)

  if (include_data_vcov) {
    vcov_object[["data"]] <- data[, call[["rhs_fe"]], with = F]
    if (!is.null(call[["rhs_cluster_se"]])) vcov_object[["data"]] <- cbind(vcov_object[["data"]],
                                                                           data[, call[["rhs_cluster_se"]], with = F])
    vcov_object[["data"]] <- cbind(vcov_object[["data"]],
                                   data[, c("score", "score_dm",
                                            call[["rhs_var"]], str_c(call[["rhs_var"]], "_dm"),
                                            "w"), with = F])
    vcov_object[["bread"]] <- vcov_object[["vcov_hessian"]]
    vcov_object[["meat"]] <- vcov_cluster_meat
    vcov_object[["grad"]] <- vcov_grad
  }

  return(structure(vcov_object, class = "glmhdfe_vcov"))
}
