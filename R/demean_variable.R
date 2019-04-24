#' Pseudo-demean variables
#'
#' @param variable Character vector of variables that should be demeaned
#' @param fixed_effects Character vector of fixed effects whose influence should be removed from variables
#' @param weights Character string of weights for demeaning
#' @param data Data.table with relevant data
#' @param demean_iterations Maximum number of demeaning iterations
#' @param demean_tolerance Demeaning tolerance
#' @param parallel_cores Number of cores to be used in demeaning. Currently fixed at 1.
#' @param suffix Suffix for demeaned variables.
#' @param trace Show some information during estimation
#' @param verbose Show more information during estimation for the impatient
#'
#' @export demean_variable
demean_variable <- function(variable, fixed_effects, data,
                           weights = NULL,
                           demean_iterations = 100,
                           demean_tolerance = 1e-8,
                           parallel_cores = 1,
                           suffix = "_dm",
                           trace = F,
                           verbose = F) {
  # define local variables
  w <- NULL
  w_sqrt <- NULL

  # keep only necessary variables
  data <- data[, c(variable, fixed_effects, weights), with = F]

  # make sure all variables are numeric
  data[, (variable) := lapply(mget(variable), as.numeric)]

  # pseudo-demean all variables
  if (is.null(weights)) {
    for (i in fixed_effects) data[, (variable) := demean_list(mget(variable)), by = c(i)]
  } else {
    criterion <- demean_criterion(as.matrix(data[, variable, with = F]), data[[weights]])
    m <- 1
    repeat {
      criterion_old <- criterion
      pretty_message(verbose | trace, str_c("Demeaning transformed variables, iteration: ", m), time = T, linebreak = T)
      for (i in fixed_effects) data[, (variable) := wdemean_list(mget(variable), get(weights)), by = c(i)]
      criterion <- demean_criterion(as.matrix(data[, variable, with = F]), data[[weights]])
      criterion_change <- crossprod(criterion - criterion_old) / crossprod(criterion_old)
      pretty_message(verbose, str_c("change: ", formatC(criterion_change, format = "e", digits = 2)), task = T, linebreak = T)
      if (m > demean_iterations) break
      if (criterion_change < demean_tolerance) break
      m <- m + 1
    }
    for (j in variable) data[, (j) := sqrt(get(weights)) * get(j)]
  }

  # return demeaned variables
  setnames(data, variable, str_c(variable, suffix))
  return(data[, str_c(variable, suffix), with = F])
}
