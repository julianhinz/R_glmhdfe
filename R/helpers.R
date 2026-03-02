#' Parse character formula
#'
#' @param t Text to be concatenated and parsed
#' @param s Separator
#' @param c Collapse
#'
parse_text <- function (t, s="", c = "+") parse(text=paste0(t, collapse = c))

#' Check if even number
#'
#' @param x Integer to check if even
is.even <- function(x) as.numeric(x %% 2 == 0)

#' Check if minimum not equal to maximum
#'
#' @param x Numeric vector
not_min_max <- function(x) min(x) != max(x)

#' List of combinations
#'
#' @importFrom utils combn
#' @param elements Vector of elements
list_combinations = function (elements) {
  return(unlist(
    lapply(
      X = lapply(
        X = seq(elements), FUN = combn, x = elements, simplify = FALSE),
      FUN = as.list
    ),
    recursive = FALSE,
    use.names = FALSE)
  )
}

#' Pretty messages for the impatient
#'
#' @param verbose Be verbose and display messages?
#' @param description Message to be displayed
#' @param time Timestamp?
#' @param task Make task, i.e. dash before description
#' @param linebreak Break line after message?
#' @param color Colorful messages
#' @param checkmark Display checkmark at the end of message?
#'
#' @import crayon
#'
pretty_message <- function(verbose = FALSE, description = NULL, time = FALSE, task = FALSE, linebreak = FALSE, color = black, checkmark = FALSE) {
  if (verbose) {
    if (time) message(color(paste0(Sys.time(), ": ")), appendLF = FALSE)
    if (task) message(color("- "), appendLF = FALSE)
    if (!is.null(description)) message(color(description), appendLF = FALSE)
    if (checkmark) message(green(" \u2713\n"), appendLF = FALSE)
    if (linebreak) message("\n", appendLF = FALSE)
  }
}

#' Predict convergence time
#'
#' @param change Vector of relative change in variable of interest
#' @param time Vector of time, relative or absolute
#' @param tolerance Tolerance for which convergence is assumed to be achieved
#' @importFrom utils tail
predict_convergence_time <- function(time, change, tolerance) {
  y <- log(tail(change,3))
  time <- tail(time,3)
  start <- time[1]
  time <- (time - min(time))
  X <- matrix(c(rep(1,3), time), ncol = 2)
  b <- solve(t(X) %*% X) %*% t(X) %*% y
  x <- start + (log(tolerance) - b[1]) / b[2]
  return(x)
}
