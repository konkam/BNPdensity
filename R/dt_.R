#' Non-standard student-t density
#'
#' Computes the density.
#'
#' For internal use
#'
#' ## The function is currently defined as
#' function(x, df, mean, sd) {
#'   dt((x - mean) / sd, df, ncp = 0) / sd
#' }
dt_ <-
  function(x, df, mean, sd) {
    dt((x - mean) / sd, df, ncp = 0) / sd
  }
