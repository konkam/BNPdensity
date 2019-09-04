#' Non-standard student-t densit
#' 
#' Computes the density.
#' 
#' For internal use
#' 
#' @keywords internal
#' @examples
#' 
#' ## The function is currently defined as
#' function (x, df, mean, sd) 
#' {
#'     dt((x - mean)/sd, df, ncp = 0)/sd
#'   }
#' 
dt_ <-
  function(x, df, mean, sd) {
    dt((x - mean) / sd, df, ncp = 0) / sd
  }
