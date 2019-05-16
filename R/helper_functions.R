#' Tests if a fit is a semi parametric or nonparametric model.
#'
#' @param fit The result of the fit, obtained through the function MixNRMI1 or MixNRMI2.
#' @return TRUE if the fit is a semiparametric model
#' @examples
#' set.seed(150520)
#' data(acidity)
#' x <- enzyme
#' out <- MixNRMI1(enzyme, extras = TRUE)
#' is_semiparametric(out)
is_semiparametric = function(fit){
  return(!is.null(fit$S))
}

#' Repeat the common scale parameter of a semiparametric model to match the dimension of the location parameters.
#'
#' @param semiparametric_fit The result of the fit, obtained through the function MixNRMI1.
#' @return an adequate list of vectors of sigmas
fill_sigmas = function(semiparametric_fit){
  mapply(FUN = function(means, sigma){
    rep(sigma, length(means))
  }, semiparametric_fit$means, semiparametric_fit$S)
}

#' Create a plotting grid from non-censored data.
#'
#' @param data Input data from which to compute the grid.
#' @param npoints Number of points on the grid.
#' @return a vector containing the plotting grid
grid_from_data = function(data, npoints = 100){
  data_range = max(data) - min(data)
  return(seq(min(data) - 0.1*data_range, max(data) + 0.1*data_range, length.out = 100))
}