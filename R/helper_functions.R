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
#' @param fit The result of the fit, obtained through the function MixNRMI1 or MixNRMI2.
#' @return TRUE if the fit is a semiparametric model
#' @examples
#' set.seed(150520)
#' data(acidity)
#' x <- enzyme
#' out <- MixNRMI1(enzyme, extras = TRUE)
#' is_semiparametric(out)
fill_sigmas = function(semiparametric_fit){
  mapply(FUN = function(means, sigma){
    rep(sigma, length(means))
  }, semiparametric_fit$means, semiparametric_fit$S)
}