fill_sigmas = function(fit){
  mapply(FUN = function(means, sigma){
    rep(sigma, length(means))
  }, fit$means, fit$sigmas)
}


pmixnorm_vec_loop = function(xs, means_list, sigmas_list, weights_list){
  nit = length(means_list)
  mapply(function(means, sigmas, weights){weights*pnorm(q = xs, mean = means, sd = sigmas)}, means_list, sigmas_list, weights_list)
}

get_CDF_full_BNPdensity = function(fit, xs = seq(-5,5, length.out = 100)){

  pmixnorm_vec_loop(xs = xs, means_list = fit$means, sigmas_list = fit$sigmas, weights_list = fit$weights)

}

get_CDF_semi_BNPdensity = function(fit, xs = seq(-5,5, length.out = 100)){

  fit$sigmas_filled = fill_sigmas(fit)

  pmixnorm_vec_loop(xs = xs, means_list = fit$means, sigmas_list = fit$sigmas_filled, weights_list = fit$weights)

}

#' Plot the empirical and fitted CDF for non censored data.
#'
#' @param data The non censored dataset.
#' @param fit The result of the fit, obtained through the function MixNRMI1 or MixNRMI2.
#' @return Plot of the empirical and fitted CDF for non censored data.
#' @examples
#' set.seed(150520)
#' data(acidity)
#' x <- enzyme
#' out <- MixNRMI1(enzyme, extras = TRUE)
#' plotCDF_noncensored(out)
plotCDF_noncensored = function(data, fit){

  data_range = max(data) - min(data)
  grid = seq(min(data)-0.1*data_range, max(data)+0.1*data_range, length.out = 100)




}