#' @import stats
#' @import ggplot2

pmix_vec_loop = function(xs, locations_list, scales_list, weights_list, distr.k){
  res = 0.0*xs
  nit = length(locations_list)
  for (it in 1:nit){
    for (cmp in seq_along(locations_list[[it]]))
      res = res + weights_list[[it]][cmp] * pk(q = xs, distr = distr.k, mu = locations_list[[it]][cmp], sigma = scales_list[[it]][cmp])
  }
  return(res/nit)
}

get_CDF_full_BNPdensity = function(fit, xs = seq(-5,5, length.out = 100)){

  pmix_vec_loop(xs = xs, locations_list = fit$means, scales_list = fit$sigmas, weights_list = fit$weights, distr.k = fit$distr.k)

}

get_CDF_semi_BNPdensity = function(fit, xs = seq(-5,5, length.out = 100)){

  fit$sigmas_filled = fill_sigmas(fit)

  pmix_vec_loop(xs = xs, locations_list = fit$means, scales_list = fit$sigmas_filled, weights_list = fit$weights, distr.k = fit$distr.k)

}

#' Plot the empirical and fitted CDF for non censored data.
#'
#' @param data The non censored dataset.
#' @param fit The result of the fit, obtained through the function MixNRMI1 or MixNRMI2.
#' @return Plot of the empirical and fitted CDF for non censored data.
#' @examples
#' set.seed(150520)
#' data(acidity)
#' out <- MixNRMI1(acidity, extras = TRUE)
#' plotCDF_noncensored(out)
plotCDF_noncensored = function(fit){

  data = fit$data

  data_range = max(data) - min(data)
  grid = seq(min(data)-0.1*data_range, max(data)+0.1*data_range, length.out = 100)

  if(is_semiparametric(fit)){
    cdf = get_CDF_semi_BNPdensity(fit = fit, xs = grid)
  }
  else{
    cdf = get_CDF_full_BNPdensity(fit = fit, xs = grid)
  }
  p = ggplot2::ggplot(data = data.frame(data = grid, CDF = cdf), ggplot2::aes(x = data, y = CDF)) +
    geom_line(colour= 'red') +
    theme_classic() +
    stat_ecdf(data = data.frame(data), aes(y = NULL), geom = "step")
}