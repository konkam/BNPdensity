pmixnorm_vec_loop = function(xs, means_list, sigmas_list, weights_list){
  res = 0.0*xs
  nit = length(means_list)
  for (it in 1:nit){
    for (cmp in seq_along(means_list[[it]]))
    res = res + weights_list[[it]][cmp] * pnorm(q = xs, mean = means_list[[it]][cmp], sd = sigmas_list[[it]][cmp])
  }
  return(res/nit)
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
  p = ggplot2::ggplot(data = data.frame(data = grid, CDF = cdf), aes(x = data, y = CDF)) +
    geom_line(colour= 'red') +
    theme_classic() +
    stat_ecdf(data = data.frame(data), aes(y = NULL), geom = "step")
}