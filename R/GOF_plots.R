#' @import stats
#' @import ggplot2

pmix_vec_loop = function(xs, locations_list, scales_list, weights_list, distr.k){
  additive_mix_vec_loop(xs, locations_list, scales_list, weights_list, distr.k, pk)
}

dmix_vec_loop = function(xs, locations_list, scales_list, weights_list, distr.k){
  additive_mix_vec_loop(xs, locations_list, scales_list, weights_list, distr.k, dk)
}

additive_mix_vec_loop = function(xs, locations_list, scales_list, weights_list, distr.k, distfun){
  res = 0.0*xs
  nit = length(locations_list)
  for (it in 1:nit){
    for (cmp in seq_along(locations_list[[it]]))
      res = res + weights_list[[it]][cmp] * distfun(xs, distr = distr.k, mu = locations_list[[it]][cmp], sigma = scales_list[[it]][cmp])
  }
  return(res/nit)
}


get_CDF_full_BNPdensity = function(fit, xs = seq(-5,5, length.out = 100)){
  pmix_vec_loop(xs = xs, locations_list = fit$means, scales_list = fit$sigmas, weights_list = fit$weights, distr.k = fit$distr.k)
}
get_PDF_full_BNPdensity = function(fit, xs = seq(-5,5, length.out = 100)){
  dmix_vec_loop(xs = xs, locations_list = fit$means, scales_list = fit$sigmas, weights_list = fit$weights, distr.k = fit$distr.k)
}

get_CDF_semi_BNPdensity = function(fit, xs = seq(-5,5, length.out = 100)){
  fit$sigmas_filled = fill_sigmas(fit)
  pmix_vec_loop(xs = xs, locations_list = fit$means, scales_list = fit$sigmas_filled, weights_list = fit$weights, distr.k = fit$distr.k)
}
get_PDF_semi_BNPdensity = function(fit, xs = seq(-5,5, length.out = 100)){
  fit$sigmas_filled = fill_sigmas(fit)
  dmix_vec_loop(xs = xs, locations_list = fit$means, scales_list = fit$sigmas_filled, weights_list = fit$weights, distr.k = fit$distr.k)
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

  grid = grid_from_data(data)

  if(is_semiparametric(fit)){
    cdf = get_CDF_semi_BNPdensity(fit = fit, xs = grid)
  }
  else{
    cdf = get_CDF_full_BNPdensity(fit = fit, xs = grid)
  }
  ggplot2::ggplot(data = data.frame(data = grid, CDF = cdf), aes(x = data, y = CDF)) +
    geom_line(colour= 'red') +
    theme_classic() +
    stat_ecdf(data = data.frame(data), aes(y = NULL), geom = "step") +
    xlab("Data")
}

#' Plot the density and a histogram for non censored data.
#'
#' @param data The non censored dataset.
#' @param fit The result of the fit, obtained through the function MixNRMI1 or MixNRMI2.
#' @return Plot of the ensity and a histogram for non censored data.
#' @examples
#' set.seed(150520)
#' data(acidity)
#' out <- MixNRMI1(acidity, extras = TRUE, Nit = 100)
#' plotPDF_noncensored(out)
plotPDF_noncensored = function(fit){

  data = fit$data

  grid = grid_from_data(data)

  if(is_semiparametric(fit)){
    pdf = get_PDF_semi_BNPdensity(fit = fit, xs = grid)
  }
  else{
    pdf = get_PDF_full_BNPdensity(fit = fit, xs = grid)
  }
  ggplot2::ggplot(data = data.frame(data = grid, PDF = pdf), aes(x = data, y = PDF)) +
    geom_histogram(data = data.frame(data), aes(y = ..density..)) +
    geom_line(colour= 'red') +
    theme_classic() +
    xlab("Data")

}

#' Plot the density and a histogram for non censored data.
#'
#' @param data The non censored dataset.
#' @param fit The result of the fit, obtained through the function MixNRMI1 or MixNRMI2.
#' @return Plot of the ensity and a histogram for non censored data.
#' @examples
#' set.seed(150520)
#' data(acidity)
#' out <- MixNRMI1(acidity, extras = TRUE, Nit = 100)
#' plotPDF_noncensored(out)
pp_plot_noncensored = function(fit){

  data = fit$data

  # grid = grid_from_data(data)

  if(is_semiparametric(fit)){
    cdf = get_CDF_semi_BNPdensity(fit = fit, xs = data)
  }
  else{
    cdf = get_CDF_full_BNPdensity(fit = fit, xs = data)
  }
  ggplot2::ggplot(data = data.frame(x = cdf, y = ecdf(data)(data)), aes(x = x, y = y)) +
    geom_point() +
    geom_abline(slope = 1, intercept = 0, colour= 'red') +
    theme_classic() +
    xlab("Theoretical cumulative distribution") +
    ylab("Empirical cumulative distribution")

}

#' Plot Goodness of fits graphical checks for non censored data
#'
#' @param fit
#'
#' @return A density plot with histogram, a cumulative density plot with the empirical cumulative distribution, and a percentile-percentile plot.
#' @export
#'
#' @examples
#' set.seed(150520)
#' data(acidity)
#' out <- MixNRMI1(acidity, extras = TRUE, Nit = 100)
#' plotGOF_noncensored(out)
plotGOF_noncensored = function(fit){

  CDFplot = plotCDF_noncensored(fit)
  PDFplot = plotPDF_noncensored(fit)
  pplot = pp_plot_noncensored(fit)
  gridExtra::grid.arrange(PDFplot, CDFplot, pplot)
}