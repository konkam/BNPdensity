pmix_vec_loop <-
  function(xs,
             locations_list,
             scales_list,
             weights_list,
             distr.k) {
    additive_mix_vec_loop(xs, locations_list, scales_list, weights_list, distr.k, pk)
  }

dmix_vec_loop <-
  function(xs,
             locations_list,
             scales_list,
             weights_list,
             distr.k) {
    additive_mix_vec_loop(xs, locations_list, scales_list, weights_list, distr.k, dk)
  }

additive_mix_vec_loop <-
  function(xs,
             locations_list,
             scales_list,
             weights_list,
             distr.k,
             distfun) {
    res <- 0.0 * xs
    nit <- length(locations_list)
    for (it in 1:nit) {
      res <- res + mixdistfun(
        xs,
        locations_list[[it]],
        scales_list[[it]],
        weights_list[[it]],
        distr.k,
        distfun
      )
    }
    return(res / nit)
  }


mixdistfun <-
  function(xs,
             locations,
             scales,
             weights,
             distr.k,
             distfun) {
    res <- 0.0 * xs
    for (cmp in seq_along(locations)) {
      res <-
        res + weights[cmp] * distfun(xs,
          distr = distr.k,
          mu = locations[cmp],
          sigma = scales[cmp]
        )
    }
    return(res)
  }

dmix <- function(xs, locations, scales, weights, distr.k) {
  mixdistfun(xs, locations, scales, weights, distr.k, dk)
}
pmix <- function(xs, locations, scales, weights, distr.k) {
  mixdistfun(xs, locations, scales, weights, distr.k, pk)
}

mixdistfun_cens <-
  function(xlefts,
             xrights,
             c_code_filters,
             locations,
             scales,
             weights,
             distr.k,
             distfun) {
    res <- 0.0 * xlefts
    for (cmp in seq_along(locations)) {
      res <-
        res + weights[cmp] * distfun(xlefts,
          xrights,
          c_code_filters,
          distr = distr.k,
          mu = locations[cmp],
          sigma = scales[cmp]
        )
    }
    return(res)
  }
dmixcens <- function(xlefts,
                     xrights,
                     c_code_filters,
                     locations,
                     scales,
                     weights,
                     distr.k) {
  mixdistfun_cens(xlefts, xrights, c_code_filters, locations, scales, weights, distr.k, dkcens2)
}
