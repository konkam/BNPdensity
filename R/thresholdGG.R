# This function uses the M_array which provides the threshold which ensures
# a moment match of 5% for NGG parameters alpha, kappa, gama


thresholdGG <-
  function(alpha = 1, kappa = 1, gama = 1 / 2, max_threshold = 200) {
    alpha_vect <- c(.1, 1, 5, 20) # mass param
    kappa_vect <- c(.1, 1, 5, 20)
    gama_vect <- c(0, .2, .4, .6)
    alpha_index <- which.max(alpha_vect >= alpha)
    kappa_index <- which.max(kappa_vect >= kappa)
    gama_index <- which.max(gama_vect >= gama)
    M <- M_array[alpha_index, kappa_index, gama_index]
    # if we are out of the grid, we assign the max_threshold
    out_of_grid <- (prod(1 - (alpha_vect >= alpha)) |
      prod(1 - (kappa_vect >= kappa)) |
      prod(1 - (gama_vect >= gama)))
    if (out_of_grid) {
      M <- max_threshold
    }
    return(M)
  }

globalVariables(names = c("M_array"))
