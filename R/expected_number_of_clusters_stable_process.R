log_Vnk_PY <- function(n, k, Alpha, Gama) {
  if (k == 1) {
    lognum <- 0
  }
  else {
    lognum <- sum(log(Alpha + Gama * 1:(k - 1)))
  }
  return(lognum - sum(log(Alpha + 1 + 0:(n - 2))))
}


# log_Vnk_PY(7, 6, 0.5, 0.01)
# log_Vnk_PY(6, 5, 0.5, 0.001)

# Pochhammer = function(x, n){
#   x_mpfr = as.bigq(x)
#   gamma(x_mpfr+n)/gamma(x_mpfr)
# }

Cnk <- function(n, k, Gama) {
  # factor_k <- gmp::factorialZ(k)
  # Using more precise sums (package PreciseSums did not afford any improvement)
  # This function still risks underflow/overflow in spite of the arbitrary precision packages
  # sum((-1)^(1:k) * gmp::chooseZ(n = k, k = 1:k) * Rmpfr::pochMpfr(-(1:k) * Gama, n) / factor_k)
  (-1)^(n-k)*noncentral_generalised_factorial_coefficient(n = n, k = k, s = Gama, r = 0)
}
log_Cnk = function(n, k, Gama){
  log(Cnk(n, k, Gama))
}

# log_Cnk(100, 64, 0.4)
#
#
# log_Cnk(6, 5, 0.5)
#
# library(Rmpfr)
# library(gmp)


Pkn_PY <- function(k, n, Alpha, Gama) {
  exp(log_Vnk_PY(n = n, k = k, Alpha = Alpha, Gama = Gama) - k*log(Gama) + log_Cnk(n = n, k = k, Gama = Gama))
  # Using this form, the inaccuracies in Cnk (i.e. getting a negative number) do not give NaN.
  # This error might be cancelled when computing the expected number of components.
  # exp(log_Vnk_PY(n = n, k = k, Alpha = Alpha, Gama = Gama) - k * log(Gama)) * Cnk(n = n, k = k, Gama = Gama)
}

# Pkn_PY(3, 5, 0.2, 0.4)

Pkn_Dirichlet <- function(k, n, Alpha) {
  exp(k * log(Alpha) + log(abs(gmp::Stirling1(n, k))) - sum(log(Alpha + 0:(n - 1))))
}

# Pkn_Dirichlet(3, 5, 0.2)


expected_number_of_components_PY <- function(n, Alpha, Gama, ntrunc = NULL) {
  if (is.null(ntrunc)) {
    ntrunc <- n
  } else if (ntrunc > n) ntrunc <- n
  res <- 0
  for (k in 1:ntrunc) {
    # print(k)
    # print(k*Pkn_PY(k, n, Alpha, Gama))
    res <- res + k * Pkn_PY(k, n, Alpha, Gama)
  }
  return(res)
}



#' Computes the expected number of components for a Dirichlet process.
#'
#' @export
#'
#' @param n Number of data points
#' @param Alpha Numeric constant. Total mass of the centering measure.
#' @param ntrunc Level of truncation when computing the expectation. Defaults to n. If greater than n, it is fixed to n.
#'
#' @return A real value which approximates the expected number of components
#'
#' Reference: P. De Blasi, S. Favaro, A. Lijoi, R. H. Mena, I. Prünster, and M. Ruggiero, “Are gibbs-type priors the most natural generalization of the dirichlet process?,” IEEE Trans. Pattern Anal. Mach. Intell., vol. 37, no. 2, pp. 212–229, 2015.
#' @examples
#' expected_number_of_components_Dirichlet(100, 1.2)
expected_number_of_components_Dirichlet <- function(n, Alpha, ntrunc = NULL) {
  if (is.null(ntrunc)) {
    ntrunc <- n
  } else if (ntrunc > n) ntrunc <- n
  res <- 0
  for (k in 1:ntrunc) {
    # print(k)
    # print(k*Pkn_PY(k, n, Alpha, Gama))
    res <- res + k * Pkn_Dirichlet(k, n, Alpha)
  }
  return(res)
}

# expected_number_of_components_PY(10, 0., 0.4)

#' Computes the expected number of components for a stable process.
#' @export
#'
#' @param n Number of data points
#' @param Gama Numeric constant. 0 <= Gama <=1.
#' @param ntrunc Level of truncation when computing the expectation. Defaults to n. If greater than n, it is fixed to n.
#'
#' @return A real value of type mpfr1 which approximates the expected number of components
#'
#' In spite of the high precision arithmetics packages used for in function, it can be numerically unstable for small values of Gama.
#' This is because evaluating a sum with alternated signs, in the generalised factorial coefficients, is tricky.
#' Reference: P. De Blasi, S. Favaro, A. Lijoi, R. H. Mena, I. Prünster, and M. Ruggiero, “Are gibbs-type priors the most natural generalization of the dirichlet process?,” IEEE Trans. Pattern Anal. Mach. Intell., vol. 37, no. 2, pp. 212–229, 2015.
#'
#' @examples
#' expected_number_of_components_stable(100, 0.8)
expected_number_of_components_stable <- function(n, Gama, ntrunc = NULL) {
  if (!requireNamespace("Rmpfr", quietly = TRUE) && !requireNamespace("gmp", quietly = TRUE)) {
    stop("Packages Rmpfr and gmp are needed for this function to work. Please install them.",
      call. = FALSE
    )
  }
  expected_number_of_components_PY(n, 0, Gama, ntrunc = ntrunc)
}

#' Plot the prior number of components for a stable process and for a Dirichlet process with Alpha = 1.
#'
#' This plots the prior distribution on the number of components for the stable process. The Dirichlet process is provided for comparison.
#'
#' @param n Number of data points
#' @param Gama Numeric constant. 0 <= Gama <=1.
#' @param Alpha Numeric constant. Total mass of the centering measure for the Dirichlet process.
#' @param grid Level of truncation when computing the expectation. Defaults to n. If greater than n, it is fixed to n.
#'
#' @return A plot with the prior distribution on the number of components.
#' @export
#'
#' @examples
#' plot_prior_number_of_components(50, 0.4)
plot_prior_number_of_components <- function(n, Gama, Alpha = 1, grid = NULL) {
  if (!requireNamespace("Rmpfr", quietly = TRUE) && !requireNamespace("gmp", quietly = TRUE)) {
    stop("Packages Rmpfr and gmp are needed for this function to work. Please install them.",
      call. = FALSE
    )
  }
  if (is.null(grid)) grid <- 1:n
  grid <- unique(round(grid)) # Make sure it is a grid of integers
  to_plot <- rbind(
    data.frame(K = grid, Pk = Rmpfr::toNum(Vectorize(Pkn_PY, vectorize.args = "k")(grid, n, 0, Gama)), Process = "Stable"),
    data.frame(K = grid, Pk = Vectorize(Pkn_Dirichlet, vectorize.args = "k")(grid, n, Alpha), Process = "Dirichlet")
  )
  ggplot(data = to_plot, aes_string(x = "K", y = "Pk", colour = "factor(Process)", group = "Process")) +
    geom_point() +
    geom_line() +
    theme_classic() +
    viridis::scale_colour_viridis(discrete = T, name = "Process") +
    ylab(expression(P[K]))
}
