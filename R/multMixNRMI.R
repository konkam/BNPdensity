#' Multiple chains of MixNRMI1
#'
#'
#' @param x Numeric vector. Data set to which the density is fitted.
#' @param probs Numeric vector. Desired quantiles of the density estimates.
#' @param Alpha Numeric constant. Total mass of the centering measure. See details.
#' @param Kappa Numeric positive constant. See details.
#' @param Gama Numeric constant. 0 <= Gama <= 1. See details.
#' @param distr.k Integer number identifying the mixture kernel: 1 = Normal; 2 = Gamma; 3 = Beta; 4 = Double Exponential; 5 = Lognormal.
#' @param distr.p0 Integer number identifying the centering measure: 1 = Normal; 2 = Gamma; 3 = Beta.
#' @param asigma Numeric positive constant. Shape parameter of the gamma prior on the standard deviation of the mixture kernel distr.k.
#' @param bsigma Numeric positive constant. Rate parameter of the gamma prior on the standard deviation of the mixture kernel distr.k.
#' @param delta Numeric positive constant. Metropolis-Hastings proposal variation coefficient for sampling sigma.
#' @param Delta Numeric positive constant. Metropolis-Hastings proposal variation coefficient for sampling the latent U.
#' @param Meps Numeric constant. Relative error of the jump sizes in the continuous component of the process. Smaller values imply larger number of jumps.
#' @param Nx Integer constant. Number of grid points for the evaluation of the density estimate.
#' @param Nit Integer constant. Number of MCMC iterations.
#' @param Pbi Numeric constant. Burn-in period proportion of Nit.
#' @param epsilon Numeric constant. Extension to the evaluation grid range. See details.
#' @param printtime Logical. If TRUE, prints out the execution time.
#' @param extras Logical. If TRUE, gives additional objects: means, weights and Js (the jump sizes).
#' @param nchains The number of chains to run.
#' @param parallel Whether to run the chains in parallel. Only works on UNIX-like systems as it rests on Fork parallelism
#' @param ncores Number of cores for the parallel run. Defaults to parallel::detectCores(), i.e. the maximum number of cores detected by R on your system.
#'
#' @return a list containing the multiple fits.
#' @export
#'
#' @examples
#' data(acidity)
#' multMixNRMI1(acidity, parallel = TRUE, Nit = 10, ncores = 2)
multMixNRMI1 <- function(x, probs = c(0.025, 0.5, 0.975), Alpha = 1, Kappa = 0,
                         Gama = 0.4, distr.k = 1, distr.p0 = 1, asigma = 0.5, bsigma = 0.5,
                         delta = 3, Delta = 2, Meps = 0.01, Nx = 150, Nit = 1500,
                         Pbi = 0.1, epsilon = NULL, printtime = TRUE, extras = TRUE,
                         nchains = 4, parallel = TRUE, ncores = parallel::detectCores()) {
  if (Sys.info()[["sysname"]] == "Windows") parallel <- FALSE
  parallel::mclapply(
    X = 1:nchains,
    FUN = function(chainID) {
      MixNRMI1(
        x, probs, Alpha, Kappa,
        Gama, distr.k, distr.p0, asigma, bsigma,
        delta, Delta, Meps, Nx, Nit, Pbi,
        epsilon, printtime, extras
      )
    },
    mc.cores = ifelse(test = parallel, yes = ncores, no = 1),
    mc.set.seed = TRUE
  )
}



#' Multiple chains of MixNRMI2
#'
#' @param x Numeric vector. Data set to which the density is fitted.
#' @param probs Numeric vector. Desired quantiles of the density estimates.
#' @param Alpha Numeric constant. Total mass of the centering measure. See details.
#' @param Kappa Numeric positive constant. See details.
#' @param Gama Numeric constant. 0 <= Gama <= 1. See details.
#' @param distr.k Integer number identifying the mixture kernel: 1 = Normal; 2 = Gamma; 3 = Beta; 4 = Double Exponential; 5 = Lognormal.
#' @param distr.py0	Integer number identifying the centering measure for locations: 1 = Normal; 2 = Gamma; 3 = Beta.
#' @param distr.pz0 Integer number identifying the centering measure for scales: 2 = Gamma. For more options use MixNRMI2cens.
#' @param mu.pz0 Numeric constant. Prior mean of the centering measure for scales.
#' @param sigma.pz0 Numeric constant. Prior standard deviation of the centering measure for scales.
#' @param delta Numeric positive constant. Metropolis-Hastings proposal variation coefficient for sampling the scales.
#' @param kappa Numeric positive constant. Metropolis-Hastings proposal variation coefficient for sampling the location parameters.
#' @param Delta Numeric positive constant. Metropolis-Hastings proposal variation coefficient for sampling the latent U.
#' @param Meps Numeric constant. Relative error of the jump sizes in the continuous component of the process. Smaller values imply larger number of jumps.
#' @param Nx Integer constant. Number of grid points for the evaluation of the density estimate.
#' @param Nit Integer constant. Number of MCMC iterations.
#' @param Pbi Numeric constant. Burn-in period proportion of Nit.
#' @param epsilon Numeric constant. Extension to the evaluation grid range. See details.
#' @param printtime Logical. If TRUE, prints out the execution time.
#' @param extras Logical. If TRUE, gives additional objects: means, sigmas, weights and Js.
#' @param nchains The number of chains to run.
#' @param parallel Whether to run the chains in parallel. Only works on UNIX-like systems as it rests on Fork parallelism
#' @param ncores Number of cores for the parallel run. Defaults to parallel::detectCores(), i.e. the maximum number of cores detected by R on your system.
#'
#' @return a list containing the multiple fits.
#' @export
#'
#' @examples
#' data(acidity)
#' multMixNRMI2(acidity, parallel = TRUE, Nit = 10, ncores = 2)
multMixNRMI2 <- function(x, probs = c(0.025, 0.5, 0.975), Alpha = 1, Kappa = 0,
                         Gama = 0.4, distr.k = 1, distr.py0 = 1, distr.pz0 = 2, mu.pz0 = 3,
                         sigma.pz0 = sqrt(10), delta = 4, kappa = 2, Delta = 2, Meps = 0.01,
                         Nx = 150, Nit = 1500, Pbi = 0.1, epsilon = NULL, printtime = TRUE, extras = TRUE,
                         nchains = 4, parallel = FALSE, ncores = parallel::detectCores()) {
  if (Sys.info()[["sysname"]] == "Windows") parallel <- FALSE

  parallel::mclapply(
    X = 1:nchains,
    FUN = function(chainID) {
      MixNRMI2(
        x, probs, Alpha, Kappa,
        Gama, distr.k, distr.py0, distr.pz0, mu.pz0,
        sigma.pz0, delta, kappa, Delta, Meps,
        Nx, Nit, Pbi, epsilon, printtime, extras
      )
    },
    mc.cores = ifelse(test = parallel, yes = ncores, no = 1),
    mc.set.seed = TRUE
  )
}


#' Multiple chains of MixNRMI1cens
#'
#'
#' @inherit multMixNRMI1
#' @param xleft Numeric vector. Lower limit of interval censoring. For exact data the same as xright
#' @param xright Numeric vector. Upper limit of interval censoring. For exact data the same as xleft.
#'
#' @return a list containing the multiple fits.
#' @export
#'
#' @examples
#' data(salinity)
#' multMixNRMI1cens(salinity$left, salinity$right, parallel = TRUE, Nit = 10, ncores = 2)
multMixNRMI1cens <- function(xleft, xright, probs = c(0.025, 0.5, 0.975), Alpha = 1, Kappa = 0,
                             Gama = 0.4, distr.k = 1, distr.p0 = 1, asigma = 0.5, bsigma = 0.5,
                             delta = 3, Delta = 2, Meps = 0.01, Nx = 150, Nit = 1500,
                             Pbi = 0.1, epsilon = NULL, printtime = TRUE, extras = TRUE,
                             nchains = 4, parallel = TRUE, ncores = parallel::detectCores()) {
  if (Sys.info()[["sysname"]] == "Windows") parallel <- FALSE

  parallel::mclapply(
    X = 1:nchains,
    FUN = function(chainID) {
      MixNRMI1cens(
        xleft, xright, probs, Alpha, Kappa,
        Gama, distr.k, distr.p0, asigma, bsigma,
        delta, Delta, Meps, Nx, Nit, Pbi,
        epsilon, printtime, extras
      )
    },
    mc.cores = ifelse(test = parallel, yes = ncores, no = 1),
    mc.set.seed = TRUE
  )
}

#' Multiple chains of MixNRMI2cens
#'
#'
#' @inherit multMixNRMI2
#' @param xleft Numeric vector. Lower limit of interval censoring. For exact data the same as xright
#' @param xright Numeric vector. Upper limit of interval censoring. For exact data the same as xleft.
#'
#' @return a list containing the multiple fits.
#' @export
#'
#' @examples
#' data(salinity)
#' \dontrun{
#' multMixNRMI2cens(salinity$left, salinity$right, parallel = TRUE, Nit = 2, ncores = 2)
#' }
multMixNRMI2cens <- function(xleft, xright, probs = c(0.025, 0.5, 0.975), Alpha = 1,
                             Kappa = 0, Gama = 0.4, distr.k = 1, distr.py0 = 1, distr.pz0 = 2,
                             mu.pz0 = 3, sigma.pz0 = sqrt(10), delta = 4, kappa = 2, Delta = 2,
                             Meps = 0.01, Nx = 150, Nit = 1500, Pbi = 0.1, epsilon = NULL,
                             printtime = TRUE, extras = TRUE,
                             nchains = 4, parallel = TRUE, ncores = parallel::detectCores()) {
  if (Sys.info()[["sysname"]] == "Windows") parallel <- FALSE

  parallel::mclapply(
    X = 1:nchains,
    FUN = function(chainID) {
      MixNRMI2cens(
        xleft, xright, probs, Alpha, Kappa,
        Gama, distr.k, distr.py0, distr.pz0, mu.pz0,
        sigma.pz0, delta, kappa, Delta, Meps,
        Nx, Nit, Pbi, epsilon, printtime, extras
      )
    },
    mc.cores = ifelse(test = parallel, yes = ncores, no = 1),
    mc.set.seed = TRUE
  )
}
