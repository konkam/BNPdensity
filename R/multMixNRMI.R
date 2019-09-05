#' Multiple chains of MixNRMI1
#'
#' @param x Numeric vector. Data set to which the density is fitted.
#' @param probs Numeric vector. Desired quantiles of the density estimates.
#' @param Alpha Numeric constant. Total mass of the centering measure. See
#' details.
#' @param Kappa Numeric positive constant. See details.
#' @param Gama Numeric constant. 0 <= Gama <= 1. See details.
#' @param distr.k Integer number identifying the mixture kernel: 1 = Normal; 2
#' = Gamma; 3 = Beta; 4 = Double Exponential; 5 = Lognormal.
#' @param distr.p0 Integer number identifying the centering measure: 1 =
#' Normal; 2 = Gamma; 3 = Beta.
#' @param asigma Numeric positive constant. Shape parameter of the gamma prior
#' on the standard deviation of the mixture kernel distr.k.
#' @param bsigma Numeric positive constant. Rate parameter of the gamma prior
#' on the standard deviation of the mixture kernel distr.k.
#' @param delta Numeric positive constant. Metropolis-Hastings proposal
#' variation coefficient for sampling sigma.
#' @param Delta Numeric positive constant. Metropolis-Hastings proposal
#' variation coefficient for sampling the latent U.
#' @param Meps Numeric constant. Relative error of the jump sizes in the
#' continuous component of the process. Smaller values imply larger number of
#' jumps.
#' @param Nx Integer constant. Number of grid points for the evaluation of the
#' density estimate.
#' @param Nit Integer constant. Number of MCMC iterations.
#' @param Pbi Numeric constant. Burn-in period proportion of Nit.
#' @param epsilon Numeric constant. Extension to the evaluation grid range. See
#' details.
#' @param printtime Logical. If TRUE, prints out the execution time.
#' @param extras Logical. If TRUE, gives additional objects: means, weights and
#' Js (the jump sizes).
#' @param nchains The number of chains to run.
#' @param parallel Whether to run the chains in parallel. Only works on UNIX-like systems as it rests on Fork parallelism
#' @param ncores Number of cores for the parallel run. Defaults to parallel::detectCores(), i.e. the maximum number of cores detected by R on your system.
#'
#' @return a list containing the multiple fits.
#' @seealso \code{\link{MixNRMI2}}, \code{\link{MixNRMI1cens}},
#' \code{\link{MixNRMI2cens}}
#' @examples
#'
#' data(acidity)
#' multMixNRMI1(acidity, parallel = TRUE, Nit = 10, ncores = 2)
#' @export multMixNRMI1
multMixNRMI1 <- function(x, probs = c(0.025, 0.5, 0.975), Alpha = 1, Kappa = 0,
                         Gama = 0.4, distr.k = 1, distr.p0 = 1, asigma = 0.5, bsigma = 0.5,
                         delta = 3, Delta = 2, Meps = 0.01, Nx = 150, Nit = 1500,
                         Pbi = 0.1, epsilon = NULL, printtime = TRUE, extras = TRUE,
                         nchains = 4, parallel = TRUE, ncores = parallel::detectCores()) {
  if (Sys.info()[["sysname"]] == "Windows") parallel <- FALSE
  res <- parallel::mclapply(
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
  return(structure(res, class = c("multNRMI", "NRMI1")))
}





#' Multiple chains of MixNRMI2
#'
#' @param x Numeric vector. Data set to which the density is fitted.
#' @param probs Numeric vector. Desired quantiles of the density estimates.
#' @param Alpha Numeric constant. Total mass of the centering measure. See
#' details.
#' @param Kappa Numeric positive constant. See details.
#' @param Gama Numeric constant. 0 <= Gama <= 1. See details.
#' @param distr.k Integer number identifying the mixture kernel: 1 = Normal; 2
#' = Gamma; 3 = Beta; 4 = Double Exponential; 5 = Lognormal.
#' @param distr.py0 Integer number identifying the centering measure for
#' locations: 1 = Normal; 2 = Gamma; 3 = Beta.
#' @param distr.pz0 Integer number identifying the centering measure for
#' scales: 2 = Gamma. For more options use MixNRMI2cens.
#' @param mu.pz0 Numeric constant. Prior mean of the centering measure for
#' scales.
#' @param sigma.pz0 Numeric constant. Prior standard deviation of the centering
#' measure for scales.
#' @param delta Numeric positive constant. Metropolis-Hastings proposal
#' variation coefficient for sampling the scales.
#' @param kappa Numeric positive constant. Metropolis-Hastings proposal
#' variation coefficient for sampling the location parameters.
#' @param Delta Numeric positive constant. Metropolis-Hastings proposal
#' variation coefficient for sampling the latent U.
#' @param Meps Numeric constant. Relative error of the jump sizes in the
#' continuous component of the process. Smaller values imply larger number of
#' jumps.
#' @param Nx Integer constant. Number of grid points for the evaluation of the
#' density estimate.
#' @param Nit Integer constant. Number of MCMC iterations.
#' @param Pbi Numeric constant. Burn-in period proportion of Nit.
#' @param epsilon Numeric constant. Extension to the evaluation grid range. See
#' details.
#' @param printtime Logical. If TRUE, prints out the execution time.
#' @param extras Logical. If TRUE, gives additional objects: means, sigmas,
#' weights and Js.
#' @param nchains The number of chains to run.
#' @param parallel Whether to run the chains in parallel. Only works on UNIX-like systems as it rests on Fork parallelism
#' @param ncores Number of cores for the parallel run. Defaults to parallel::detectCores(), i.e. the maximum number of cores detected by R on your system.
#'
#' @return a list containing the multiple fits.
#' @seealso \code{\link{MixNRMI2}}, \code{\link{MixNRMI1cens}},
#' \code{\link{MixNRMI2cens}}, \code{\link{multMixNRMI1}}
#' @examples
#'
#' data(acidity)
#' multMixNRMI2(acidity, parallel = TRUE, Nit = 10, ncores = 2)
#' @export multMixNRMI2
multMixNRMI2 <- function(x, probs = c(0.025, 0.5, 0.975), Alpha = 1, Kappa = 0,
                         Gama = 0.4, distr.k = 1, distr.py0 = 1, distr.pz0 = 2, mu.pz0 = 3,
                         sigma.pz0 = sqrt(10), delta = 4, kappa = 2, Delta = 2, Meps = 0.01,
                         Nx = 150, Nit = 1500, Pbi = 0.1, epsilon = NULL, printtime = TRUE, extras = TRUE,
                         nchains = 4, parallel = FALSE, ncores = parallel::detectCores()) {
  if (Sys.info()[["sysname"]] == "Windows") parallel <- FALSE

  res <- parallel::mclapply(
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
  return(structure(res, class = c("multNRMI", "NRMI2")))
}




#' Multiple chains of MixNRMI1cens
#'
#' @param xleft Numeric vector. Lower limit of interval censoring. For exact
#' data the same as xright
#' @param xright Numeric vector. Upper limit of interval censoring. For exact
#' data the same as xleft.
#' @param probs Numeric vector. Desired quantiles of the density estimates.
#' @param Alpha Numeric constant. Total mass of the centering measure. See
#' details.
#' @param Kappa Numeric positive constant. See details.
#' @param Gama Numeric constant. 0 <= Gama <= 1. See details.
#' @param distr.k Integer number identifying the mixture kernel: 1 = Normal; 2
#' = Gamma; 3 = Beta; 4 = Double Exponential; 5 = Lognormal.
#' @param distr.p0 Integer number identifying the centering measure: 1 =
#' Normal; 2 = Gamma; 3 = Beta.
#' @param asigma Numeric positive constant. Shape parameter of the gamma prior
#' on the standard deviation of the mixture kernel distr.k.
#' @param bsigma Numeric positive constant. Rate parameter of the gamma prior
#' on the standard deviation of the mixture kernel distr.k.
#' @param delta Numeric positive constant. Metropolis-Hastings proposal
#' variation coefficient for sampling sigma.
#' @param Delta Numeric positive constant. Metropolis-Hastings proposal
#' variation coefficient for sampling the latent U.
#' @param Meps Numeric constant. Relative error of the jump sizes in the
#' continuous component of the process. Smaller values imply larger number of
#' jumps.
#' @param Nx Integer constant. Number of grid points for the evaluation of the
#' density estimate.
#' @param Nit Integer constant. Number of MCMC iterations.
#' @param Pbi Numeric constant. Burn-in period proportion of Nit.
#' @param epsilon Numeric constant. Extension to the evaluation grid range. See
#' details.
#' @param printtime Logical. If TRUE, prints out the execution time.
#' @param extras Logical. If TRUE, gives additional objects: means, weights and
#' Js (the jump sizes).
#' @param nchains The number of chains to run.
#' @param parallel Whether to run the chains in parallel. Only works on
#' UNIX-like systems as it rests on Fork parallelism
#' @param ncores Number of cores for the parallel run. Defaults to
#' parallel::detectCores(), i.e. the maximum number of cores detected by R on
#' your system.
#' @return a list containing the multiple fits.
#' @seealso \code{\link{MixNRMI2}}, \code{\link{MixNRMI1cens}},
#' \code{\link{MixNRMI2cens}}, \code{\link{multMixNRMI1}}
#' @examples
#'
#' data(salinity)
#' multMixNRMI1cens(salinity$left, salinity$right, parallel = TRUE, Nit = 10, ncores = 2)
#' @export multMixNRMI1cens
multMixNRMI1cens <- function(xleft, xright, probs = c(0.025, 0.5, 0.975), Alpha = 1, Kappa = 0,
                             Gama = 0.4, distr.k = 1, distr.p0 = 1, asigma = 0.5, bsigma = 0.5,
                             delta = 3, Delta = 2, Meps = 0.01, Nx = 150, Nit = 1500,
                             Pbi = 0.1, epsilon = NULL, printtime = TRUE, extras = TRUE,
                             nchains = 4, parallel = TRUE, ncores = parallel::detectCores()) {
  if (Sys.info()[["sysname"]] == "Windows") parallel <- FALSE

  res <- parallel::mclapply(
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
  return(structure(res, class = c("multNRMI", "NRMI1cens")))
}



#' Multiple chains of MixNRMI2cens
#'
#' @param xleft Numeric vector. Lower limit of interval censoring. For exact
#' data the same as xright
#' @param xright Numeric vector. Upper limit of interval censoring. For exact
#' data the same as xleft.
#' @param probs Numeric vector. Desired quantiles of the density estimates.
#' @param Alpha Numeric constant. Total mass of the centering measure. See
#' details.
#' @param Kappa Numeric positive constant. See details.
#' @param Gama Numeric constant. 0 <= Gama <= 1. See details.
#' @param distr.k Integer number identifying the mixture kernel: 1 = Normal; 2
#' = Gamma; 3 = Beta; 4 = Double Exponential; 5 = Lognormal.
#' @param distr.py0 Integer number identifying the centering measure for
#' locations: 1 = Normal; 2 = Gamma; 3 = Beta.
#' @param distr.pz0 Integer number identifying the centering measure for
#' scales: 2 = Gamma. For more options use MixNRMI2cens.
#' @param mu.pz0 Numeric constant. Prior mean of the centering measure for
#' scales.
#' @param sigma.pz0 Numeric constant. Prior standard deviation of the centering
#' measure for scales.
#' @param delta Numeric positive constant. Metropolis-Hastings proposal
#' variation coefficient for sampling the scales.
#' @param kappa Numeric positive constant. Metropolis-Hastings proposal
#' variation coefficient for sampling the location parameters.
#' @param Delta Numeric positive constant. Metropolis-Hastings proposal
#' variation coefficient for sampling the latent U.
#' @param Meps Numeric constant. Relative error of the jump sizes in the
#' continuous component of the process. Smaller values imply larger number of
#' jumps.
#' @param Nx Integer constant. Number of grid points for the evaluation of the
#' density estimate.
#' @param Nit Integer constant. Number of MCMC iterations.
#' @param Pbi Numeric constant. Burn-in period proportion of Nit.
#' @param epsilon Numeric constant. Extension to the evaluation grid range. See
#' details.
#' @param printtime Logical. If TRUE, prints out the execution time.
#' @param extras Logical. If TRUE, gives additional objects: means, sigmas,
#' weights and Js.
#' @param nchains The number of chains to run.
#' @param parallel Whether to run the chains in parallel. Only works on
#' UNIX-like systems as it rests on Fork parallelism
#' @param ncores Number of cores for the parallel run. Defaults to
#' parallel::detectCores(), i.e. the maximum number of cores detected by R on
#' your system.
#' @return a list containing the multiple fits.
#' @seealso \code{\link{MixNRMI2}}, \code{\link{MixNRMI1cens}},
#' \code{\link{MixNRMI2cens}}, \code{\link{multMixNRMI1}}
#' @examples
#'
#' data(salinity)
#' \dontrun{
#' multMixNRMI2cens(salinity$left, salinity$right, parallel = TRUE, Nit = 20, ncores = 2)
#' }
#'
#' @export multMixNRMI2cens
multMixNRMI2cens <- function(xleft, xright, probs = c(0.025, 0.5, 0.975), Alpha = 1,
                             Kappa = 0, Gama = 0.4, distr.k = 1, distr.py0 = 1, distr.pz0 = 2,
                             mu.pz0 = 3, sigma.pz0 = sqrt(10), delta = 4, kappa = 2, Delta = 2,
                             Meps = 0.01, Nx = 150, Nit = 1500, Pbi = 0.1, epsilon = NULL,
                             printtime = TRUE, extras = TRUE,
                             nchains = 4, parallel = TRUE, ncores = parallel::detectCores()) {
  if (Sys.info()[["sysname"]] == "Windows") parallel <- FALSE

  res <- parallel::mclapply(
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
  return(structure(res, class = c("multNRMI", "NRMI2cens")))
}



#' Convert the output of multMixNRMI into a coda mcmc object
#'
#' Convert the output of multMixNRMI into a coda mcmc object
#'
#' @importFrom  coda as.mcmc
#' @param fitlist Output of multMixNRMI.
#' @param thinning_to Final length of the chain after thinning.
#' @return a coda::mcmc object
#' @method as.mcmc multNRMI
#' @export
as.mcmc.multNRMI <- function(fitlist, thinning_to = 1000) {
  res <- coda::as.mcmc(lapply(Convert_to_matrix_list(fitlist, thinning_to = thinning_to), coda::mcmc))
  class(res) <- c("multNRMI", class(res))
  return(res)
}

#' Plot the density estimate and the 95\% credible interval
#'
#' The density estimate is the mean posterior density computed on the data
#' points.
#'
#'
#' @param x An object of class multNRMI
#' @param ... Further arguments to be passed to generic functions, ignored at the moment
#' @return A graph with the density estimate, the 95\% credible interval.
#' Includes a histogram if the data is non censored.
#' @export
#' @examples
#'
#' fit <- multMixNRMI2cens(salinity$left, salinity$right, parallel = TRUE, Nit = 20, ncores = 2)
#' plot(fit)
plot.multNRMI <- function(x, ...) {
  # This assumes that chains have the same length and can be given equal weight when combining
  res <- x[[1]]
  nchains <- length(x)
  m <- ncol(res$qx)
  res$qx[, 1] <- 1 / nchains * Reduce(f = add, lapply(X = x, FUN = function(x) x$qx[, 1]))
  res$qx[, 2] <- 1 / nchains * Reduce(f = add, lapply(X = x, FUN = function(x) x$qx[, 2]))
  res$qx[, m] <- 1 / nchains * Reduce(f = add, lapply(X = x, FUN = function(x) x$qx[, m]))
  plot(res)
}

#' S3 method for class 'multNRMI'
#'
#' @param x An object of class multNRMI
#' @param ... Further arguments to be passed to generic functions, ignored at the moment
#'
#' @return A visualisation of the important information about the object
#' @export
#'
#' @examples
#' data(salinity)
#' out <- multMixNRMI2cens(salinity$left, salinity$right, parallel = TRUE, Nit = 20, ncores = 2)
#' print(out)
print.multNRMI <- function(x, ...) {
  print(x[[1]])
  writeLines(paste(length(x), "independent MCMC chains were run in parallel"))
}

#' S3 method for class 'multNRMI'
#'
#' @param object A fitted object of class NRMI1cens
#' @param number_of_clusters Whether to compute the optimal number of clusters, which can be a time-consuming operation (see \code{\link{compute_optimal_clustering}})
#' @param ... Further arguments to be passed to generic function, ignored at the moment
#'
#' @return Prints out the text for the summary S3 methods
#' @export
#'
#' @examples
#' data(salinity)
#' out <- multMixNRMI2cens(salinity$left, salinity$right, parallel = TRUE, Nit = 20, ncores = 2)
#' summary(out)
summary.multNRMI <- function(object, number_of_clusters = FALSE, ...) {
  kernel_name <- tolower(give_kernel_name(object[[1]]$distr.k))
  NRMI_comment <- paste("Density estimation using a", comment_on_NRMI_type(object[[1]]$NRMI_params))
  kernel_comment <- paste("A nonparametric", kernel_name, "mixture model was used.")
  ndata <- ifelse(is_censored(object[[1]]$data), nrow(object[[1]]$data), length(object[[1]]$data))
  data_comment <- paste("There were", ndata, "data points.")
  n_chains <- length(object)
  MCMC_comment <- paste(n_chains, " MCMC chains were run for ", object[[1]]$Nit, " iterations with ", 100 * object[[1]]$Pbi, "% discarded for burn-in.", sep = "")
  if (number_of_clusters) {
    collected_allocs <- list("Allocs" = Reduce(c, lapply(object, function(x) x$Allocs)))
    estimated_clustering <- compute_optimal_clustering(collected_allocs)
    clustering_comment <- paste("The estimated number of clusters in the data is ", length(unique(estimated_clustering)), ".", sep = "")
  }
  else {
    clustering_comment <- "To obtain information on the estimated number of clusters, please use summary(object, number_of_clusters = TRUE)."
  }
  writeLines(paste(NRMI_comment, "\n", kernel_comment, "\n", data_comment, "\n", MCMC_comment, "\n", clustering_comment, sep = ""))
}
