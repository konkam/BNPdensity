

#' Acidity Index Dataset
#'
#' Concerns an acidity index measured in a sample of 155 lakes in north-central
#' Wisconsin.
#'
#'
#' @name acidity
#' @docType data
#' @format A real vector with 155 observations.
#' @references Crawford, S. L., DeGroot, M. H., Kadane, J. B. and Small, M. J.
#' (1992). Modeling lake chemistry distributions: approximate Bayesian methods
#' for estimating a finite mixture model. Technometrics, 34, 441-453.
#' @keywords datasets
#' @examples
#'
#' data(acidity)
#' hist(acidity)
NULL





#' Bayesian nonparametric density estimation
#'
#' This package performs Bayesian nonparametric density estimation for exact
#' and censored data via a normalized random measure mixture model. The package
#' allows the user to specify the mixture kernel, the mixing normalized measure
#' and the choice of performing fully nonparametric mixtures on locations and
#' scales, or semiparametric mixtures on locations only with common scale
#' parameter. Options for the kernels are: two kernels with support in the real
#' line (gaussian and double exponential), two more kernels in the positive
#' line (gamma and lognormal) and one with bounded support (beta). The options
#' for the normalized random measures are members of the class of normalized
#' generalized gamma, which include the Dirichlet process, the normalized
#' inverse gaussian process and the normalized stable process.  The type of
#' censored data handled by the package is right, left and interval.
#'
#' \tabular{ll}{ Package: \tab BNPdensity\cr Type: \tab Package\cr Version:
#' \tab 2016.10\cr Date: \tab 2016-10-14\cr License: \tab GPL version 2 or
#' later\cr LazyLoad: \tab yes\cr } The package includes four main functions:
#' MixNRMI1, MixNRMI2, MixNRMI1cens and MixNRMI2cens which implement
#' semiparametric and fully nonparametric mixtures for exact data, and
#' semiparametric and fully nonparametric mixtures for censored data
#' respectively. Additionally, the package includes several other functions
#' required for sampling from conditional distributions in the MCMC
#' implementation. These functions are intended for internal use only.
#'
#' @name BNPdensity-package
#' @aliases BNPdensity-package BNPdensity
#' @docType package
#' @author Barrios, E., Lijoi, A., Nieto-Barajas, L. E. and Prüenster, I.;
#' Contributor: Guillaume Kon Kam King.; Maintainer: Ernesto Barrios <ebarrios
#' at itam.mx>
#' @seealso \code{\link{MixNRMI1}}, \code{\link{MixNRMI2}},
#' \code{\link{MixNRMI1cens}}, \code{\link{MixNRMI2cens}}
#' @references Barrios, E., Lijoi, A., Nieto-Barajas, L. E. and Prüenster, I.
#' (2013). Modeling with Normalized Random Measure Mixture Models. Statistical
#' Science. Vol. 28, No. 3, 313-334.
#'
#' Kon Kam King, G., Arbel, J. and Prüenster, I. (2016). Species Sensitivity
#' Distribution revisited: a Bayesian nonparametric approach. In preparation.
#' @keywords package
#' @examples
#'
#' example(MixNRMI1)
#' example(MixNRMI2)
#' example(MixNRMI1cens)
#' example(MixNRMI2cens)
NULL





#' Fit of MixNRMI1 function to the enzyme dataset
#'
#' This object contains the output when setting set.seed(150520) and running
#' the function MixNRMI1(x, Alpha = 1, Kappa = 0.007, Gama = 0.5, distr.k = 2,
#' distr.p0 = 2, asigma = 1, bsigma = 1, Meps=0.001, Nit = 5000, Pbi = 0.2)
#'
#' See function MixNRMI1
#'
#' @name Enzyme1.out
#' @docType data
#' @keywords datasets
#' @examples
#'
#' data(Enzyme1.out)
NULL





#' Fit of MixNRMI2 function to the enzyme dataset
#'
#' This object contains the output when setting set.seed(150520) and running
#' the function Enzyme2.out <- MixNRMI2(x, Alpha = 1, Kappa = 0.007, Gama =
#' 0.5, distr.k = 2, distr.py0 = 2, distr.pz0 = 2, mu.pz0 = 1, sigma.pz0 = 1,
#' Meps=0.001, Nit = 5000, Pbi = 0.2)
#'
#' See function MixNRMI2
#'
#' @name Enzyme2.out
#' @docType data
#' @keywords datasets
#' @examples
#'
#' data(Enzyme2.out)
NULL





#' Enzyme Dataset
#'
#' Concerns the distribution of enzymatic activity in the blood, for an enzyme
#' involved in the metabolism of carcinogenetic substances, among a group of
#' 245 unrelated individuals.
#'
#'
#' @name enzyme
#' @docType data
#' @format A data frame with 244 observations on the following variable:
#' \describe{ \item{list("enzyme")}{A numeric vector.} }
#' @references Bechtel, Y. C., Bonaiti-Pellie, C., Poisson, N., Magnette, J.
#' and Bechtel, P.R. (1993). A population and family study of
#' N-acetyltransferase using caffeine urinary metabolites. Clin. Pharm. Therp.,
#' 54, 134-141.
#' @keywords datasets
#' @examples
#'
#' data(enzyme)
#' hist(enzyme)
NULL





#' Fit of MixNRMI1 function to the galaxy dataset
#'
#' This object contains the output when setting set.seed(150520) and running
#' the function MixNRMI1(x, Alpha = 1, Kappa = 0.015, Gama = 0.5, distr.k = 1,
#' distr.p0 = 2, asigma = 1, bsigma = 1, Meps=0.001, Nit = 5000, Pbi = 0.2)
#'
#' See function MixNRMI1.
#'
#' @name Galaxy1.out
#' @docType data
#' @keywords datasets
#' @examples
#'
#' data(Galaxy1.out)
NULL





#' Fit of MixNRMI2 function to the galaxy dataset
#'
#' This object contains the output when setting set.seed(150520) and running
#' the function MixNRMI2(x, Alpha = 1, Kappa = 0.015, Gama = 0.5, distr.k = 1,
#' distr.py0 = 2, distr.pz0 = 2, mu.pz0 = 1, sigma.pz0 = 1, Meps=0.001, Nit =
#' 5000, Pbi = 0.2)
#'
#' See function MixNRMI2.
#'
#' @name Galaxy2.out
#' @docType data
#' @keywords datasets
#' @examples
#'
#' data(Galaxy2.out)
NULL





#' Galaxy Data Set
#'
#' Velocities of 82 galaxies diverging from our own galaxy.
#'
#'
#' @name galaxy
#' @docType data
#' @format A data frame with 82 observations on the following variable:
#' \describe{ \item{list("velocity")}{A numeric vector.} }
#' @references Roeder, K. (1990) "Density estimation with confidence sets
#' exemplified by superclusters and voids in the galaxies". Journal of the
#' American Statistical Association. 85, 617-624.
#' @keywords datasets
#' @examples
#'
#' data(galaxy)
#' hist(galaxy)
NULL

#' Salinity tolerance
#'
#' 72-hour acute salinity tolerance (LC50 values) of riverine
#' macro-invertebrates.
#'
#'
#' @name salinity
#' @docType data
#' @format A data frame with 108 observations on the following two variables:
#' \describe{
#' \item{left}{A numeric vector.}
#' \item{right}{A
#' numeric vector.} }
#' @references Kefford, B.J., Nugegoda, D., Metzeling, L., Fields, E. 2006.
#' Validating species sensitivity distributions using salinity tolerance of
#' riverine macroinvertebrates in the southern Murray-darling Basin (Victoria,
#' Australia). Canadian Journal of Fisheries and Aquatic Science, 63,
#' 1865-1877.
#' @source \code{fitdistrplus} R-package
#' @keywords datasets
#' @examples
#'
#' data(salinity)
#' hist(salinity$left)
NULL
