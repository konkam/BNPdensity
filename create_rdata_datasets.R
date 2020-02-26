library(BNPdensity)
library(usethis)

set.seed(150520)
data(enzyme)
x <- enzyme
Enzyme1.out <- MixNRMI1(x,
  Alpha = 1, Kappa = 0.007, Gama = 0.5,
  distr.k = 2, distr.p0 = 2, asigma = 1, bsigma = 1, Meps = 0.005,
  Nit = 5000, Pbi = 0.2
)
usethis::use_data(Enzyme1.out, overwrite = T)

set.seed(150520)
data(galaxy)
x <- galaxy
Galaxy1.out <- MixNRMI1(x, Alpha = 1, Kappa = 0.015, Gama = 0.5, distr.k = 1, distr.p0 = 2, asigma = 1, bsigma = 1, Meps = 0.005, Nit = 5000, Pbi = 0.2)
usethis::use_data(Galaxy1.out, overwrite = T)


set.seed(150520)
data(enzyme)
x <- enzyme
Enzyme2.out <- MixNRMI2(x,
  Alpha = 1, Kappa = 0.007, Gama = 0.5,
  distr.k = 2, distr.py0 = 2,
  distr.pz0 = 2, mu.pz0 = 1, sigma.pz0 = 1, Meps = 0.005,
  Nit = 5000, Pbi = 0.2
)
usethis::use_data(Enzyme2.out, overwrite = T)

set.seed(150520)
data(galaxy)
x <- galaxy
Galaxy2.out <- MixNRMI2(x, Alpha = 1, Kappa = 0.015, Gama = 0.5,
                          distr.k = 1, distr.py0 = 2,
                          distr.pz0 = 2, mu.pz0 = 1, sigma.pz0 = 1,  Meps=0.005,
                          Nit = 5000, Pbi = 0.2)
usethis::use_data(Galaxy2.out, overwrite = T)
