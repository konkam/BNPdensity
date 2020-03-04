library(BNPdensity)
library(usethis)

set.seed(150520)
data(enzyme)
Enzyme1.out <- MixNRMI1(enzyme, Alpha = 1, Kappa = 0.007, Gama = 0.5, distr.k = "gamma", distr.p0 = "gamma", asigma = 1, bsigma = 1, Meps = 0.005, Nit = 5000, Pbi = 0.8)

# # This is for diagnostic
# attach(Enzyme1.out)
# # Plotting density estimate + 95% credible interval
# m <- ncol(qx)
# ymax <- max(qx[, m])
# par(mfrow = c(1, 1))
# hist(enzyme, probability = TRUE, breaks = 20, col = grey(.9), ylim = c(0, ymax))
# lines(xx, qx[, 1], lwd = 2)
# lines(xx, qx[, 2], lty = 3, col = 4)
# lines(xx, qx[, m], lty = 3, col = 4)
# # Plotting number of clusters
# par(mfrow = c(2, 1))
# plot(R, type = "l", main = "Trace of R")
# hist(R, breaks = min(R - 0.5):max(R + 0.5), probability = TRUE)
# # Plotting sigma
# par(mfrow = c(2, 1))
# plot(S, type = "l", main = "Trace of sigma")
# hist(S, nclass = 20, probability = TRUE, main = "Histogram of sigma")
# print(length(unique(S))/length(S))

usethis::use_data(Enzyme1.out, overwrite = T)

set.seed(150520)
data(enzyme)
Enzyme2.out <- MixNRMI2(enzyme, Alpha = 1, Kappa = 0.007, Gama = 0.5, distr.k = "gamma", distr.py0 = "gamma", distr.pz0 = "gamma", mu.pz0 = 1, sigma.pz0 = 1, Meps = 0.005, Nit = 5000, Pbi = 0.8)

# # This is for diagnostic
# attach(Enzyme2.out)
# # Plotting density estimate + 95% credible interval
# m <- ncol(qx)
# ymax <- max(qx[, m])
# par(mfrow = c(1, 1))
# hist(enzyme, probability = TRUE, breaks = 20, col = grey(.9), ylim = c(0, ymax))
# lines(xx, qx[, 1], lwd = 2)
# lines(xx, qx[, 2], lty = 3, col = 4)
# lines(xx, qx[, m], lty = 3, col = 4)
# # Plotting number of clusters
# par(mfrow = c(2, 1))
# plot(R, type = "l", main = "Trace of R")
# hist(R, breaks = min(R - 0.5):max(R + 0.5), probability = TRUE)

usethis::use_data(Enzyme2.out, overwrite = T)

set.seed(150520)
data(galaxy)
Galaxy1.out <- MixNRMI1(galaxy, Alpha = 1, Kappa = 0.015, Gama = 0.5, distr.k = "normal", distr.p0 = "gamma", asigma = 1, bsigma = 1, delta_S = 7, Meps = 0.005, Nit = 5000, Pbi = 0.8)

# # This is for diagnostic
# attach(Galaxy1.out)
# # Plotting density estimate + 95% credible interval
# m <- ncol(qx)
# ymax <- max(qx[, m])
# par(mfrow = c(1, 1))
# hist(galaxy, probability = TRUE, breaks = 20, col = grey(.9), ylim = c(0, ymax))
# lines(xx, qx[, 1], lwd = 2)
# lines(xx, qx[, 2], lty = 3, col = 4)
# lines(xx, qx[, m], lty = 3, col = 4)
# # Plotting number of clusters
# par(mfrow = c(2, 1))
# plot(R, type = "l", main = "Trace of R")
# hist(R, breaks = min(R - 0.5):max(R + 0.5), probability = TRUE)
# # Plotting sigma
# par(mfrow = c(2, 1))
# plot(S, type = "l", main = "Trace of sigma")
# hist(S, nclass = 20, probability = TRUE, main = "Histogram of sigma")
# print(length(unique(S))/length(S))

usethis::use_data(Galaxy1.out, overwrite = T)

set.seed(150520)
data(galaxy)
Galaxy2.out <- MixNRMI2(galaxy, Alpha = 1, Kappa = 0.015, Gama = 0.5, distr.k = "normal", distr.py0 = "gamma", distr.pz0 = "gamma", mu.pz0 = 1, sigma.pz0 = 1,  Meps = 0.005, Nit = 5000, Pbi = 0.8)

# # This is for diagnostic
# attach(Galaxy2.out)
# # Plotting density estimate + 95% credible interval
# m <- ncol(qx)
# ymax <- max(qx[, m])
# par(mfrow = c(1, 1))
# hist(galaxy, probability = TRUE, breaks = 20, col = grey(.9), ylim = c(0, ymax))
# lines(xx, qx[, 1], lwd = 2)
# lines(xx, qx[, 2], lty = 3, col = 4)
# lines(xx, qx[, m], lty = 3, col = 4)
# # Plotting number of clusters
# par(mfrow = c(2, 1))
# plot(R, type = "l", main = "Trace of R")
# hist(R, breaks = min(R - 0.5):max(R + 0.5), probability = TRUE)

usethis::use_data(Galaxy2.out, overwrite = T)
