#' @export
MixNRMI2cens <-
function (xleft, xright, probs = c(0.025, 0.5, 0.975), Alpha = 1,
    Kappa = 0, Gama = 0.4, distr.k = 1, distr.py0 = 1, distr.pz0 = 2,
    mu.pz0 = 3, sigma.pz0 = sqrt(10), delta = 4, kappa = 2, Delta = 2,
    Meps = 0.01, Nx = 150, Nit = 1500, Pbi = 0.1, epsilon = NULL,
    printtime = TRUE, extras = TRUE)
{
    if (is.null(distr.k))
        stop("Argument distr.k is NULL. Should be provided. See help for details.")
    if (is.null(distr.py0))
        stop("Argument distr.py0 is NULL. Should be provided. See help for details.")
    tInit <- proc.time()
    cens_data_check(xleft, xright)
    xpoint <- as.numeric(na.omit(0.5 * (xleft + xright)))
    npoint <- length(xpoint)
    censor_code <- censor_code_rl(xleft, xright)
    censor_code_filters <- lapply(0:3, FUN = function(x) censor_code ==
        x)
    names(censor_code_filters) <- 0:3
    n <- length(xleft)
    y <- seq(n)
    xsort = sort(xpoint)
    y[seq(n/2)] <- mean(xsort[seq(npoint/2)])
    y[-seq(n/2)] <- mean(xsort[-seq(npoint/2)])
    z <- rep(1, n)
    u <- 1
    if (is.null(epsilon))
        epsilon <- sd(xpoint)/4
    xx <- seq(min(xpoint) - epsilon, max(xpoint) + epsilon, length = Nx)
    Fxx <- matrix(NA, nrow = Nx, ncol = Nit)
    fx <- matrix(NA, nrow = n, ncol = Nit)
    R <- seq(Nit)
    U <- seq(Nit)
    Nmt <- seq(Nit)
    Allocs <- vector(mode = "list", length = Nit)
    if (extras) {
        means <- vector(mode = "list", length = Nit)
        sigmas <- vector(mode = "list", length = Nit)
        weights <- vector(mode = "list", length = Nit)
        Js <- vector(mode = "list", length = Nit)
    }
    mu.py0 = mean(xpoint)
    sigma.py0 = sd(xpoint)
    for (j in seq(Nit)) {
        if (floor(j/500) == ceiling(j/500))
            cat("MCMC iteration", j, "of", Nit, "\n")
        tt <- comp2(y, z)
        ystar <- tt$ystar
        zstar <- tt$zstar
        nstar <- tt$nstar
        rstar <- tt$rstar
        idx <- tt$idx
        Allocs[[max(1, j - 1)]] <- idx
        if (Gama != 0)
            u <- gs3(u, n = n, r = rstar, alpha = Alpha, beta = Kappa,
                gama = Gama, delta = Delta)
        JiC <- MvInv(eps = Meps, u = u, alpha = Alpha, beta = Kappa,
            gama = Gama, N = 50001)
        Nm <- length(JiC)
        TauyC <- rk(Nm, distr = distr.py0, mu = mu.py0, sigma = sigma.py0)
        TauzC <- rk(Nm, distr = distr.pz0, mu = mu.pz0, sigma = sigma.pz0)
        if (distr.pz0 == 2) {
            tt <- gsYZstarcens2(ystar = ystar, zstar = zstar,
                nstar = nstar, rstar = rstar, idx = idx, xleft = xleft,
                xright = xright, censor_code = censor_code, delta = delta,
                kappa = kappa, distr.k = distr.k, distr.py0 = distr.py0,
                mu.py0 = mu.py0, sigma.py0 = sigma.py0, distr.pz0 = distr.pz0,
                mu.pz0 = mu.pz0, sigma.pz0 = sigma.pz0)
            ystar <- tt$ystar
            zstar <- tt$zstar
        }
        tt <- gsHP(ystar, rstar, distr.py0)
        mu.py0 <- tt$mu.py0
        sigma.py0 <- tt$sigma.py0
        Jstar <- rgamma(rstar, nstar - Gama, Kappa + u)
        Tauy <- c(TauyC, ystar)
        Tauz <- c(TauzC, zstar)
        J <- c(JiC, Jstar)
        tt <- fcondYZXAcens2(xleft = xleft, xright = xright,
            censor_code_filters = censor_code_filters, distr = distr.k,
            Tauy = Tauy, Tauz = Tauz, J = J)
        y <- tt[, 1]
        z <- tt[, 2]
        Fxx[, j] <- fcondXA2(xx, distr = distr.k, Tauy, Tauz,
            J)
        fx[, j] <- fcondXA2cens2(xleft = xleft, xright = xright,
            censor_code_filters = censor_code_filters, distr = distr.k,
            Tauy, Tauz, J)
        R[j] <- rstar
        U[j] <- u
        Nmt[j] <- Nm
        if (extras) {
            means[[j]] <- Tauy
            sigmas[[j]] <- Tauz
            weights[[j]] <- J/sum(J)
            Js[[j]] <- J
        }
    }
    tt <- comp2(y, z)
    Allocs[[Nit]] <- tt$idx
    biseq <- seq(floor(Pbi * Nit))
    Fxx <- Fxx[, -biseq]
    qx <- as.data.frame(t(apply(Fxx, 1, quantile, probs = probs)))
    names(qx) <- paste("q", probs, sep = "")
    qx <- cbind(mean = apply(Fxx, 1, mean), qx)
    R <- R[-biseq]
    U <- U[-biseq]
    Allocs <- Allocs[-biseq]
    if (extras) {
        means <- means[-biseq]
        sigmas <- sigmas[-biseq]
        weights <- weights[-biseq]
        Js <- Js[-biseq]
    }
    cpo <- 1/apply(1/fx[, -biseq], 1, mean)
    if (printtime) {
        cat(" >>> Total processing time (sec.):\n")
        print(procTime <- proc.time() - tInit)
    }
    res = list(xx = xx, qx = qx, cpo = cpo, R = R, U = U,
               Allocs = Allocs, Nm = Nmt, Nx = Nx, Nit = Nit, Pbi = Pbi,
               procTime = procTime, distr.k = distr.k, data = data.frame(left = xleft, right = xright))
    if (extras) {
      res$means = means
      res$sigmas = sigmas
      res$weights = weights
      res$Js = Js
    }
    return(res)
}
