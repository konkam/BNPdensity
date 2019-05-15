MixNRMI1cens <-
function (xleft, xright, probs = c(0.025, 0.5, 0.975), Alpha = 1,
    Beta = 0, Gama = 0.4, distr.k = 1, distr.p0 = 1, asigma = 0.5,
    bsigma = 0.5, delta = 3, Delta = 2, Meps = 0.01, Nx = 150,
    Nit = 1500, Pbi = 0.1, epsilon = NULL, printtime = TRUE,
    extras = FALSE)
{
    if (is.null(distr.k))
        stop("Argument distr.k is NULL. Should be provided. See help for details.")
    if (is.null(distr.p0))
        stop("Argument distr.p0 is NULL. Should be provided. See help for details.")
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
    u <- 1
    sigma <- sd(xpoint)
    if (is.null(epsilon))
        epsilon <- sd(xpoint)/4
    xx <- seq(min(xpoint) - epsilon, max(xpoint) + epsilon, length = Nx)
    Fxx <- matrix(NA, nrow = Nx, ncol = Nit)
    fx <- matrix(NA, nrow = n, ncol = Nit)
    R <- seq(Nit)
    S <- seq(Nit)
    U <- seq(Nit)
    Nmt <- seq(Nit)
    Allocs <- vector(mode = "list", length = Nit)
    if (extras) {
        means <- vector(mode = "list", length = Nit)
        weights <- vector(mode = "list", length = Nit)
        Js <- vector(mode = "list", length = Nit)
    }
    mu.p0 = mean(xpoint)
    sigma.p0 = sd(xpoint)
    for (j in seq(Nit)) {
        if (floor(j/500) == ceiling(j/500))
            cat("MCMC iteration", j, "of", Nit, "\n")
        tt <- comp1(y)
        ystar <- tt$ystar
        nstar <- tt$nstar
        r <- tt$r
        idx <- tt$idx
        Allocs[[max(1, j - 1)]] <- idx
        if (Gama != 0)
            u <- gs3(u, n = n, r = r, alpha = Alpha, beta = Beta,
                gama = Gama, delta = Delta)
        JiC <- MvInv(eps = Meps, u = u, alpha = Alpha, beta = Beta,
            gama = Gama, N = 50001)
        Nm <- length(JiC)
        TauiC <- rk(Nm, distr = distr.p0, mu = mu.p0, sigma = sigma.p0)
        ystar <- gs4cens2(ystar = ystar, xleft = xleft, xright = xright,
            censor_code = censor_code, idx = idx, distr.k = distr.k,
            sigma.k = sigma, distr.p0 = distr.p0, mu.p0 = mu.p0,
            sigma.p0 = sigma.p0)
        Jstar <- rgamma(r, nstar - Gama, Beta + u)
        Tau <- c(TauiC, ystar)
        J <- c(JiC, Jstar)
        tt <- gsHP(ystar, r, distr.p0)
        mu.p0 <- tt$mu.py0
        sigma.p0 <- tt$sigma.py0
        y <- fcondYXAcens2(xleft = xleft, xright = xright, censor_code_filters = censor_code_filters,
            distr = distr.k, Tau = Tau, J = J, sigma = sigma)
        sigma <- gs5cens2(sigma = sigma, xleft = xleft, xright = xright,
            censor_code = censor_code, y = y, distr = distr.k,
            asigma = asigma, bsigma = bsigma, delta = delta)
        Fxx[, j] <- fcondXA(xx, distr = distr.k, Tau = Tau, J = J,
            sigma = sigma)
        fx[, j] <- fcondYXAcens2(xleft = xleft, xright = xright,
            censor_code_filters = censor_code_filters, distr = distr.k,
            Tau = Tau, J = J, sigma = sigma)
        R[j] <- r
        S[j] <- sigma
        U[j] <- u
        Nmt[j] <- Nm
        if (extras) {
            means[[j]] <- Tau
            weights[[j]] <- J/sum(J)
            Js[[j]] <- J
        }
    }
    tt <- comp1(y)
    Allocs[[Nit]] <- tt$idx
    biseq <- seq(floor(Pbi * Nit))
    Fxx <- Fxx[, -biseq]
    qx <- as.data.frame(t(apply(Fxx, 1, quantile, probs = probs)))
    names(qx) <- paste("q", probs, sep = "")
    qx <- cbind(mean = apply(Fxx, 1, mean), qx)
    R <- R[-biseq]
    S <- S[-biseq]
    U <- U[-biseq]
    Allocs <- Allocs[-biseq]
    if (extras) {
        means <- means[-biseq]
        weights <- weights[-biseq]
        Js <- Js[-biseq]
    }
    cpo <- 1/apply(1/fx[, -biseq], 1, mean)
    if (printtime) {
        cat(" >>> Total processing time (sec.):\n")
        print(procTime <- proc.time() - tInit)
    }
    res = list(xx = xx, qx = qx, cpo = cpo, R = R, S = S,
               U = U, Allocs = Allocs, Nm = Nmt, Nx = Nx, Nit = Nit,
               Pbi = Pbi, procTime = procTime, distr.k = distr.k, data_left = xleft, data_right = xright)
    if (extras) {
      res$means = means
      res$weights = weights
      res$Js = Js
    }
    return(res)
}
