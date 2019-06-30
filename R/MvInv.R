MvInv <-
function (eps, u = 0.5, alpha = 1, beta = 1, gama = 1/2, N = 3001) # eps no longer required
{
    x <- -log(seq(from = exp(-1e-05), to = exp(-10), length = N))
    f <- alpha/gamma(1 - gama) * x^(-(1 + gama)) * exp(-(u +
        beta) * x)
    dx <- diff(x)
    h <- (f[-1] + f[-N])/2
    Mv <- c(rev(cumsum(rev(dx[-N] * h[-N]))), 0)
    # err <- 1
    # w <- 0
    # v <- NULL
    # while (err > eps) {
    #     w <- w + rgamma(1, 1, 1)
    #     v <- c(v, x[which.min(Mv > w)])
    #     err <- min(v)/sum(v)
    # }
    M <- thresholdGG(alpha = alpha,
                     kappa = beta+u,
                     gama = gama) # upper bound defined via the grid
    W <- rexp(M)
    W <- cumsum(W)
    x_which_min = function(w){ # I guess that this function could be defined outside of MvInV
      x[which.min(Mv > w)]
    }
    x_which_min = Vectorize(x_which_min, vectorize.args = "w")
    v <- x_which_min(W)
    return(v)
}
