MvInv <-
  function(eps, u = 0.5, alpha = 1, beta = 1, gama = 1 / 2, N = 3001) # eps no longer required
  {
    x <- -log(seq(from = exp(-1e-05), to = exp(-10), length = N))
    f <- alpha / gamma(1 - gama) * x^(-(1 + gama)) * exp(-(u +
      beta) * x)
    dx <- diff(x)
    h <- (f[-1] + f[-N]) / 2
    Mv <- c(rev(cumsum(rev(dx[-N] * h[-N]))), 0)

    M <- ceiling(thresholdGG(
      alpha = alpha,
      kappa = beta + u,
      gama = gama
    )) # upper bound defined via the grid
    M <- max(10, M) # We wish to make sure we at least use a few jumps
    W <- rexp(n = M)
    W <- cumsum(W)
    # x_which_min = function(w){ # I guess that this function could be defined outside of MvInV
    #   x[which.min(Mv > w)]
    # }
    # x_which_min = Vectorize(x_which_min, vectorize.args = "w")
    # v <- x_which_min(W)
    if (M < 25) {
      ## This version is faster because it has no loop, but it involves many passes over Mv which has 3001 elements
      return(fill_v1(M, Mv, W, x))
    }
    else {
      return(fill_v2(M, Mv, W, N, x)) # Faster for large values of N
    }
  }

fill_v1 <- function(M, Mv, W, x) {
  v <- rep(NA, M)
  for (j in seq(M)) v[j] <- x[which.min(Mv > W[j])]
  return(v)
}

fill_v2 <- function(M, Mv, W, N, x) {
  v <- rep(NA, M)
  iMv <- N
  for (i in seq(M)) {
    while (iMv > 0 && Mv[iMv] < W[i]) {
      iMv <- iMv - 1
      # print(paste(iMv, Mv[iMv], W[i]))
    }
    v[i] <- x[iMv + 1] # This index shift is to keep consistency  with previous version of the function, not necessary.
  }
  # for (j in seq(M)) v[j] <- x[which.min(Mv > W[j])]
  return(v)
}

fill_v3 <- function(M, Mv, W, N, x) {
  v <- rep(NA, M)
  i <- 1
  for (iMv in N:1) {
    if (i > M) break()
    if (Mv[iMv] > W[i]) {
      v[i] <- x[iMv + 1] # This index shift is to keep consistency  with previous version of the function, not necessary.
      i <- i + 1
    }
  }
  if (i < M) v[i:M] <- x[1]
  # for (j in seq(M)) v[j] <- x[which.min(Mv > W[j])]
  return(v)
}

# fill_v4 = function(M, Mv, W, N, x){
#   v <- rep(NA, M)
#   iMv = N
#   for (i in seq(M)){
#     while (Mv[iMv] < W[i]){
#       iMv = iMv-1
#     }
#     print(paste("now index", i))
#     print(which.min(Mv > W[i]))
#     print(iMv+1)
#     v[i] = x[iMv+1]
#   }
#   # for (j in seq(M)) v[j] <- x[which.min(Mv > W[j])]
#   return(v)
# }
#

# library(microbenchmark)
#
# microbenchmark(fill_v1(M, Mv, W, x), fill_v2(M, Mv, W, N, x), fill_v3(M, Mv, W, N, x))
