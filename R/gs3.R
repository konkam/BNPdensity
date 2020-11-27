#' Conditional posterior distribution of latent U
#'
#' This function simulates from the conditional posterior distribution of the
#' latent U.
#'
#' For internal use.
#'
#' @keywords internal
#' @examples
#'
#' ## The function is currently defined as
#' function(ut, n = 200, r = 20, alpha = 1, beta = 1, gama = 1 / 2,
#'          delta = 2) {
#'   w <- ut
#'   ratio <- NaN
#'   while (is.nan(ratio)) {
#'     v <- ustar <- rgamma(1, shape = delta, rate = delta / ut)
#'     vw <- v / w
#'     vb <- v + beta
#'     wb <- w + beta
#'     A <- vw^(n - 2 * delta)
#'     B <- (vb / wb)^(r * gama - n)
#'     D <- vb^gama - wb^gama
#'     E <- 1 / vw - vw
#'     ratio <- A * B * exp(-alpha / gama * D - delta * E)
#'   }
#'   p <- min(1, ratio)
#'   u <- ifelse(runif(1) <= p, ustar, ut)
#'   return(u)
#' }
gs3 <-
  function(ut, n, r, alpha, beta, gama, delta) {
    w <- ut
    ratio <- NaN
    while (is.nan(ratio)) {
      v <- ustar <- rgamma(1, shape = delta, rate = delta / ut)
      vw <- v / w
      vb <- v + beta
      wb <- w + beta
      A <- vw^(n - 2 * delta)
      B <- (vb / wb)^(r * gama - n)
      D <- vb^gama - wb^gama
      E <- 1 / vw - vw
      ratio <- A * B * exp(-alpha / gama * D - delta * E)
    }
    p <- min(1, ratio)
    u <- ifelse(runif(1) <= p, ustar, ut)
    return(u)
  }

gs3_adaptive <- function(ut, n, r, alpha, beta, gama, delta, U, iter, adapt = FALSE) {
  target_acc_rate <- 0.44
  batch_size <- 100
  if (adapt && (iter %% batch_size == 0)) {
    acc_rate <- length(unique(U[(iter - batch_size + 1):iter])) / batch_size
    logincrement <- 2*min(0.25, 1 / sqrt(iter))
    # increment = min(0.5, 5 / sqrt(iter))
    if (acc_rate < 0.44) {
      delta_i <- delta * exp(logincrement)
    }
    else {
      delta_i <- delta * exp(-logincrement)
    }
  }
  else {
    delta_i <- delta
  }
  # print(delta_i)
  u_prime <- gs3(ut, n, r, alpha, beta, gama, delta_i)
  return(list(u_prime = u_prime, delta = delta_i))
}
