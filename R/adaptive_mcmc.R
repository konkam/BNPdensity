logf_u_cond_y <- function(u, n, r, gamma, kappa, a) {
  (n - 1) * log(u) + (r * gamma - n) * log(u + kappa) - a / gamma * (u + kappa)^gamma
}

logdprop_u <- function(u_prime, u, delta) {
  dgamma(x = u_prime, shape = delta, rate = delta / u, log = T)
}

rprop_u <- function(u, delta) {
  rgamma(n = 1, shape = delta, rate = delta / u)
}

logacceptance_ratio_u <- function(u, u_prime, n, r, gamma, kappa, a, delta) {
  log_ratio <- logf_u_cond_y(u_prime, n, r, gamma, kappa, a) - logf_u_cond_y(u, n, r, gamma, kappa, a) + logdprop_u(u, u_prime, delta) - logdprop_u(u_prime, u, delta)
  return(min(0, log_ratio))
}

#' Title
#'
#' @param n_iter
#' @param u0
#' @param n
#' @param r
#' @param gamma
#' @param kappa
#' @param a
#' @param delta
#'
#' @return
#' @export
#'
#' @examples
non_adaptive_mcmc_sampler <- function(n_iter, u0, n, r, gamma, kappa, a, delta) {
  res <- rep(NA, n_iter)
  i <- 0
  res[1] <- u0
  while (i < n_iter) {
    i <- i + 1
    u_prime <- rprop_u(u = res[i], delta = delta)
    logq1 <- logacceptance_ratio_u(u = res[i], u_prime = u_prime, n = n, r = r, gamma = gamma, kappa = kappa, a = a, delta = delta)

    if (log(runif(n = 1)) < logq1) {
      res[i + 1] <- u_prime
    }
    else {
      res[i + 1] <- res[i]
    }
  }
  return(res)
}

#' Title
#'
#' @param n_iter
#' @param u0
#' @param n
#' @param r
#' @param gamma
#' @param kappa
#' @param a
#' @param delta
#' @param debug
#'
#' @return
#' @export
#'
#' @examples
#' library(tidyverse)
#' out <- adaptive_mcmc_sampler(n_iter = 2^17, u0 = 1.2, n = 100, r = 10, gamma = 0.5, kappa = 1.2, a = 1., delta = 2, debug = T)
#' tibble(acc_rate = out$acc_rates, delta = out$deltas[2:length(out$deltas)]) %>%
#'   mutate(batch_number = seq_along(acc_rate)) %>%
#'   gather(variable, value, -batch_number) %>%
#'   ggplot(aes(x = batch_number, y = value)) +
#'   theme_bw() +
#'   facet_wrap(~variable, scales = "free_y") +
#'   geom_line()
adaptive_mcmc_sampler <- function(n_iter, u0, n, r, gamma, kappa, a, delta, debug = FALSE) {
  batch_size <- 100
  res <- rep(NA, n_iter)
  i <- 0
  res[1] <- u0
  batch_index <- 0
  delta_i <- delta

  if (debug) {
    deltas <- rep(NA, floor(n_iter / batch_size))
    deltas[1] <- delta
    deltas_idx <- 1
    acc_rates <- rep(NA, floor(n_iter / batch_size))
  }

  while (i < n_iter) {
    i <- i + 1

    if (i %% batch_size == 0) {
      batch_index <- batch_index + 1

      acc_rate <- length(unique(res[(i - batch_size + 1):i])) / batch_size
      if (acc_rate < 0.44) {
        delta_i <- delta_i * exp(2*min(0.01, 1 / sqrt(i)))
      }
      else {
        delta_i <- delta_i * exp(-2*min(0.01, 1 / sqrt(i)))
      }

      if (debug) {
        deltas_idx <- deltas_idx + 1
        deltas[deltas_idx] <- delta_i
        acc_rates[deltas_idx - 1] <- acc_rate
      }
    }

    u_prime <- rprop_u(u = res[i], delta = delta_i)
    logq1 <- logacceptance_ratio_u(u = res[i], u_prime = u_prime, n = n, r = r, gamma = gamma, kappa = kappa, a = a, delta = delta)

    if (log(runif(n = 1)) < logq1) {
      res[i + 1] <- u_prime
    }
    else {
      res[i + 1] <- res[i]
    }
  }
  if (debug) {
    return(list(res = res, deltas = deltas, acc_rates = acc_rates))
  }
  else {
    return(res)
  }
}
