#' Title
#'
#' @param x
#' @param Alpha
#' @param Gama
#' @param Nx
#' @param Nit
#' @param Pbi
#' @param epsilon
#' @param printtime
#'
#' @return
#' @export
#'
#' @examples
MixPY1 <- function(x, Alpha = 1, Gama = 0.4, Nx = 150, Nit = 1500, Pbi = 0.5, epsilon = NULL, printtime = TRUE, extras = TRUE) {
  if (!requireNamespace("BNPmix", quietly = TRUE)) {
    stop("Package \"BNPmix\" needed for this function to work. Please install it.",
      call. = FALSE
    )
  }

  nburn <- floor(Pbi * Nit)
  if (is.null(epsilon)) {
    epsilon <- sd(x) / 4
  }
  xx <- seq(min(x) - epsilon, max(x) + epsilon, length = Nx)

  restmp <- BNPmix::PYdensity(
    y = x,
    mcmc = list(niter = Nit, nburn = nburn, model = "L", print_message = printtime),
    prior = list(strength = Alpha, discount = Gama),
    output = list(grid = xx, out_param = TRUE, out_type = "FULL")
  )

  res <- list(
    xx = xx, R = unlist(lapply(X = restmp$mean, FUN = length)), S = sqrt(restmp$sigma2)
    )
  if (extras) {
    res$means <- restmp$mean
    res$weights <- weights
    res$Js <- Js
  }
  return(structure(res, class = "PY1"))
}
