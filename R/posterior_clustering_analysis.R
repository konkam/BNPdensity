compute_optimal_clustering <- function(fit, silent = FALSE, loss_function = "VI") {
  if (!requireNamespace("mcclust.ext", quietly = TRUE)) {
    stop("Package mcclust.ext is needed for this function to work. Please install it.",
      call. = FALSE
    )
  }
  if (!requireNamespace("mcclust", quietly = TRUE)) {
    stop("Package mcclust is needed for this function to work. Please install it.",
      call. = FALSE
    )
  }

  if (!silent) print("Estimating the optimal clustering minimising the expected VI loss is time consuming, please have a little patience. For a faster result, you may use Binder's loss function instead.")

  fit.draw <- Reduce(rbind, fit$Allocs)
  psm <- mcclust::comp.psm(fit.draw)
  fit_VI <- mcclust.ext::minVI(psm, fit.draw, method = ("greedy"))
  return(fit_VI$cl)
}

# clustering = compute_optimal_clustering(out)

plot_clustering_and_CDF_noncensored <- function(fit, clustering, label_vector = NULL) {
  data <- fit$data

  grid <- grid_from_data(data)

  if (is_semiparametric(fit)) {
    cdf <- get_CDF_semi_BNPdensity(fit = fit, xs = data)
  }
  else {
    cdf <- get_CDF_full_BNPdensity(fit = fit, xs = data)
  }
  p <- ggplot(data.frame(data = data, cluster_id = clustering, cdf = cdf)) +
    theme_bw() +
    geom_step(aes(x = data, y = ecdf(data)(data))) +
    geom_point(aes(x = data, y = cdf, colour = factor(cluster_id))) +
    viridis::scale_colour_viridis(discrete = TRUE) +
    theme(legend.position = "none") +
    ylab("CDF") +
    xlab("Data")

  if (!is.null(label_vector)) {
    p + geom_text(data = data.frame(txt = label_vector, x = data, y = cdf + 0.05, cluster_id = clustering), aes(x = x, y = y, colour = factor(cluster_id), label = txt))
  }
  else {
    return(p)
  }
}

decide_abscissa <- function(censored_data, clustering) {
  df <- cbind(censored_data, data.frame(
    cluster_id = clustering,
    loc = rowMeans(censored_data),
    is_censored = is.na(rowMeans(censored_data))
  ))
  df$loc <- unlist(mapply(FUN = function(is_cens, cluster_id, loc) {
    if (is_cens) {
      mean(df$loc[df$cluster_id == cluster_id], na.rm = TRUE)
    }
    else {
      loc
    }
  }, df$is_censored, df$cluster_id, df$loc, SIMPLIFY = FALSE))
  return(df)
}

plot_clustering_and_CDF_censored <- function(fit, clustering, label_vector = NULL) {
  data <- fit$data

  # grid <- grid_from_data(data)
  grid <- decide_abscissa(data, clustering)$loc

  Survival_object <- survival::survfit(formula = survival::Surv(data$left, data$right, type = "interval2") ~ 1)

  if (is_semiparametric(fit)) {
    cdf <- get_CDF_semi_BNPdensity(fit = fit, xs = grid[!is.na(grid)])
  }
  else {
    cdf <- get_CDF_full_BNPdensity(fit = fit, xs = grid[!is.na(grid)])
  }
  ggplot2::ggplot(data = data.frame(data = grid[!is.na(grid)], CDF = cdf, cluster_id = clustering[!is.na(grid)]), aes(x = data, y = CDF)) +
    geom_point(aes(colour = factor(cluster_id))) +
    theme_classic() +
    geom_step(data = data.frame(x = c(Survival_object$time, max(grid)), y = c(1 - Survival_object$surv, 1)), aes(x = x, y = y)) +
    viridis::scale_colour_viridis(discrete = TRUE) +
    theme(legend.position = "none") +
    ylab("CDF") +
    xlab("Data")

  if (!is.null(label_vector)) {
    p + geom_text(data = data.frame(txt = label_vector[!is.na(grid)], x = grid[!is.na(grid)], y = cdf + 0.05, cluster_id = clustering[!is.na(grid)]), aes(x = x, y = y, colour = factor(cluster_id), label = txt))
  }
  else {
    return(p)
  }
}
