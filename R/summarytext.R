#' Common text for the summary S3 methods
#'
#' @param fit NRMIx object
#' @param kernel_comment Text specific to the parametric and nonparametric nature of the model
#' @param number_of_clusters Flag to decide whether to compute the optimal clustering
#'
#' @return Prints out the text for the summary S3 methods
summarytext = function(fit, kernel_comment, number_of_clusters = FALSE){
  NRMI_comment <- paste("Density estimation using a", comment_on_NRMI_type(fit$NRMI_params))
  data_comment <- paste("There were", length(fit$data), "data points.")
  MCMC_comment <- paste("The MCMC algorithm was run for ", fit$Nit, " iterations with ", 100 * fit$Pbi, "% discarded for burn-in.", sep = "")
  if (number_of_clusters) {
    estimated_clustering <- compute_optimal_clustering(fit)
    clustering_comment <- paste("The estimated number of clusters in the data is ", length(unique(estimated_clustering)), ".", sep = "")
  }
  else {
    clustering_comment <- "To obtain information on the estimated number of clusters, please use summary(object, number_of_clusters = TRUE)."
  }
  writeLines(paste(NRMI_comment, "\n", kernel_comment, "\n", data_comment, "\n", MCMC_comment, "\n", clustering_comment, sep = ""))

}
