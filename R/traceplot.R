#' Draw a traceplot for multiple chains
#'
#' @inherit convert_to_mcmc
#'
#' @returna A traceplot for multiple chains.
#' @export
#'
#' @examples
traceplot = function(fitlist){
  mcmc_object = convert_to_mcmc(fitlist)
  to_plot = lapply(seq_along(mcmc_object), function(chain_id){dplyr::mutate(data.frame(mcmc_object[[chain_id]]), chain_id = chain_id) %>%
      rowid_to_column("iteration")}) %>%
    dplyr::bind_rows() %>%
    tidyr::gather(param, value, -chain_id, -iteration)

  ggplot(to_plot, aes(x = iteration, y = value, colour = factor(chain_id), group = chain_id)) +
    geom_line() +
    facet_wrap(~param, scales = 'free') +
    theme_classic() +
    ylab("") +
    theme(legend.position = 'none') +
    xlab("Iteration")
}