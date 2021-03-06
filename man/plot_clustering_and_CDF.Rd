% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/posterior_clustering_analysis.R
\name{plot_clustering_and_CDF}
\alias{plot_clustering_and_CDF}
\title{Plot the clustering and the Cumulative Distribution Function}
\usage{
plot_clustering_and_CDF(fit, clustering, label_vector = NULL)
}
\arguments{
\item{fit}{The fitted object, obtained from one of the MixNRMIx functions}

\item{clustering}{A vector of integers with the same length as the data,
representing the allocation variable for data each point.}

\item{label_vector}{A vector of data labels to be plotted, to provide some
identification to each point.}
}
\value{
A plot of the Cumulative Distribution Function (or Turnbull estimate
for censored data) with data points whose color denotes the cluster
allocation. For censored data, right or left censored data points are not
represented, while interval censored data points are represented at the
middle of the censoring interval.
}
\description{
This is a function to visualize the clustering induced by the BNP model. The
data points are plotted with a color reflecting their cluster.
}
