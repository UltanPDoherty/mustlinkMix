% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mustlink_estep.R
\name{mustlink_estep}
\alias{mustlink_estep}
\title{E-Step for Must-link Constrained GMM.}
\usage{
mustlink_estep(
  data,
  constraints_info,
  params,
  event_num = nrow(data),
  var_num = ncol(params$mu),
  clust_num = nrow(params$mu)
)
}
\arguments{
\item{data}{Dataset in \code{matrix} or \code{data.frame} format.}

\item{constraints_info}{.}

\item{params}{Model parameters, e.g., output from \code{initial_parameters}.}

\item{event_num}{Number of observations in the data set.}

\item{var_num}{Number of variables in the data set.}

\item{clust_num}{Number of clusters / components.}
}
\value{
List:
\itemize{
\item loglike: log likelihood value
\item postprob_sets: posterior probability matrix
\item postprob_events: posterior probability matrix
}
}
\description{
Implement E-step of EM algorithm for GMM with positive / must-link
constraints.
}
