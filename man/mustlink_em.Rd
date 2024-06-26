% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mustlink_em.R
\name{mustlink_em}
\alias{mustlink_em}
\title{Must-Link / Positive Constraint EM GMM.}
\usage{
mustlink_em(
  data,
  constraints_unique,
  params,
  clust_num,
  zone_num,
  burnin = 2,
  maxit = 10000,
  eps = 1e-10,
  print_freq = 1,
  drop_cluster = FALSE
)
}
\arguments{
\item{data}{Dataset in \code{matrix} or \code{data.frame} format.}

\item{constraints_unique}{vector of constrained set labels,
(unconstrained events have unique labels).}

\item{params}{Model parameters, e.g., output from \code{initial_parameters}.}

\item{clust_num}{Number of clusters / components.}

\item{zone_num}{Number of zones.}

\item{burnin}{Controls how many loops are completed before testing for
likelihood convergence.}

\item{maxit}{Maximum number of EM iterations.}

\item{eps}{Convergence criterion for relative difference in log-likelihood.}

\item{print_freq}{Controls how frequently the log-likelihood and time are
printed in the EM loop.}

\item{drop_cluster}{Whether empty clusters should be dropped.}
}
\value{
List:
\itemize{
\item postprob_block:
\item params:
\item loglike: vector of log-likelihood values for each iteration
}
}
\description{
EM algorithm
}
