% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_cooc_ints.R
\name{get_cooc_ints}
\alias{get_cooc_ints}
\title{Derive parameter confidence intervals}
\usage{
get_cooc_ints(draws, species_names, quantile_n = "mean")
}
\arguments{
\item{draws}{Set of mcmc draws.}

\item{species_names}{Names of species to analyse.}

\item{quantile_n}{0.025, "mean" or 0.975. Default of "mean".}
}
\value{
An n x n table of interval values for our species co-occurrence parameters.
}
\description{
Function extracts confidence intervals for co-occurrence values.
}
