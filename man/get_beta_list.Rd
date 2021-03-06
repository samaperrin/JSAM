% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_beta_list.R
\name{get_beta_list}
\alias{get_beta_list}
\title{Creates CI matrix.}
\usage{
get_beta_list(draws, beta_shared, beta_ssv, species_names, env_names,
  ssv = FALSE)
}
\arguments{
\item{draws}{A set of mcmc draws produced by the GRETA function draws.}

\item{beta_shared}{GRETA array of environmental parameters x species.}

\item{beta_ssv}{GRETA array of species x species with species-speecific parameters defined along the diagonal. Not necessary if ssv=FALSE.}

\item{species_names}{Names of species in data.}

\item{env_names}{Names of environmental parameters.}

\item{ssv}{Whether or not you are calculating species specific parameter values.}
}
\value{
A matrix with lower, mean and upper confidence intervals for beta effects of parmaeters on species presence/absence.
}
\description{
Creates matrix which shows beta effect confidence intervals for each parameter on each species.
}
