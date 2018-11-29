#' Function extracts confidence intervals for our parameters.
#'
#' @title Derive parameter confidence intervals
#' @param parameter Parameter for which we want to derive credible intervals.
#' @param draws Set of mcmc draws.
#' @param species_names Names of species to analyse.
#' @param env_names Names of environmental variables.
#' @param ssv Whether or not we are trying to extract species-specific variables.
#' @param quantile_n 0.025, "mean" or 0.975. Default of "mean".
#' 
#' @return A table of interval values for our parameters.
#'
#' @export


get_param_ints <- function(parameter,draws,species_names,env_names,ssv=FALSE,quantile_n="mean") {
  n_species <- length(species_names)
  n_env_shared <- length(env_names)
  initial_matrix <- as.matrix(calculate(parameter,draws))
  if (quantile_n=="mean") { result_matrix <- apply(initial_matrix, 2, mean)
  } else {
    result_matrix <- apply(initial_matrix, 2, quantile,
                           probs = quantile_n)
  }
  if (ssv==TRUE) {
    ssv_matrix <- matrix(result_matrix,nrow = n_species, ncol=n_species)
    betas <- matrix(diag(ssv_matrix),nrow = n_species,ncol=1)
    colnames(betas) <- "ssv"
  } else {
    betas <- matrix(result_matrix,nrow = n_species, ncol=n_env_shared)
    colnames(betas) <- env_names
  }
  return(betas)
}

