#' Function extracts confidence intervals for co-occurrence values.
#'
#' @title Derive parameter confidence intervals
#' @param draws Set of mcmc draws.
#' @param species_names Names of species to analyse.
#' @param quantile_n 0.025, "mean" or 0.975. Default of "mean".
#' 
#' @return An n x n table of interval values for our species co-occurrence parameters.
#'
#' @export


get_cooc_ints <- function(draws,species_names,quantile_n="mean") {
  n_species <- length(species_names)
  initial_matrix <- as.matrix(calculate(R,draws))
  if (quantile_n=="mean") { result_matrix <- apply(initial_matrix, 2, mean)
  } else {
    result_matrix <- apply(initial_matrix, 2, quantile,
                           probs = quantile_n)
  }
  cooc <- matrix(result_matrix,nrow = n_species, ncol=n_species)
  rownames(cooc) <- species_names
  colnames(cooc) <- species_names
  return(cooc)
}

