#' Creates sppecies cooccurrence matrix at certain value of environmental parameter.
#'
#' @title Predictive species association matrix
#' @param darws Series of mcmc draws produced by the GRETA function draws.
#' @param X_assoc_pred Vector of previously defined values, which represent standard deviations to or from the mean of the environmental parameter you are monitoring changes over.
#' @param species_names Names of species in model.
#'
#' @return A species x species matrix with assocation means for species associations when an environmental parameter is of a certain value.
#'
#' @export


pred_correlation <- function(draws,X_assoc_pred,species_names) {
  n_species <- length(species_names)
  returned <- list()
  for (i in 1:length(X_assoc_pred)) {
    temp_R <- colMeans(as.matrix(calculate(R_pred[[i]],draws)))
    temp_R_matrix <- matrix(temp_R,n_species,n_species)
    colnames(temp_R_matrix) <- species_names
    rownames(temp_R_matrix) <- species_names
    diag(temp_R_matrix) <- 0
    returned[[i]] <- temp_R_matrix
  }
  return(returned)
}
