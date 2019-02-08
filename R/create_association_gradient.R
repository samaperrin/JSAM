#' Takes list of association matrices and creates a new matrix showing the change in species associations over a change in environmental variable for one species.
#'
#' @title Species associations gradient.
#' @param association_matrix_list A list of N x N co-occurrence matrices, one for each level on the environmental variable that was emasured.
#' @param focal_species The single species you want to observe the changes in associations for
#'
#' @return A species x covariate level matrix.
#'
#' @export
#'
create_association_gradient <- function(association_matrix_list,focal_species) {
  ascmat_length <- length(association_matrix_list)
  n_species <- nrow(association_matrix_list[[1]])
  asc_map <- matrix(NA,n_species,ascmat_length)
  for(i in 1:ascmat_length) {
    asc_ind <- association_matrix_list[[i]][focal_species,]
    asc_map[,i] <- asc_ind
  }
  rownames(asc_map) <- colnames(association_matrix_list[[1]])
  return(asc_map)
}
