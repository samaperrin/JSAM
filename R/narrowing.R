#' Takes co-occurrence matrix and narrows down to only show relationships between species of interest.
#'
#' @title Subsetting of co-occurrence matrix.
#' @param correlation_matrix_list An N x N co-occurrence matrix.
#' @param focal_species_list A vector of M characters which correspond to species names.
#'
#' @return An M x M co-occurrence matrix.
#'
#' @export
#'
narrowing <- function(correlation_matrix_list, focal_species_list) {
  narrowed_correlation_matrix_list <- vector("list", length(correlation_matrix_list))
  for(i in 1:length(correlation_matrix_list)) {
    narrowed_correlation_matrix_list[[i]] <- correlation_matrix_list[[i]][focal_species_list,focal_species_list]
  }
  return(narrowed_correlation_matrix_list)
}
