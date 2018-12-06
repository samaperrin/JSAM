#' Calculates conditional likelihood of one species presences based on environmental values and presence of other species at a range of sites.
#'
#' @title Conditional prediction of species at range of sites site.
#' @param draws A set of mcmc draws produced by the GRETA function draws.
#' @param betas_shared Non species-specific enrivonmental parameters.
#' @param betas_ssv Species-specific environmental parameters.
#' @param occupancy_matrix Occupancy matrix for all species.
#' @param temp_vector Vector of temperatures at different sites, used to determine cooccurrence matrices.
#' @param focal_species Integer which dictates which species in the occupancy vector you want to test likelihood of.
#' @param ssv Whether or not the model contains species specific variables.
#'
#' @return Vector of length site giving likelihoods of species presence at every site.
#'
#' @export
#'
#'
pred_cond_entire <- function(draws,betas_shared,betas_ssv=NULL,site_matrix,occupancy_matrix,temp_vector,focal_species,ssv=FALSE) {
  presences <- matrix(NA,ncol=1,nrow=nrow(site_matrix))
  for (i in 1:nrow(site_matrix)) {
    site_env <- matrix(site_matrix[i,],nrow=ncol(site_matrix),ncol=1)
    occupancy <- as.vector(occupancy_matrix[i,])
    cooccurrence <- pred_correlation(draws,temp_vector[i],colnames(occupancy_matrix))
    diag(cooccurrence[[1]]) <- 1

    presences[i] <- pred_cond(betas_shared,betas_ssv,site_env,occupancy,cooccurrence[[1]],focal_species,ssv)
    if (i %% 10 == 0) {print(paste0(i," runs are complete.")) }
  }
  return(presences)
}
