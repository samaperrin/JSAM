#' Calculates likelihood of species presences based solely on environmental values (marignal prediction).
#'
#' @title Marginal prediction of species likelihoods at sites.
#' @param draws A set of mcmc draws produced by the GRETA function draws.
#' @param X_newdata Data for new sites to predict for.
#' @param site_id_newdata Catchment IDs for new sites.
#' @param X_assoc_newdata Avergae association correlation matrix.
#' @param species_list List of species names.
#'
#' @return Matrix of site x species with likelihoods of presence at every site.
#'
#' @export



predict_species_perc <- function(draws,X_newdata, site_id_newdata = NULL, X_assoc_newdata = NULL,species_list) {

  eta_new <- pred_eta(X_newdata, site_id_newdata,X_assoc_newdata)
  p_new <- ilogit(eta_new)

  fill_preds <- matrix(NA,nrow=nrow(p_new),ncol=ncol(p_new))

  for (i in 1:length(species_list)) {
    p_new1_draws1 <- calculate(p_new[, i], draws)
    p_new1_mn1 <- c(colMeans(as.matrix(p_new1_draws1)))
    fill_preds[,i] <- p_new1_mn1
  }
  colnames(fill_preds) <- species_list
  return(fill_preds)
}
