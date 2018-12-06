#' Creates matrix which shows beta effect confidence intervals for each parameter on each species.
#'
#' @title Creates CI matrix.
#' @param draws A set of mcmc draws produced by the GRETA function draws.
#' @param beta_shared GRETA array of environmental parameters x species.
#' @param beta_ssv GRETA array of species x species with species-speecific parameters defined along the diagonal. Not necessary if ssv=FALSE.
#' @param species_names Names of species in data.
#' @param env_names Names of environmental parameters.
#' @param ssv Whether or not you are calculating species specific parameter values.
#'
#' @return A matrix with lower, mean and upper confidence intervals for beta effects of parmaeters on species presence/absence.
#'
#' @export


get_beta_list <- function(draws,beta_shared,beta_ssv,species_names,env_names,ssv=FALSE){
  demo_betas_upper <- as.data.frame(get_param_ints(beta_shared,draws,species_names,env_names,quantile_n = 0.975))
  demo_betas_lower <- as.data.frame(get_param_ints(beta_shared,draws,species_names,env_names,quantile_n = 0.025))
  demo_betas_mean <- as.data.frame(get_param_ints(beta_shared,draws,species_names,env_names))
  if (ssv==TRUE) {
    demo_betas_ssv_upper <- get_param_ints(beta_ssv,draws,species_names,env_names,ssv=TRUE,quantile_n = 0.975)
    demo_betas_upper <- as.data.frame(cbind(demo_betas_upper,demo_betas_ssv_upper))
    demo_betas_ssv_lower <- get_param_ints(beta_ssv,draws,species_names,env_names,ssv=TRUE,quantile_n = 0.025)
    demo_betas_lower <- as.data.frame(cbind(demo_betas_lower,demo_betas_ssv_lower))
    demo_betas_ssv_mean <- get_param_ints(beta_ssv,draws,species_names,env_names,ssv=TRUE)
    demo_betas_mean <- as.data.frame(cbind(demo_betas_mean,demo_betas_ssv_mean))
  }
  demo_betas_mean$species <- species_names
  demo_betas_lower$species <- species_names
  demo_betas_upper$species <- species_names
  library(reshape2)
  upper_betas_reshaped <- melt(demo_betas_upper, id.vars = "species")
  lower_betas_reshaped <- melt(demo_betas_lower, id.vars = "species")
  mean_betas_reshaped <- melt(demo_betas_mean, id.vars = "species")
  full_betas1 <- merge(lower_betas_reshaped,mean_betas_reshaped,by=c("species","variable"))
  full_betas <- merge(full_betas1,upper_betas_reshaped,by=c("species","variable"))
  names(full_betas) <- c("species","variable","lower","mean","upper")
  return(full_betas)
}
