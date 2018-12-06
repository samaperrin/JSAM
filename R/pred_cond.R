#' Calculates conditional likelihood of one species presences based on environmental values and presence of other species.
#'
#' @title Conditional prediction of species at one site.
#' @param betas_shared Non species-specific enrivonmental parameters.
#' @param betas_ssv Species-specific environmental parameters.
#' @param occupancy Occupancy vector for all species.
#' @param cooccurrence Cooccurrence matrix for the site (detrmined by the temperature at the site, use pred_correlation to derive this).
#' @param focal_species Integer which dictates which species in the occupancy vector you want to test likelihood of.
#' @param ssv Whether or not the model contains species specific variables.
#'
#' @return Value with likelihood of species at that site.
#'
#' @export

pred_cond <- function(betas_shared,betas_ssv=NULL,site_env,occupancy,cooccurrence,focal_species,ssv=FALSE){
  library(tmvtnorm)
  if (ssv==TRUE) {
    betas_ssv_matrix <- matrix(0,length(betas_ssv),length(betas_ssv))
    diag(betas_ssv_matrix) <- betas_ssv
    betas <- cbind(betas_shared,betas_ssv_matrix)
  } else {betas <- betas_shared}
  mean1 <- c(betas %*% site_env)
  sigma1 <- cooccurrence # Need to edit this so it generates cooccurence matrix based on our temeprature

  lower1 <- c(ifelse(occupancy==0,-Inf,0))
  lower1[focal_species] <- -Inf
  upper1 <- c(ifelse(occupancy==0,0,Inf))

  lowerx1 <- c(lower1)
  lowerx1[focal_species] <- 0
  upperx1 <- c(upper1)
  upperx1[focal_species] <- Inf

  result <- ptmvnorm(mean=mean1,sigma=sigma1,lower=lower1,upper=upper1,lowerx=lowerx1,upperx=upperx1)
  return(result[1])
}
