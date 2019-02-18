#    ---------------------------------------    #
#### 1. Package loading and data preparation ####
#    ---------------------------------------    #

library(devtools)
#install_github("samaperrin/JSAM")
library(JSAM)
library(greta)
library(truncnorm)
library(corrplot)
library(purrr)



# The following script contains several functions which are used to construct the model below.
source("./pred_g_ftn.R")

# Next we set the seed.
set.seed(123)

# We define our simulated dimensions (Note: I have formulated this for a system without a spatial hierarchy)
n_sites_sim <- 3200
n_env_sim <- 5
n_latent_sim <- 3
n_species_sim <- 20

# Simulate environmental variables
X_sim <- matrix(rnorm(n_sites_sim * n_env_sim), nrow = n_sites_sim)

# Define our temperature variable and the temperature variations that we want to observe association matrices for
X_assoc_sim <- X_sim[,1]
X_assoc_pred_sim <- c(-2,-1,0,1,2)

# Let's simulate some random species names
sim_species_names <- rep(NA, n_species_sim)
for(i in 1:n_species_sim) {
  sim_species_names[i] <- paste0("species",i)
}

# We now simulate random alpha and beta values. I've reduced the alpha values slightly to ensure lower occurences, hopefully mirroring a more natural environment.
alpha_sim <- rnorm(n_species_sim, -0.4, 1)

beta_sim <- matrix(rnorm(n_env_sim*n_species_sim,0,1), n_env_sim, n_species_sim)
# Lambda * z will produce our error matrix
z_sim <- matrix(rnorm(n_latent_sim*n_sites_sim,0,1), n_latent_sim, n_sites_sim)

# The following is a small function used to create a lambda in which the upper diagonal is comprised of zeroes and the diagonal is positive. It utilises the package truncnorm, loaded above.
create_sim_lambda <- function(n_species,n_latent) {
  lambda_sim <- matrix(0, n_species, n_latent)
  diag(lambda_sim) <- rtruncnorm(n_latent, 0, Inf, 0, 1)
  lower_idx <- lower.tri(lambda_sim)
  lambda_sim[lower_idx] <- rnorm(sum(lower_idx),0,1)
  return(lambda_sim)
}
lambda_sim_int <- create_sim_lambda(n_species_sim,n_latent_sim)
lambda_sim_coef <- create_sim_lambda(n_species_sim,n_latent_sim)

# Now produce our prediction function for our later correlation matrix
R_sim_pred <- pred_R(X_assoc_pred_sim,lambda_int=lambda_sim_int, lambda_coef=lambda_sim_coef)
corrplot(R_sim_pred[[3]], method = "color",type="lower",
         col = colorRampPalette(c("blue", "white", "red"))(200),order="FPC")

# Now we construct the presences and absences
eta_sim <- pred_eta(X_sim,alpha=alpha_sim,beta=beta_sim)

# Define the distribution
p_true <- pnorm(eta_sim)
Y_sim <- rbernoulli(1,p_true)

Y_sim[Y_sim==TRUE] <- 1

# You now have simulated presence/absences based on the above species associations and covariate effects,
