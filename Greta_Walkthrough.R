#    ---------------------------------------    #
#### 1. Package loading and data preparation ####
#    ---------------------------------------    #

install.packages("JSAM")


load("./Data/MDB_Data.rda")
source("./pred_g_ftn.R")
library(greta)



## Transform objects into greta arrays ##

# Take mini sample of full training set for demo
demo_sample <- unlist(tapply(1:nrow(MDB_Data$Random1_train), MDB_Data$Random1_train$valley, sample, 20))

X_shared <- as_data(as.matrix(MDB_Data$X_train)[demo_sample,])
Y <- as_data(MDB_Data$Y_train[demo_sample,1:12])
Spatial <- as_data(as.numeric(MDB_Data$Random1_train[demo_sample,1]))

demo_species_names <- colnames(MDB_Data$Y_train[demo_sample,1:12])
demo_env_names <- colnames(MDB_Data$X_train)

X_assoc <- X_shared[,2]
X_assoc_pred <- c(0,1,2)

# Define dimensions
n_sites<- nrow(X_shared)
n_env_shared <- ncol(X_shared)
n_latent <- 3
n_valleys <- 23
n_species <- ncol(Y)

# These are our species specific environmental variables
X_ssv <- array(runif(n_sites*n_sites,0,1),dim=c(n_sites,n_species))

X_full <- cbind(X_shared,X_ssv)
n_env <- ncol(X_full)

#    ---------------------    #
#### 2. Model construction ####
#    ---------------------    #

alpha <- normal(0, 1, dim = n_species)

# Theory for traits
# traits is a n_species x (n_traits+1) matrix (needs an intercept))
# g <- normal(0, 10, dim = c(n_env, n_traits))
# beta <- g %*% t(traits)
# Theory for phylogenies
# Same as above but decompose the phylogenic correlation matrix into coordinates (like you did for dsiitancae sbetween lakes)

beta_shared <- normal(0, 10, dim = c(n_env_shared, n_species))

# Create a 40 x 40 matrix with all non-diagonals set to zero
beta_ssv <- zeros(n_env-n_env_shared,n_species)
beta_ssv_diag <- normal(0,10, dim=n_species)
diag(beta_ssv) <- beta_ssv_diag

beta <- rbind(beta_ssv,beta_shared)

# Lambda * z will produce our error matrix
z <- normal(0, 1, dim = c(n_latent, n_sites))

# Define our priors for lambda. Lambda_int will create the correlation matrix later on.
# create_lambda
lambda_int <- create_lambda(n_species,n_latent)
lambda_coef <- create_lambda(n_species,n_latent)

# Now produce our prediction function for our later correlation matrix
R_pred <- pred_R(X_assoc_pred,lambda_int=lambda_int)

# We now add in the hierarchical component
valley_sd <- normal(0, 3, truncation = c(0,Inf))
valley_offset <- normal(0, valley_sd, dim = n_valleys)
valley_effect <- valley_offset[Spatial]

# Now just put the model together
eta <- pred_eta(X_full,Spatial,X_assoc)

# Define the distribution
p <- ilogit(eta)
distribution(Y) <- bernoulli(p)

# Define our correlation matrix
R <- cov2cor(lambda_int %*% t(lambda_int))
lower_idx_R <- lower.tri(R)
R_lower <- R[lower_idx_R]

# Define and plot  model
model_v1 <- model(p,R_lower,beta)
plot(model_v1)

# Aaaaand perform MCMC
demo_draws <- mcmc(model_v1,n_samples = 1000,warmup = 1000)

# Take extra draws if we want to
demo_draws_extra <- extra_samples(demo_draws,n_samples = 2000)


#    -----------    #
#### 3. Analysis ####
#    -----------    #


### Check convergence

# Plot mcmc
plot(calculate(beta_shared,demo_draws_extra))
plot(calculate(beta_ssv,demo_draws_extra))
plot(calculate(R_lower,demo_draws_extra))
plot(calculate(p[1:12,],demo_draws_extra))

# Rhat and effective size checks
coda::gelman.diag(calculate(R_lower,demo_draws_extra),multivariate = FALSE)
summary(coda::effectiveSize(calculate(R_lower,demo_draws_extra)))

### Analyse output


get_beta_list(demo_draws_extra,beta_shared,beta_ssv,species_names=demo_species_names,env_names=demo_env_names,ssv=TRUE)

demo_cooccurrence <- get_cooc_ints(demo_draws_extra,demo_species_names)


# Look at changes in species associations
# Vector as defined in first section gives an increase of 0, 1 and 2 standard deviations in temperature


par(mfrow=c(1,1))
demo_correlation_temp <- pred_correlation(demo_draws_extra,X_assoc_pred,demo_species_names)
corrplot(demo_correlation_temp[[1]], method = "color",type="lower",
         col = colorRampPalette(c("blue", "white", "red"))(200),order="original")
corrplot(demo_correlation_temp[[2]], method = "color",type="lower",
         col = colorRampPalette(c("blue", "white", "red"))(200),order="original")
corrplot(demo_correlation_temp[[3]], method = "color",type="lower",
         col = colorRampPalette(c("blue", "white", "red"))(200),order="original")

### Calculating deviance

# Now we look at deviance
Y_dev <- MDB_Data$Y_train[demo_sample,1:12]

# intercept-only NULL deviance
p_null <- colMeans(Y_dev)
p_null <- kronecker(rep(1, n_sites), t(p_null))
deviance(as.matrix(Y_dev),p_null)

# environment-only deviance
dev_env <- 0
for (i in seq_len(ncol(Y))) {
  dat <- data.frame(y = Y_dev[, i],
                    MDB_Data$X_train[demo_sample, ])
  m <- glm(y ~ ., data = dat,
           family = stats::binomial("probit"))
  pred <- predict(m, type = "response")
  dev_env <- dev_env + deviance(dat$y, pred)
}
dev_env

# full model deviance - should be at least better!
p_dev <- array(colMeans(as.matrix(calculate(p, demo_draws_extra))),dim=dim(p))

deviance(as.matrix(Y_dev),p_dev)

#    --------------    #
#### 4. Predictions ####
#    --------------    #

# Making predictions based on beta values alone

X_shared_new <- as_data(as.matrix(MDB_Data$X_val))
X_ssv_new <- array(runif(nrow(X_shared_new)*12,0,1),dim=c(nrow(X_shared_new),12))
X_full_new <- cbind(X_shared_new,X_ssv_new)
site_no_new <- as_data(as.numeric(MDB_Data$Random1_val[,1]))

predict_species_perc(demo_draws_extra,X_full_new, site_no_new, species_list = demo_species_names)

# Predict conditional probability of one species (given others) at one site

conditional_p_value <- pred_cond(demo_betas_shared,demo_betas_ssv,site_env=site_env,occupancy=occupancy,cooccurrence=demo_cooccurrence,ssv=TRUE,focal_species=12)

# Predict conditional probability of one species (given others) at range of sites

conditional_entire_value <- pred_cond_entire(demo_draws,demo_mean_betas_shared,demo_mean_betas_ssv,site_matrix=as.matrix(MDB_Data$X_train)[demo_sample,],occupancy_matrix=as.matrix(MDB_Data$Y_train)[demo_sample,1:12],temp_vector = as.matrix(MDB_Data$X_train)[demo_sample,2],focal_species=4)
