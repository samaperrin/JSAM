# 
library(greta)
X_newdata <- X_full
 X_assoc_newdata=X_assoc
site_id_newdata <- as_data(Spatial)

pred_eta <- function(X_newdata, site_id_newdata = NULL, X_assoc_newdata = NULL) {
  
  eta <- pred_env(X_newdata)
  
  if (!is.null(site_id_newdata)) {
    eta <- eta + pred_site(site_id_newdata)
  }
  
  if (!is.null(X_assoc_newdata)) {
    eta <- eta + pred_assoc(X_assoc_newdata)
  }
  
  return(eta)
}

# pred()


pred_env <- function(X_newdata) {
  XB <- X_newdata %*% beta
  eta <- sweep(XB, 2, alpha, "+")
  return(eta)
}


pred_site <- function(site_id_newdata) {
  valley_effect <- valley_offset[site_id_newdata]
  
  valleyB <- sweep(ones(nrow(valley_effect),n_species), 1, valley_effect, "+")
  return(valleyB)
  }

pred_assoc <- function(X_assoc) {
  lambda_effect_rep <- kronecker(X_assoc, lambda_coef) # (site-lambdas stacked on top of one another)
  lambda_int_rep <- kronecker(ones(n_sites, 1), lambda_int)
  lambda_rep <- lambda_int_rep + lambda_effect_rep
  
  #Do the same thing for our zs so we can perform multiplication
  z_rep <- kronecker(ones(n_species, 1), t(z))
  
  t_e <- rowSums(lambda_rep * z_rep)
  dim(t_e) <- c(n_species, n_sites)
  e <- t(t_e)
  return(e)
}

# This SHOULD give us our new R matrix, but it needs work.
#X_assoc_pred <- X_assoc_pred[1,]
# lambda_int <- create_lambda(n_species,n_latent)

pred_one_R <- function(X_assoc_pred = NULL, lower_only = FALSE, lambda_int){
  
  lambda <- lambda_int
  if (!is.null(X_assoc_pred)) {
    lambda <- lambda + lambda_coef * X_assoc_pred
  }
  
  R <- greta::cov2cor(lambda %*% t(lambda))
  
  if (lower_only) {
    lower_idx_R <- lower.tri(R)
    R <- R[lower_idx_R]
  }
  
  return(R)
}


# This SHOULD give us our new R matrix, but it needs work. new versioN

pred_R <- function(X_assoc_pred, lower_only = FALSE, lambda_int) {

  if (!inherits(X_assoc_pred, "greta_array")) {
    X_assoc_pred <- as.matrix(X_assoc_pred)
  }
  
  n <- nrow(X_assoc_pred)
  
  result <- list()
  
  for (i in seq_len(n)) {
    result[[i]] <- pred_one_R(X_assoc_pred[i, ], lower_only = lower_only, lambda_int)
  }
  
  result
  
}
