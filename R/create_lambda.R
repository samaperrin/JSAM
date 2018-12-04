#' Creates restrained lambda matrix for us in JSAM model.
#'
#' @title Creates lambda matrix.
#' @param n_species Number of species involved in model.
#' @param n_latent Number of latent variables chosen.
#'
#' @return A species x latent variable GRETA matrix with a uniformly distributed diagonal and an upper right triangle of zeroes.
#'
#' @export

create_lambda <- function(n_species,n_latent) {
  lambda <- zeros(n_species, n_latent)
  diag(lambda) <- normal(0, 1, dim = n_latent, truncation = c(0,Inf))
  lower_idx <- lower.tri(lambda)
  lambda[lower_idx] <- normal(0, 1, dim = sum(lower_idx))
  return(lambda)
}
