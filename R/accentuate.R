#' Takes co-occurrence matrix and multiplies paffects so that they are easier to view, with the absolute maximum covariate being increased to 1.
#'
#' @title Magnification of co-occurrence matrix.
#' @param correlation_matrix_list An N x N co-occurrence matrix.
#'
#' @return An N x N co-occurrence matrix with increased covariate effects.
#'
#' @export
#'
accentuate <- function(correlation_matrix_list) {
  unlisted <- unlist(correlation_matrix_list)

  times_factor <- 1/max(abs(unlisted))
  accentuated <- lapply(correlation_matrix_list, "*",times_factor)
  return(list(accentuated=accentuated,factor=times_factor))
}
