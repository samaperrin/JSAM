#' Calculates deviance of a given model.
#'
#' @title Calculate deviance
#' @param y Acutal presences/absence of species
#' @param p Predicted likelihood of presence or absence of species.
#'
#' @return Deviance value for model.
#'
#' @export

calc_deviance <- function (y, p, log=TRUE) {
  -2 * sum(dbinom(y, 1, p, log = log))
}
