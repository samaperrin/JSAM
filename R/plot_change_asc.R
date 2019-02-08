#' Takes association change matrix and plots it, visualizing the changes in species association between one focal species and others.
#'
#' @title Plot species association changes.
#' @param association_gradient_change_matrix A species by environmental covariate level matrix.
#'
#' @return A plot showing changes in species associations.
#'
#' @export
#'
plot_change_asc <- function(association_gradient_change_matrix,legend_position=c(2,1)) {
  n_temps <- ncol(association_gradient_change_matrix)
  n_species <- nrow(association_gradient_change_matrix)
  plot(seq(-1,1,length.out=n_temps) ~ c(1:n_temps),type='n',xlab="Temperature",ylab="Association")
  cl <- brewer.pal(n_species, "Set1")
  for(i in 1:n_temps) {
    lines(seq(1,n_temps,1) ,association_gradient_change_matrix[i,],col=cl[i])
  }
  legend(legend_position[1], legend_position[2], legend=rownames(association_gradient_change_matrix), fill = cl)
}
