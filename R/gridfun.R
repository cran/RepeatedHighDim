#' Specifies a k-dimensional array as grid for the calculation of the
#' halfspace location depths.
#'
#' D must have at least three columns. If D has three columns,
#' automatically a 3-dimensional grid is generated. If D has more
#' than three columns, k must be specified.
#' @title Specifies grid for the calculation of the halfspace location depths
#' @param D Data set with rows representing the individuals and
#'     columns representing the features. In the case of three
#'     dimensions, the colnames of D must be c("x", "y", "z").
#' @param grid.size Number of grid points in each dimension.
#' @param k Number of dimensions of the grid. Needs only be specified
#'     if D has more than columns.
#' @return
#' A list containing the following elements:
#' \describe{
#' \item{\emph{H}}{The k-dimensional array.}
#' }
#' In the case of a 3-dimensional array, additional elements are:
#' \describe{
#' \item{\emph{grid.x, grid.y, grid.z}}{The coordinates of the grid points at each dimension.}
#' }
#' In the case that the array has more than three dimensions, additional elements are:
#' \describe{
#' \item{\emph{grid.k}}{A matrix with the coordinates of the grid. Row represents dimensions and columns represent grid points.}
#' }
#' @author Jochen Kruppa, Klaus Jung
#' @export
gridfun <- function (D, grid.size, k=4)
{
  if (dim(D)[2]<3) message("Error: number of dimensions must be greater or equal 3.")
  if (dim(D)[2]==3) {
    grid.x <- seq(min(D$x)[1], max(D$x)[1], length.out = grid.size)
    grid.y <- seq(min(D$y)[1], max(D$y)[1], length.out = grid.size)
    grid.z <- seq(min(D$z)[1], max(D$z)[1], length.out = grid.size)
    H <- array(NA, c(grid.size, grid.size, grid.size))
    return(list(grid.x = grid.x, grid.y = grid.y, grid.z = grid.z, H = H))
  }
  if (dim(D)[2]>3) {
    grid.k <- matrix(NA, k, grid.size)
    for (i in 1:k) grid.k[i,] <- seq(min(D[,i])[1], max(D[,i])[1], length.out = grid.size)
    H <- array(NA, rep(grid.size, k))
    return(list(grid.k = grid.k, H = H))
  }
}
