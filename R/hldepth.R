#' Calculates the halfspace location depth for each point in a given grid.
#'
#' Calculation of the halfspace location depth at each grid point is
#' mandatory before calculating the depth median
#' (\code{\link{depmed}}), the bag (\code{\link{bag}}) and the loop
#' (\code{\link{loop}}). Ideally, the output is assigned to the array
#' H produced by \code{\link{gridfun}}.
#' @title Calculates the halfspace location depth
#' @param D Data set with rows representing the individuals and
#'     columns representing the features. In the case of three
#'     dimensions, the colnames of D must be c("x", "y", "z").
#' @param G List containing the grid information produced by
#'     \code{\link{gridfun}}.
#' @param verbose Logical. Indicates whether progress information is
#'     printed during calculation.
#' @return
#' \describe{
#' \item{\emph{H}}{An array of the same dimension as the array in argument G. The elements contain the halfspace location depth at the related grid location.}
#' }
#'
#' @references Rousseeuw, P. J., Ruts, I., & Tukey, J. W. (1999). The
#'     bagplot: a bivariate boxplot. The American Statistician,
#'     53(4), 382-387.
#' @author Jochen Kruppa, Klaus Jung
#' @importFrom ddalpha depth.halfspace
#' @export
#' @examples
#' ## Attention: calculation is currently time-consuming.
#' ## Remove #-Symbols to run examples
#'
#' ## A 3-dimensional example data set D1
#'# n <- 200
#'# x1 <- rnorm(n, 0, 1)
#'# y1 <- rnorm(n, 0, 1)
#'# z1 <- rnorm(n, 0, 1)
#'# D1 <- data.frame(cbind(x1, y1, z1))
#'# colnames(D1) <- c("x", "y", "z")
#'
#' ## Specification of the grid and calculation of the halfspace location depth at each grid location.
#'# G <- gridfun(D1, grid.size=20)
#'# G$H <- hldepth(D1, G, verbose=TRUE)
hldepth = function (D, G, verbose = TRUE)
{
  if (dim(D)[2]==3) {
    n <- dim(D)[1]
    grid.size <- length(G$grid.x)
    perc <- 10
    H <- G$H
    for (i in 1:grid.size) {
      for (j in 1:grid.size) {
        for (k in 1:grid.size) {
          u <- c(G$grid.x[i], G$grid.y[j], G$grid.z[k])
          H[i, j, k] <- n * depth.halfspace(x=u,data= D,exact=T)
        }
      }
      if (100 * i/grid.size >= perc && verbose == TRUE) {
        message(paste("Calculation of halfspace location depths: ",
                      round(100 * i/grid.size, 0), " % of grid points done",
                      sep = ""))
        perc <- perc + 10
      }
    }
  }
  if (dim(D)[2]>3) {
    n <- dim(D)[1]
    k = dim(G$grid.k)[1]
    grid.size <- dim(G$grid.k)[2]
    grid.kk <- grid.size^k
    H <- G$H
    H2 <- H
    H2[1:grid.kk] <- array(1:grid.kk, rep(grid.size, k))
    perc <- 10
    for (j in 1:grid.kk) {
      index <- which(H2==j, arr.ind=TRUE)
      u <- rep(NA, k)
      for (i in 1:k) u[i] <- G$grid.k[i,index[i]]
      H[j] <- n * depth.halfspace(x=u, data=D,exact=F)
      if (100 * j/grid.kk >= perc && verbose == TRUE) {
        message(paste("Calculation of halfspace location depths: ",
                      round(100 * j/grid.kk, 0), " % of grid points done",
                      sep = ""))
        perc <- perc + 10
      }
    }
  }
  return(H)
}
