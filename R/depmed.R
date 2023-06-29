#' Calculates the depth median.
#'
#' Calculates the depth median in a specified grid array with given
#' halfspace location depth at each grid location.
#' @title Calculates the depth median.
#' @param G List containing the grid information produced by
#'     \code{\link{gridfun}} and the halfspace location depths
#'     produced by \code{\link{hldepth}}.
#' @return An vector with a length equal to the number of dimension
#'     of the array in G, containing the coordinates of the depth
#'     median.
#' @author Jochen Kruppa, Klaus Jung
#' @references Rousseeuw, P. J., Ruts, I., & Tukey, J. W. (1999). The
#'     bagplot: a bivariate boxplot. The American Statistician,
#'     53(4), 382-387.
#' @importFrom stats median
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
#'# dm <- depmed(G) ## Calculation of the depth median
depmed = function (G)
{
  if (length(G$grid.k)==0) {
    grid.size <- length(G$grid.x)
    maxi <- max(G$H)[1]
    H2 <- (G$H == maxi)
    i2 <- rep(0, grid.size)
    j2 <- rep(0, grid.size)
    k2 <- rep(0, grid.size)
    for (i in 1:grid.size) {
      for (j in 1:grid.size) {
        for (k in 1:grid.size) {
          if (H2[i, j, k] == TRUE) {
            i2[i] <- i
            j2[j] <- j
            k2[k] <- k
          }
        }
      }
    }
    med.i <- median(G$grid.x[which(i2 != 0)])
    med.j <- median(G$grid.y[which(j2 != 0)])
    med.k <- median(G$grid.z[which(k2 != 0)])
    return(c(med.i, med.j, med.k))
  }
  if (length(G$grid.k)>0) {
    k = dim(G$grid.k)[1]
    maxi = max(G$H)[1]
    index = which(G$H==maxi, arr.ind=TRUE)
    med = matrix(NA, k, dim(index)[1])
    for (i in 1:k) med[i,] = G$grid.k[i,index[i]]
    med = apply(med, 1, median)
    return(med)
  }
}
