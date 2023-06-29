#' Calculates the bag of a gemplot (i.e. the inner gemstone).
#'
#' Determines those grid points that belong to the bag, i.e. a convex
#' hull that contains 50 percent of the data. In the case of a
#' 3-dimensional data set, the bag can be visualized by an inner
#' gemstone that can be accompanied by an outer gemstone (\code{\link{loop}}).
#'
#' @title  Calculates the bag
#' @param D Data set with rows representing the individuals and
#'     columns representing the features. In the case of three
#'     dimensions, the colnames of D must be c("x", "y", "z").
#' @param G List containing the grid information produced by
#'     \code{\link{gridfun}} and the halfspace location depths calculated by
#'     \code{\link{hldepth}}.
#'
#' @return A list containg the following elements:
#' \describe{
#' \item{\emph{coords}}{Coordinates of the grid points that belong to
#'     the bag. Each row represents a grid point and each column
#'     represents one dimension.}
#' \item{\emph{hull}}{A data matrix that
#'     contains the indices of the margin grid points of the bag that
#'     cover the convex hull by triangles. Each row represents one
#'     triangle. The indices correspond to the rows of coords.}
#' }
#'
#' @references
#' Rousseeuw, P. J., Ruts, I., & Tukey, J. W. (1999). The bagplot: a bivariate boxplot. \emph{The American Statistician}, \strong{53(4)}, 382-387. \doi{10.1080/00031305.1999.10474494}

#'
#' Kruppa, J., & Jung, K. (2017). Automated multigroup outlier identification in molecular high-throughput data using bagplots and gemplots. \emph{BMC bioinformatics}, \strong{18(1)}, 1-10. \url{https://link.springer.com/article/10.1186/s12859-017-1645-5}
#' @author Jochen Kruppa, Klaus Jung
#' @importFrom  rgl material3d bg3d points3d text3d spheres3d axes3d
#' @importFrom geometry convhulln
#' @export
#' @examples
#' ## Attention: calculation is currently time-consuming.
#' ## Remove #-Symbols to run examples
#'
#' ## Two 3-dimensional example data sets D1 and D2
#'# n <- 200
#'# x1 <- rnorm(n, 0, 1)
#'# y1 <- rnorm(n, 0, 1)
#'# z1 <- rnorm(n, 0, 1)
#'# D1 <- data.frame(cbind(x1, y1, z1))
#'# x2 <- rnorm(n, 1, 1)
#'# y2 <- rnorm(n, 1, 1)
#'# z2 <- rnorm(n, 1, 1)
#'# D2 <- data.frame(cbind(x2, y2, z2))
#'# colnames(D1) <- c("x", "y", "z")
#'# colnames(D2) <- c("x", "y", "z")
#'
#' ## Placing outliers in D1 and D2
#'# D1[17,] = c(4, 5, 6)
#'# D2[99,] = -c(3, 4, 5)
#'
#' ## Grid size and graphic parameters
#'# grid.size <- 20
#'# red <- rgb(200, 100, 100, alpha = 100, maxColorValue = 255)
#'# blue <- rgb(100, 100, 200, alpha = 100, maxColorValue = 255)
#'# yel <- rgb(255, 255, 102, alpha = 100, maxColorValue = 255)
#'# white <- rgb(255, 255, 255, alpha = 100, maxColorValue = 255)
#'# require(rgl)
#'# material3d(color=c(red, blue, yel, white),
#'# alpha=c(0.5, 0.5, 0.5, 0.5), smooth=FALSE, specular="black")
#'
#' ## Calucation and visualization of gemplot for D1
#'# G <- gridfun(D1, grid.size=20)
#'# G$H <- hldepth(D1, G, verbose=TRUE)
#'# dm <- depmed(G)
#'# B <- bag(D1, G)
#'# L <- loop(D1, B, dm=dm)
#'# bg3d(color = "gray39" )
#'# points3d(D1[L$outliers==0,1], D1[L$outliers==0,2], D1[L$outliers==0,3], col="green")
#'# text3d(D1[L$outliers==1,1], D1[L$outliers==1,2],D1[L$outliers==1,3],
#'# as.character(which(L$outliers==1)), col=yel)
#'# spheres3d(dm[1], dm[2], dm[3], col=yel, radius=0.1)
#'# material3d(1,alpha=0.4)
#'# gem(B$coords, B$hull, red)
#'# gem(L$coords.loop, L$hull.loop, red)
#'# axes3d(col="white")
#'
#' ## Calucation and visualization of gemplot for D2
#'# G <- gridfun(D2, grid.size=20)
#'# G$H <- hldepth(D2, G, verbose=TRUE)
#'# dm <- depmed(G)
#'# B <- bag(D2, G)
#'# L <- loop(D2, B, dm=dm)
#'# points3d(D2[L$outliers==0,1], D2[L$outliers==0,2], D2[L$outliers==0,3], col="green")
#'# text3d(D2[L$outliers==1,1], D2[L$outliers==1,2],D2[L$outliers==1,3],
#'# as.character(which(L$outliers==1)), col=yel)
#'# spheres3d(dm[1], dm[2], dm[3], col=yel, radius=0.1)
#'# gem(B$coords, B$hull, blue)
#'# gem(L$coords.loop, L$hull.loop, blue)
bag <- function (D, G)
{
  if (dim(D)[2]==3) {
    grid.size = dim(G$H)[1]
    n <- dim(D)[1]
    D.k <- rep(NA, n)
    for (i in 1:n) {
      I <- matrix(NA, 8, 3)
      I[1, ] <- c(G$grid.x[max(which(G$grid.x <= D[i, 1]))],
                  G$grid.y[max(which(G$grid.y <= D[i, 2]))],
                  G$grid.z[max(which(G$grid.z <= D[i, 3]))])
      I[2, ] <- c(G$grid.x[max(which(G$grid.x <= D[i, 1]))],
                  G$grid.y[max(which(G$grid.y <= D[i, 2]))],
                  G$grid.z[min(which(G$grid.z >= D[i, 3]))])
      I[3, ] <- c(G$grid.x[max(which(G$grid.x <= D[i, 1]))],
                  G$grid.y[min(which(G$grid.y >= D[i, 2]))],
                  G$grid.z[max(which(G$grid.z <= D[i, 3]))])
      I[4, ] <- c(G$grid.x[max(which(G$grid.x <= D[i, 1]))],
                  G$grid.y[min(which(G$grid.y >= D[i, 2]))],
                  G$grid.z[min(which(G$grid.z >= D[i, 3]))])
      I[5, ] <- c(G$grid.x[min(which(G$grid.x >= D[i, 1]))],
                  G$grid.y[max(which(G$grid.y <= D[i, 2]))],
                  G$grid.z[max(which(G$grid.z <= D[i, 3]))])
      I[6, ] <- c(G$grid.x[min(which(G$grid.x >= D[i, 1]))],
                  G$grid.y[max(which(G$grid.y <= D[i, 2]))],
                  G$grid.z[min(which(G$grid.z >= D[i, 3]))])
      I[7, ] <- c(G$grid.x[min(which(G$grid.x >= D[i, 1]))],
                  G$grid.y[min(which(G$grid.y >= D[i, 2]))],
                  G$grid.z[max(which(G$grid.z <= D[i, 3]))])
      I[8, ] <- c(G$grid.x[min(which(G$grid.x >= D[i, 1]))],
                  G$grid.y[min(which(G$grid.y >= D[i, 2]))],
                  G$grid.z[min(which(G$grid.z >= D[i, 3]))])
      I <- cbind(I, NA)
      for (t in 1:8) {
        index1 <- match(I[t, 1], G$grid.x)
        index2 <- match(I[t, 2], G$grid.y)
        index3 <- match(I[t, 3], G$grid.z)
        I[, 4] <- G$H[index1, index2, index3]
      }
      D.k[i] <- min(I[, 4])
    }
    H2 <- (G$H >= max(which(cumsum(table(D.k)) <= (n/2))))
    BAG <- matrix(NA, 0, 3)
    for (i in 1:grid.size) {
      for (j in 1:grid.size) {
        for (k in 1:grid.size) {
          if (H2[i, j, k] == TRUE) {
            BAG <- rbind(BAG, c(G$grid.x[i], G$grid.y[j],
                                G$grid.z[k]))
          }
        }
      }
    }
    convH <- convhulln(BAG)
  }
  if (dim(D)[2]>3) {
    n = dim(D)[1]
    d = dim(D)[2]
    k = dim(G$grid.k)[2]
    I = matrix(NA, 2^d, d)
    D.k = rep(NA, n)
    for (i in 1:n) {
      U = cbind(G$grid.k, as.numeric(D[i,]))
      U2 = U
      for (t in 1:d) {
        U2[t,] = rank(U[t,], ties.method="first")
        dimnames(G$H)[[t]] = 1:k
      }
      I = t(matrix(U2[,k+1]-1, nrow=d, ncol=2^d))
      for (t in 1:d) I[,t] = I[,t] + as.numeric(gl(2, 2^(d-t), 2^d)) - 1
      I[which(I>k)] = k
      I = cbind(I, NA)
      for (s in 1:(2^d)) I[s,d+1] = G$H[t(matrix(as.character(I[s,1:d])))]
      D.k[i] = min(I[,d+1])
    }
    H2 <- (G$H >= max(which(cumsum(table(D.k)) <= (n/2))))
    H3 <- H2
    H3[1:(k^d)] <- array(1:(k^d), rep(k, d))
    BAG = matrix(NA, table(H2)[2], d)
    index = which(H2==TRUE)
    for (t in 1:length(index)) {
      index2 = which(H3==index[t], arr.ind=TRUE)
      for (j in 1:d) BAG[t,j] = G$grid.k[j,index2[j]]
    }
    convH <- convhulln(BAG)
  }
  return(list(coords = BAG, hull = convH))
}
