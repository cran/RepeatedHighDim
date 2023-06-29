#' Calculate lower and upper the bounds for pairwise correlations
#'
#' The function calculates upper and lower bounds for pairwise
#' correlations given a vector of marginal probabilities as detailed
#' in Emrich and Piedmonte (1991).
#' @title Calculate lower and upper the bounds for pairwise
#'   correlations
#' @param p Vector of marginal frequencies
#' @param R Correlation matrix
#' @return A list with three entries:
#'  \describe{
#' \item{\emph{L}}{Matrix of lower bounds}
#' \item{\emph{U}}{Matrix of upper bounds}
#' \item{\emph{Z}}{Matrix that indicates whether specified
#' correlations in R are bigger or smaller than the calculated
#' bounds}
#' }
#' @references
#' Emrich, L.J., Piedmonte, M.R.: A method for generating highdimensional multivariate binary variates. \emph{The American Statistician}, \strong{45(4)}, 302 (1991). \doi{10.1080/00031305.1991.10475828}
#' @author Jochen Kruppa, Klaus Jung
#' @export
#' @examples
#' ### A simple example
#' R <- diag(4)
#' p <- c(0.1, 0.2, 0.4, 0.5)
#'
#' rho_bounds(R, p)
rho_bounds <- function(R, p) {
  n = dim(R)[1]
  m = dim(R)[2]
  q = 1 - p
  U = matrix(NA, n, n)
  L = matrix(NA, n, n)
  Z = matrix("OK", n, n)
  for (i in 1:n) {
    for (j in 1:n) {
      L[i,j] = max(c(-sqrt(p[i]*p[j]/q[i]*q[j]), -sqrt(q[i]*q[j]/p[i]*p[j])))
      U[i,j] = min(c(sqrt(p[i]*q[j]/p[j]*q[i]), sqrt(p[j]*q[i]/p[i]*q[j])))
      if (R[i,j]>U[i,j]) Z[i,j] = "big"
      if (R[i,j]<L[i,j]) Z[i,j] = "small"
    }
  }
  diag(U) = NA
  diag(L) = NA
  diag(Z) = NA
  out = list(L=L, U=U, Z=Z)
}
