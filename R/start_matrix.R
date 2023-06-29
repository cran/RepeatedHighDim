#' Generation of the start matrix with n rows and specified marginal
#' probabilities p.
#'
#' The start matrix needs to be setup for further use in the genetic
#' algorithm implemented in the function \code{\link{iter_matrix}}. For
#' high-dimensional cases or if the marginal probabilities have
#' multiple decimal places, the number k of rows should be large (up
#' to multiple thousand).
#' @title Setup of the start matrix
#' @param p Marginal probabilities of the start matrix.
#' @param k Number of rows to be generated.
#' @return A (k x p)-Matrix with with entries 0 and 1 according to
#'   the specified marginal probabilities p.
#' @author Jochen Kruppa, Klaus Jung
#' @references
#' Kruppa, J., Lepenies, B., & Jung, K. (2018). A genetic algorithm for simulating correlated binary data from biomedical research. \emph{Computers in biology and medicine}, \strong{92}, 1-8. \doi{10.1016/j.compbiomed.2017.10.023}
#' @export
#' @examples
#' X0 <- start_matrix(p = c(0.5, 0.6), k = 10000)
#'
#' ## check if p can be restored
#' apply(X0, 2, mean)
start_matrix <- function(p, k) {
  m = length(p)
  X0 = matrix(NA, k, m)
  for (j in 1:m) {
    p0 = k*(1-p[j])
    p1 = k*p[j]
    vec = c(rep(0, p0), rep(1, p1))
    if (length(vec)<k) vec = c(vec, 1)
    X0[,j] = vec
  }
  return(X0)
}
