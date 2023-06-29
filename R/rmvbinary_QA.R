#' Generation of random sample of binary correlated variables
#'
#' The function implements the algorithm proposed by Qaqish (2003) to
#' generate a random sample of d (=length(p)) correlated binary
#' variables. The sample is generated based on given marginal
#' probabilities p of the d variables and their correlation matrix
#' R. The algorithm starts by generating a data for the first
#' variable X_1 and generates succesively the data for X_2, ... based
#' on their conditional probabilities P(X_j|X_[i-1],...,X_1),
#' j=1,...,d.
#' @title Simulating correlated binary variables using the algorithm
#'   by Qaqish (2003)
#' @param n Sample size
#' @param R Correlation matrix
#' @param p Vector of marginal probabilities
#' @return Sample (n x p)-matrix representing a random sample of size n
#'   from the specified multivariate binary distribution.
#' @references
#' Qaqish, B. F. (2003) A family of multivariate binary distributions for simulating correlated binary variables with specified marginal means and correlations. \emph{Biometrika}, \strong{90(2)}, 455-463. \doi{10.1093/biomet/90.2.455}
#' @author Jochen Kruppa, Klaus Jung
#' @importFrom stats rbinom
#' @export
#' @examples
#' ## Generation of a random sample
#' rmvbinary_QA(n = 10, R = diag(2), p = c(0.5, 0.6))
rmvbinary_QA  <-  function(n, R, p) {
  d <- length(p)
  Y  <-  matrix(NA, n, d)
  for (k in 1:n) {
    y  <-  rep(NA, d)
    y[1]  <-  rbinom(1, 1, p[1])
    d  <-  dim(R)[1]
    q  <-  p * (1 - p)
    G  <-  matrix(NA, d, d)
    for (i in 1:d) {
      for (j in 1:d) {
        G[i,j]  <-  R[i,j] * sqrt(q[i] * q[j])
      }}
    for (j in 2:d) {
      Gj  <-  G[1:(j-1),1:(j-1)]
      sj  <-  G[1:(j-1),j]
      bj  <-  solve(Gj) %*% sj
      lambdaj  <-  p[j] + sum(bj * (y[1:(j-1)] - p[1:(j-1)]))
      y[j]  <-  rbinom(1, 1, lambdaj)
    }
    Y[k,]  <-  y
  }
  return(Y)
}
