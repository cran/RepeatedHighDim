#' Generation of random sample of binary correlated variables
#'
#' The function implements the algorithm proposed by Emrich and
#' Piedmonte (1991) to generate a random sample of d (=length(p))
#' correlated binary variables. The sample is generated based on
#' given marginal probabilities p of the d variables and their
#' correlation matrix R. The algorithm generates first determines an
#' appropriate correlation matrix R' for the multivariate normal
#' distribution. Next, a sample is drawn from N_d(0, R') and each
#' variable is finnaly dichotomized with respect to p.
#' @title Simulating correlated binary variables using the algorithm
#'   by Emrich and Piedmonte (1991)
#' @param n Sample size
#' @param R Correlation matrix
#' @param p Vector of marginal probabilities
#' @return Sample (n x p)-matrix with representing a random sample of size n from the specified multivariate binary distribution.
#' @references
#' Emrich, L.J., Piedmonte, M.R. (1991) A method for generating highdimensional multivariate binary variates. \emph{The American Statistician}, \strong{45(4)}, 302. \doi{10.1080/00031305.1991.10475828}
#' @author Jochen Kruppa, Klaus Jung
#' @importFrom  mvtnorm pmvnorm
#' @importFrom  mvtnorm rmvnorm
#' @importFrom stats qnorm
#' @export
#' @examples
#' ## Generation of a random sample
#' rmvbinary_EP(n = 10, R = diag(2), p = c(0.5, 0.6))
rmvbinary_EP  <-  function(n, R, p){
  s0  <-  seq(-1, 1, 0.01)
  q  <-  1 - p
  K  <-  length(s0)
  d  <-  dim(R)[1]
  S  <-  matrix(NA, d, d)
  S2  <-  matrix(NA, d, d)
  for (i in 1:d) {
    for (j in 1:d) {
      right  <-  R[i,j] * sqrt(p[i] * q[i] * p[j] * q[j]) + (p[i] * p[j])
      left  <-  rep(NA, K)
      for (k in 1:K) {
        S0  <-  diag(2)
        S0[1,2]  <-  s0[k]
        S0[2,1]  <-  s0[k]
        zi  <-  qnorm(p[i], 0, 1)
        zj  <-  qnorm(p[j], 0, 1)
        left[k]  <-  pmvnorm(lower = c(-Inf, -Inf),
                             upper=c(zi, zj),
                             mean=c(0, 0),
                             corr=S0)
      }
      difference  <-  abs(left - right)
      if (R[i,j]<0) S[i,j]  <-  s0[min(which(difference == min(difference)))]
      if (R[i,j]>=0) S[i,j]  <-  s0[max(which(difference == min(difference)))]
    }}
  X0  <-  rmvnorm(n, mean = rep(0, d), sigma = S)
  X  <-  matrix(0, n, d)
  for (j in 1:d) X[which(X0[,j] <= qnorm(p[j], 0, 1)),j]  <-  1
  return(X)
}
