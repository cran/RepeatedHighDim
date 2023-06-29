#' Detection of global group effect
#'
#' Tests a global effect for a set of molecular features (e.g. genes,
#' proteins,...) between the two groups of samples. Missing values
#' are allowd in the expression data. Samples of the two groups are
#' supposed to be unpaired.
#' @title Detection of global group effect
#' @param X1 Matrix of expression levels in first group. Rows
#'   represent features, columns represent samples.
#' @param X2 Matrix of expression levels in second group. Rows
#'   represent features, columns represent samples.
#' @param nperm Number of permutations.
#' @return The p-value of a permutation test.
#' @export
#' @author Klaus Jung
#' @references
#' Jung K, Dihazi H, Bibi A, Dihazi GH and Beissbarth T (2014): Adaption of the Global Test Idea to Proteomics Data with Missing Values. \emph{Bioinformatics}, \strong{30}, 1424-30. \doi{10.1093/bioinformatics/btu062}
#' @importFrom  nlme lme
#' @importFrom stats anova
#' @examples
#' ### Global comparison of a set of 100 proteins between two experimental groups,
#' ### where (tau * 100) percent of expression levels are missing.
#' n1 = 10
#' n2 = 10
#' d = 100
#' tau = 0.1
#' X1 = t(matrix(rnorm(n1*d, 0, 1), n1, d))
#' X2 = t(matrix(rnorm(n2*d, 0.1, 1), n2, d))
#' X1[sample(1:(n1*d), tau * (n1*d))] = NA
#' X2[sample(1:(n2*d), tau * (n2*d))] = NA
#' GlobTestMissing(X1, X2, nperm=100)
GlobTestMissing  <-  function(X1, X2, nperm=100) {
  d = dim(X1)[1]
  n1 = dim(X1)[2]
  n2 = dim(X2)[2]
  n = n1 + n2
  x = c(as.vector(X1), as.vector(X2))
  group = c(rep(1, d*n1), rep(2, d*n2))
  feature = c(rep(1:d, n1), rep(1:d, n2))
  individual = gl(n, d)

  index = which(!is.na(x))
  x2 = x[index]
  group2 = group[index]
  feature2 = feature[index]
  individual2 = individual[index]

  K = summary(lme(x2 ~ group2 * feature2, random = ~ 1 | individual2))
  A = anova(K)
  P1 = A[[4]][4]
  Pperm = rep(0, nperm)
  Z = cbind(X1, X2)
  for (t in 1:nperm) {
    s = sample(1:n, n, replace=FALSE)
    Z = Z[,s]
    Y1 = Z[,1:n1]
    Y2 = Z[,(n1+1):n]
    x = c(as.vector(Y1), as.vector(Y2))
    individual = s %x% rep(1, d)

    index = which(!is.na(x))
    x2 = x[index]
    group2 = group[index]
    feature2 = feature[index]
    individual2 = individual[index]

    K = summary(lme(x2 ~ group2 * feature2, random = ~ 1 | individual2))
    A = anova(K)
    Pperm[t] = A[[4]][4]
  }
  P = sum(P1>Pperm) / nperm
  return(list(pval=P))
}
