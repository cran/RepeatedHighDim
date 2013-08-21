\name{GlobTestMissing}
\alias{GlobTestMissing}
\title{Detection of global group effect}
\usage{
  GlobTestMissing(X1, X2, nperm = 100)
}
\arguments{
  \item{X1}{Matrix of expression levels in first group.
  Rows represent features, columns represent samples.}

  \item{X2}{Matrix of expression levels in second group.
  Rows represent features, columns represent samples.}

  \item{nperm}{Number of permutations.}
}
\value{
  The p-value of a permutation test.
}
\description{
  Tests a global effect for a set of molecular features
  (e.g. genes, proteins,...) between the two groups of
  samples. Missing values are allowd in the expression
  data. Samples of the two groups are supposed to be
  unpaired.
}
\examples{
### Global comparison of a set of 100 proteins between two experimental groups,
### where (tau * 100) percent of expression levels are missing.
n1 = 10
n2 = 10
d = 100
tau = 0.1
X1 = t(matrix(rnorm(n1*d, 0, 1), n1, d))
X2 = t(matrix(rnorm(n2*d, 0.1, 1), n2, d))
X1[sample(1:(n1*d), tau * (n1*d))] = NA
X2[sample(1:(n2*d), tau * (n2*d))] = NA
GlobTestMissing(X1, X2, nperm=100)
}
\author{
  Klaus Jung \email{Klaus.Jung@ams.med.uni-goettingen.de}
}

