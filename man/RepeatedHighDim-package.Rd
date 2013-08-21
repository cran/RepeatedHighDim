\docType{package}
\name{RepeatedHighDim-package}
\alias{RepeatedHighDim-package}
\title{Global tests for expression data of high-dimensional sets of molecular features.}
\description{
  Global tests for expression data of high-dimensional sets
  of molecular features.
}
\details{
  \tabular{ll}{ Package: \tab RepeatedHighDim\cr Type: \tab
  Package\cr Version: \tab 2.0.0\cr Date: \tab
  2013-08-21\cr License: \tab GPL (>= 2)\cr LazyLoad: \tab
  yes\cr }

  Global tests for expression data of high-dimensional sets
  of molecular features.
}
\examples{
### Global comparison of a set of 100 genes between two experimental groups.
X1 = matrix(rnorm(1000, 0, 1), 10, 100)
X2 = matrix(rnorm(1000, 0.1, 1), 10, 100)
RHD = RepeatedHighDim(X1, X2, paired=FALSE)
summary(RHD)

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
\keyword{package}

