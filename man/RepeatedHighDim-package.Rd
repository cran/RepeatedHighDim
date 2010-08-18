\name{RepeatedHighDim-package}
\title{Global test for high-dimensional data}
\description{This package provides functions to test for a global group effect in functional gene sets.}
\details{\tabular{ll}{
Package: \tab RepeatedHighDim\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2010-08-17\cr
License: \tab GPL (>= 2)\cr
LazyLoad: \tab yes\cr
}

This package provides functions to test for a global group effect in functional gene sets.}
\alias{RepeatedHighDim-package}
\docType{package}
\author{Klaus Jung \email{Klaus.Jung@ams.med.uni-goettingen.de}}
\keyword{package}
\examples{X1 = matrix(rnorm(1000, 0, 1), 10, 100)
X2 = matrix(rnorm(1000, 0.1, 1), 10, 100)
RHD = RepeatedHighDim(X1, X2, paired=FALSE)
summary(RHD)}
