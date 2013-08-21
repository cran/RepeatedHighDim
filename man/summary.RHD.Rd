\name{summary.RHD}
\alias{summary.RHD}
\title{Summary of RepeatedHighDim function}
\usage{
  \method{summary}{RHD} (object, ...)
}
\arguments{
  \item{object}{An object provided by the RepeatedHighDim
  function.}

  \item{...}{additional arguments affecting the summary
  produced.}
}
\value{
  No value
}
\description{
  Summarizes the test results obtained by the
  RepeatedHighDim function.
}
\examples{
### Global comparison of a set of 100 genes between two experimental groups.
X1 = matrix(rnorm(1000, 0, 1), 10, 100)
X2 = matrix(rnorm(1000, 0.1, 1), 10, 100)
RHD = RepeatedHighDim(X1, X2, paired=FALSE)
summary(RHD)
}
\author{
  Klaus Jung \email{Klaus.Jung@ams.med.uni-goettingen.de}
}
\references{
  \itemize{ \item{Brunner, E (2009) Repeated measures under
  non-sphericity. Proceedings of the 6th St. Petersburg
  Workshop on Simulation, 605-609.} \item{Jung K, Becker B,
  Brunner B and Beissbarth T (2011) Comparison of Global
  Tests for Functional Gene Sets in Two-Group Designs and
  Selection of Potentially Effect-causing Genes.
  Bioinformatics, 27: 1377-1383.} }
}

