\name{RepeatedHighDim}
\alias{RepeatedHighDim}
\title{Detection of global group effect}
\usage{
  RepeatedHighDim(X1, X2, paired = TRUE)
}
\arguments{
  \item{X1}{Matrix of expression levels in first group.
  Rows represent features, columns represent samples.}

  \item{X2}{Matrix of expression levels in second group.
  Rows represent features, columns represent samples.}

  \item{paired}{FALSE if samples are unpaired, TRUE if
  samples are paired.}
}
\value{
  An object that contains the test results. Contents can be
  displayed by the summary function.
}
\description{
  Global test for a set of molecular features (e.g. genes,
  proteins,...) between two experimental groups. Paired or
  unpaired design is allowed.
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

