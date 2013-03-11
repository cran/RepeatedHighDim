\name{TestStatSimple}
\alias{TestStatSimple}
\title{Calculation of test statistic}
\usage{
  TestStatSimple(Y, H)
}
\arguments{
  \item{Y}{Matrix with differences of paires. Rows
  represent genes, columns represent samples.}

  \item{H}{Hypothesis matrix.}
}
\value{
  A list containing the following items: \item{k}{Indicates
  whether the paired or unpaired case was tested.}
  \item{d}{Number of genes.} \item{n1}{Number of samples in
  group 1.} \item{n2}{Number of samples in group 2.}
  \item{Fn}{Test statistic.} \item{f}{First degree of
  freedoms.} \item{f2}{Second degree of freedom.}
  \item{p}{p-value.}
}
\description{
  Calculates the test statistic in the case of paired
  samples.
}
\examples{
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

