\name{RepeatedHighDim}
\alias{RepeatedHighDim}
\title{Detection of global group effect...}
\usage{RepeatedHighDim(X1, X2, paired=TRUE)}
\description{Detection of global group effect}
\details{Detects global group effect between the samples of two groups (paired or unpaired).}
\value{An object that contains the test results. Contents can be displayes by the summary function.}
\author{Klaus Jung \email{Klaus.Jung@ams.med.uni-goettingen.de}}
\references{Brunner, E. (2009) Repeated measures under non-sphericity. Proceedings of the 6th St. Petersburg Workshop on Simulation, 605-609.}
\arguments{\item{X1}{Matrix of expression levels in first group. Rows represent genes, columns represent samples.}
\item{X2}{Matrix of expression levels in second group. Rows represent genes, columns represent samples.}
\item{paired}{FALSE if samples are unpaired, TRUE if samples are paired.}}
\examples{X1 = matrix(rnorm(1000, 0, 1), 10, 100)
X2 = matrix(rnorm(1000, 0.1, 1), 10, 100)
RHD = RepeatedHighDim(X1, X2, paired=FALSE)
summary(RHD)}
