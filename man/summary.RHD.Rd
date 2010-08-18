\name{summary.RHD}
\alias{summary.RHD}
\title{Summary of RepeatedHighDim function...}
\usage{summary.RHD(object, ...)}
\description{Summary of RepeatedHighDim function}
\details{Summarizes the test results obtained by the RepeatedHighDim function.}
\value{No value}
\author{Klaus Jung \email{Klaus.Jung@ams.med.uni-goettingen.de}}
\references{Brunner, E. (2009) Repeated measures under non-sphericity. Proceedings of the 6th St. Petersburg Workshop on Simulation, 605-609.}
\arguments{\item{object}{An object provided by the RepeatedHighDim function.}
\item{...}{additional arguments affecting the summary produced.}}
\examples{X1 = matrix(rnorm(1000, 0, 1), 10, 100)
X2 = matrix(rnorm(1000, 0.1, 1), 10, 100)
RHD = RepeatedHighDim(X1, X2, paired=FALSE)
summary(RHD)}
