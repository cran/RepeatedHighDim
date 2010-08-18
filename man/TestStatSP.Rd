\name{TestStatSP}
\alias{TestStatSP}
\title{Calculation of test statistic...}
\usage{TestStatSP(Y1, Y2)}
\description{Calculation of test statistic}
\details{Calculates the test statistic in the case of unpaired samples.}
\value{A list containing the following items:
\item{k}{Indicates whether the paired or unpaired case was tested.}
\item{d}{Number of genes.}
\item{n1}{Number of samples in group 1.}
\item{n2}{Number of samples in group 2.}
\item{Fn}{Test statistic.}
\item{f}{First degree of freedoms.}
\item{f2}{Second degree of freedom.}
\item{p}{p-value.}}
\author{Klaus Jung \email{Klaus.Jung@ams.med.uni-goettingen.de}}
\references{Brunner, E. (2009) Repeated measures under non-sphericity. Proceedings of the 6th St. Petersburg Workshop on Simulation, 605-609.}
\arguments{\item{Y1}{Matrix of expression levels in first group. Rows represent genes, columns represent samples.}
\item{Y2}{Matrix of expression levels in second group. Rows represent genes, columns represent samples.}}
\examples{X1 = matrix(rnorm(1000, 0, 1), 10, 100)
X2 = matrix(rnorm(1000, 0.1, 1), 10, 100)
RHD = RepeatedHighDim(X1, X2, paired=FALSE)
summary(RHD)}
