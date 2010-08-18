#' This package provides functions to test for a global group effect in functional gene sets.
#'
#' \tabular{ll}{
#' Package: \tab RepeatedHighDim\cr
#' Type: \tab Package\cr
#' Version: \tab 1.0\cr
#' Date: \tab 2010-08-17\cr
#' License: \tab GPL (>= 2)\cr
#' LazyLoad: \tab yes\cr
#' }
#'
#' This package provides functions to test for a global group effect in functional gene sets.
#'
#' @name RepeatedHighDim-package
#' @aliases RepeatedHighDim-package
#' @docType package
#' @title Global test for high-dimensional data
#' @author Klaus Jung \email{Klaus.Jung@@ams.med.uni-goettingen.de}
#' @keywords package
#' @examples
#' X1 = matrix(rnorm(1000, 0, 1), 10, 100)
#' X2 = matrix(rnorm(1000, 0.1, 1), 10, 100)
#' RHD = RepeatedHighDim(X1, X2, paired=FALSE)
#' summary(RHD)
NA



#' Calculation of test statistic
#'
#' Calculates the test statistic in the case of paired samples.
#'
#' @param Y Matrix with differences of paires. Rows represent genes, columns represent samples.
#' @param H Hypothesis matrix.
#' @return A list containing the following items:
#'   \item{k}{Indicates whether the paired or unpaired case was tested.}
#'   \item{d}{Number of genes.}
#'   \item{n1}{Number of samples in group 1.}
#'   \item{n2}{Number of samples in group 2.}
#'   \item{Fn}{Test statistic.}
#'   \item{f}{First degree of freedoms.}
#'   \item{f2}{Second degree of freedom.}
#'   \item{p}{p-value.}
#' @export
#' @author Klaus Jung \email{Klaus.Jung@@ams.med.uni-goettingen.de}
#' @references Brunner, E. (2009) Repeated measures under non-sphericity. Proceedings of the 6th St. Petersburg Workshop on Simulation, 605-609.
#' @examples
#' X1 = matrix(rnorm(1000, 0, 1), 10, 100)
#' X2 = matrix(rnorm(1000, 0.1, 1), 10, 100)
#' RHD = RepeatedHighDim(X1, X2, paired=FALSE)
#' summary(RHD)
TestStatSimple = function(Y, H) {
	d = dim(Y)[1]
	n = dim(Y)[2]
	Pn = diag(1, n) - 1/n
	Z = Pn %*% t(Y)
	TraceSigma = sum(diag(Z %*% t(Z))) / (n - 1)
	T = t(H) %*% ginv(H %*% t(H)) %*% H
	barY = apply(Y, 1, mean)
	barZ = T %*% barY
	Qn = t(barZ) %*% barZ
	Z = T %*% Y
	Fn = n * Qn / TraceSigma
	A = t(Z) %*% Z
	AA = diag(A) %*% t(diag(A))
	B1 =  2 * (sum(AA[upper.tri(AA)])) / (n * (n - 1))
	AA = A * A
	B2 =  2 * (sum(AA[upper.tri(AA)])) / (n * (n - 1))
	f = B1 / B2
	p = 1 - pf(Fn, f, (n - 1) * f)
	out = list(k=1, d=d, n1=n, n2=n, Fn=Fn, f=f, f2=(f * (n-1)), p=p)
}


#' Calculation of test statistic
#'
#' Calculates the test statistic in the case of unpaired samples.
#'
#' @param Y1 Matrix of expression levels in first group. Rows represent genes, columns represent samples.
#' @param Y2 Matrix of expression levels in second group. Rows represent genes, columns represent samples.
#' @return A list containing the following items:
#'   \item{k}{Indicates whether the paired or unpaired case was tested.}
#'   \item{d}{Number of genes.}
#'   \item{n1}{Number of samples in group 1.}
#'   \item{n2}{Number of samples in group 2.}
#'   \item{Fn}{Test statistic.}
#'   \item{f}{First degree of freedoms.}
#'   \item{f2}{Second degree of freedom.}
#'   \item{p}{p-value.}
#' @export
#' @author Klaus Jung \email{Klaus.Jung@@ams.med.uni-goettingen.de}
#' @references Brunner, E. (2009) Repeated measures under non-sphericity. Proceedings of the 6th St. Petersburg Workshop on Simulation, 605-609.
#' @examples
#' X1 = matrix(rnorm(1000, 0, 1), 10, 100)
#' X2 = matrix(rnorm(1000, 0.1, 1), 10, 100)
#' RHD = RepeatedHighDim(X1, X2, paired=FALSE)
#' summary(RHD)
TestStatSP = function(Y1, Y2) {
	d = dim(Y1)[1]
	n1 = dim(Y1)[2]
	n2 = dim(Y2)[2]
	N = n1 + n2
	barZ1 = apply(Y1, 1, sum) / n1
	barZ2 = apply(Y2, 1, sum) / n2
	QN = t(barZ1 - barZ2) %*% (barZ1 - barZ2)
	Pn1 = diag(1, n1) - matrix(1, n1, n1) / n1
	Pn2 = diag(1, n2) - matrix(1, n2, n2) / n2
	Z1 = Pn1 %*% t(Y1)
	Z2 = Pn2 %*% t(Y2)
	scm1t = (Z1 %*% t(Z1)) / (n1 - 1)
	scm2t = (Z2 %*% t(Z2)) / (n2 - 1)
	TraceSigma1 = sum(diag(scm1t))
	TraceSigma2 = sum(diag(scm2t)) 
	hatSN = TraceSigma1 / n1 + TraceSigma2 / n2
	FN = QN / hatSN
	M1 = Y1
	M2 = Y2
	A1 = t(M1) %*% M1
	A2 = t(M2) %*% M2
	A12 = t(t(M1) %*% M2)
	a1 = diag(A1)
	a2 = diag(A2)
	ONEn1 = rep(1, n1)
	ONEn2 = rep(1, n2)
	K1 = a1 %*% t(ONEn1) + ONEn1 %*% t(a1) - 2 * A1
	K2 = a2 %*% t(ONEn2) + ONEn2 %*% t(a2) - 2 * A2
	tildeB1.1 = (sum(K1)^2 - 4 * sum(K1 %*% t(K1)) + 2 * sum(K1 * K1)) / 4
	tildeB1.2 = (sum(K2)^2 - 4 * sum(K2 %*% t(K2)) + 2 * sum(K2 * K2)) / 4
	B1.1 = tildeB1.1 / (n1 * (n1 - 1) * (n1 - 2) * (n1 - 3))
	B1.2 = tildeB1.2 / (n2 * (n2 - 1) * (n2 - 2) * (n2 - 3))
	C1 = (sum(K1) * sum(K2)) / (4 * n1 * (n1 - 1) * n2 * (n2 - 1))
	A01 = A1
	diag(A01) = 0
	A02 = A2
	diag(A02) = 0
	B21.1 = (n1 - 2) * (n1 - 3) * sum(A01 * A01)
	B22.1 = 2 * (n1 - 3) * sum((A01 %*% A01) * (matrix(1, n1, n1) - diag(n1)))
	B23.1 = sum(A01)^2 - 2 * sum(A01 * A01) - 4 * sum((A01 %*% A01) * (matrix(1, n1, n1) - diag(n1)))
	B2.1 =  (B21.1 - B22.1 + B23.1) / (n1 * (n1 - 1) * (n1 - 2) * (n1 - 3))
	B21.2 = (n2 - 2) * (n2 - 3) * sum(A02 * A02)
	B22.2 = 2 * (n2 - 3) * sum((A02 %*% A02) * (matrix(1, n2, n2) - diag(n2)))
	B23.2 = sum(A02)^2 - 2 * sum(A02 * A02) - 4 * sum((A02 %*% A02) * (matrix(1, n2, n2) - diag(n2)))
	B2.2 =  (B21.2 - B22.2 + B23.2) / (n2 * (n2 - 1) * (n2 - 2) * (n2 - 3))
	C21 = n1 * n2 * sum(A12 * A12)
	C22 = n1 * sum(A12 %*% t(A12))
	C23 = n2 * sum(t(A12) %*% A12)
	C24 = sum(t(M1) %*% M2)^2
	C2 = (C21 - C22 - C23 + C24) / (n1 * (n1 - 1) * n2 * (n2 - 1))
	hatf.numer = (B1.1 / n1^2) + (B1.2 / n2^2) + 2 * (C1 / (n1 * n2))
	hatf.denom = (B2.1 / n1^2) + (B2.2 / n2^2) + 2 * (C2 / (n1 * n2))
	hatf = hatf.numer / hatf.denom
	hatf0.numer = hatf.numer
	hatf0.denom = (B2.1 / (n1^2 * (n1 - 1))) + (B2.2 / (n2^2 * (n2 - 1)))
	hatf0 = hatf0.numer / hatf0.denom
	p = 1 - pf(FN, hatf, hatf0)
	out = list(k=2, d=d, n1=n1, n2=n2, Fn=FN, f=hatf, f2=hatf0, p=p)
}



#' Detection of global group effect
#'
#' Detects global group effect between the samples of two groups (paired or unpaired).
#'
#' @param X1 Matrix of expression levels in first group. Rows represent genes, columns represent samples.
#' @param X2 Matrix of expression levels in second group. Rows represent genes, columns represent samples.
#' @param paired FALSE if samples are unpaired, TRUE if samples are paired.
#' @return An object that contains the test results. Contents can be displayes by the summary function.
#' @export
#' @author Klaus Jung \email{Klaus.Jung@@ams.med.uni-goettingen.de}
#' @references Brunner, E. (2009) Repeated measures under non-sphericity. Proceedings of the 6th St. Petersburg Workshop on Simulation, 605-609.
#' @examples
#' X1 = matrix(rnorm(1000, 0, 1), 10, 100)
#' X2 = matrix(rnorm(1000, 0.1, 1), 10, 100)
#' RHD = RepeatedHighDim(X1, X2, paired=FALSE)
#' summary(RHD)
RepeatedHighDim = function(X1, X2, paired=TRUE) {
	d = dim(X1)[1]
	if (paired==TRUE) {
		n1 = dim(X1)[2]
		n2 = dim(X2)[2]
		Y = X1 - X2
		H = diag(1, d) - matrix(1, d, d) / d
		Hyp = TestStatSimple(Y, H)
		out = Hyp
		class(out) = "RHD"
	}
	if (paired==FALSE) {
		Y = cbind(X1, X2)
		d = dim(X1)[1]
		n1 = dim(X1)[2]
		n2 = dim(X2)[2]
		N = n1 + n2
		Hyp = TestStatSP(X1, X2)
		out = Hyp
		class(out) = "RHD"
	}
	out
}


#' Summary of RepeatedHighDim function
#'
#' Summarizes the test results obtained by the RepeatedHighDim function.
#'
#' @param object An object provided by the RepeatedHighDim function.
#' @param ... additional arguments affecting the summary produced.
#' @return No value
#' @export
#' @author Klaus Jung \email{Klaus.Jung@@ams.med.uni-goettingen.de}
#' @references Brunner, E. (2009) Repeated measures under non-sphericity. Proceedings of the 6th St. Petersburg Workshop on Simulation, 605-609.
#' @examples
#' X1 = matrix(rnorm(1000, 0, 1), 10, 100)
#' X2 = matrix(rnorm(1000, 0.1, 1), 10, 100)
#' RHD = RepeatedHighDim(X1, X2, paired=FALSE)
#' summary(RHD)
summary.RHD = function(object, ...) {
	A = data.frame(effect=c("Group"), F=round(object$Fn, 4), df1=round(object$f, 4), df2=round(object$f2, 4), p=round(object$p, 4))	
	cat("Number of Genes:", object$d, "\n")
	cat("Number of Samples in Group 1:", object$n1, "\n")
	cat("Number of Samples in Group 2:", object$n2, "\n")
	if (object$k==1) cat("Samples are Paired: TRUE", "\n")
	if (object$k==2) cat("Samples are Paired: FALSE", "\n")
	cat("\n")
	print(A)
}
