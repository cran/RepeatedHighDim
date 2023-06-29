#' Calculates test statistic for RHighDim
#'
#' @param Y A matrix.
#' @param H A matrix.
#' @return A list.
#' @importFrom  MASS ginv
#' @importFrom stats pf
TestStatSimple <- function(Y, H) {
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
