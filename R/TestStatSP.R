#' Calculates test statistic for RHighDim
#'
#' @importFrom stats pf
#' @param Y1 A matrix.
#' @param Y2 A matrix.
#' @return A list.
TestStatSP <- function(Y1, Y2) {
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
