#' Detection of global group effect
#'
#' Global test for a set of molecular features (e.g. genes, proteins,...) between two experimental groups. Paired or unpaired design is allowed.
#' @title Detection of global group effect
#' @param X1 Matrix of expression levels in first group. Rows represent features, columns represent samples.
#' @param X2 Matrix of expression levels in second group. Rows represent features, columns represent samples.
#' @param paired FALSE if samples are unpaired, TRUE if samples are paired.
#' @return An object that contains the test results. Contents can be displayed by the summary function.
#' @export
#' @author Klaus Jung
#' @references
#' Brunner, E (2009) Repeated measures under non-sphericity. \emph{Proceedings of the 6th St. Petersburg Workshop on Simulation}, 605-609.
#'
#' Jung K, Becker B, Brunner B and Beissbarth T (2011) Comparison of Global Tests for Functional Gene Sets in Two-Group Designs and Selection of Potentially Effect-causing Genes. \emph{Bioinformatics}, \strong{27}, 1377-1383. \doi{10.1093/bioinformatics/btr152}
#' @examples
#' ### Global comparison of a set of 100 genes between two experimental groups.
#' X1 = matrix(rnorm(1000, 0, 1), 10, 100)
#' X2 = matrix(rnorm(1000, 0.1, 1), 10, 100)
#' RHD = RHighDim(X1, X2, paired=FALSE)
#' summary_RHD(RHD)
RHighDim <- function(X1, X2, paired=TRUE) {
  d = dim(X1)[1]
  if (paired==TRUE) {
    n1 = dim(X1)[2]
    n2 = dim(X2)[2]
    Y = X1 - X2
    H = diag(1, d) - matrix(1, d, d) / d
    Hyp = TestStatSimple(Y, H)
    out = Hyp
  }
  if (paired==FALSE) {
    Y = cbind(X1, X2)
    d = dim(X1)[1]
    n1 = dim(X1)[2]
    n2 = dim(X2)[2]
    N = n1 + n2
    Hyp = TestStatSP(X1, X2)
    out = Hyp
  }
  out
}
