#' Summary of RHighDim function
#'
#' Summarizes the test results obtained by the RHighDim function.
#'
#' @title Summary of RHighDim function
#' @param object An object provided by the RHighDim function.
#' @param ... additional arguments affecting the summary produced.
#' @return No value
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
#' RHD = RHighDim (X1, X2, paired=FALSE)
#' summary_RHD(RHD)
summary_RHD <-  function(object, ...) {
  A = data.frame(effect=c("Group"), F=round(object$Fn, 4), df1=round(object$f, 4), df2=round(object$f2, 4), p=round(object$p, 4))
  cat("Number of Genes:", object$d, "\n")
  cat("Number of Samples in Group 1:", object$n1, "\n")
  cat("Number of Samples in Group 2:", object$n2, "\n")
  if (object$k==1) cat("Samples are Paired: TRUE", "\n")
  if (object$k==2) cat("Samples are Paired: FALSE", "\n")
  cat("\n")
  print(A)
}
