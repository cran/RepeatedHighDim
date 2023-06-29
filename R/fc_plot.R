#' Volcano plot of adjusted confidence intervals
#'
#' Volcano plot of adjusted confidence intervals
#' @title Volcano plot of adjusted confidence intervals
#' @param CI Object as returned from the function fc_ci
#' @param alpha 1 - confidence level (e.g., if confidence level is 0.95, alpha is 0.05)
#' @param updown Character, 'all' if CIs for all genes, 'down' if CIs for down-regulated genes, or 'up' if CIs for up-regulated genes to be plotted
#' @param xlim Vector of length 2 with the lower and upper limits for the X-axis
#' @param ylim Vector of length 2 with the lower and upper limits for the Y-axis. Please note, that p-values are usually displayed on the -log10-scale in a volcano plot
#' @return NULL
#' @export
#' @author Klaus Jung
#' @references
#' Dudoit, S., Shaffer, J. P., & Boldrick, J. C. (2003). Multiple hypothesis testing in microarray experiments. \emph{Statistical Science}, \strong{18(1)}, 71-103. \url{https://projecteuclid.org/journals/statistical-science/volume-18/issue-1/Multiple-Hypothesis-Testing-in-Microarray-Experiments/10.1214/ss/1056397487.full}
#'
#' Jung, K., Friede, T., & Beißbarth, T. (2011). Reporting FDR analogous confidence intervals for the log fold change of differentially expressed genes. \emph{BMC bioinformatics}, \strong{12}, 1-9. \url{https://link.springer.com/article/10.1186/1471-2105-12-288}
#' @importFrom graphics abline axis points
#' @export
#' @examples
#' ### Artificial microarray data
#' d = 1000 ### Number of genes
#' n = 10 ### Sample per group
#' fc = rlnorm(d, 0, 0.1)
#' mu1 = rlnorm(d, 0, 1) ### Mean vector group 1
#' mu2 = mu1 * fc ### Mean vector group 2
#' sd1 = rnorm(d, 1, 0.2)
#' sd2 = rnorm(d, 1, 0.2)
#' X1 = matrix(NA, d, n) ### Expression levels group 1
#' X2 = matrix(NA, d, n) ### Expression levels group 2
#' for (i in 1:n) {
#'   X1[,i] = rnorm(d, mu1, sd=sd1)
#'   X2[,i] = rnorm(d, mu2, sd=sd2)
#' }
#' X = cbind(X1, X2)
#' heatmap(X)
#'
#' ### Differential expression analysis with limma
#' if(check_limma()){
#' group = gl(2, n)
#' design = model.matrix(~ group)
#' fit1 = limma::lmFit(X, design)
#' fit = limma::eBayes(fit1)
#'
#' ### Calculation of confidence intervals
#' CI = fc_ci(fit=fit, alpha=0.05, method="raw")
#' head(CI)
#' CI = fc_ci(fit=fit, alpha=0.05, method="BH")
#' head(CI)
#' CI = fc_ci(fit=fit, alpha=0.05, method="BY")
#' head(CI)
#'
#' fc_plot(CI, xlim=c(-0.5, 3), ylim=-log10(c(1, 0.0001)), updown="up")
#' fc_plot(CI, xlim=c(-3, 0.5), ylim=-log10(c(1, 0.0001)), updown="down")
#' fc_plot(CI, xlim=c(-3, 3), ylim=-log10(c(1, 0.0001)), updown="all")
#' }
fc_plot = function(CI, alpha=0.05, updown="all", xlim=c(-3, 3), ylim=-log10(c(1, 0.001))) {
  logFC = CI$logFC
  d = length(logFC)
  sel = seq(1:d)
  if (updown=="up") sel = sel[-which(logFC<0)]
  if (updown=="down") sel = sel[-which(logFC>0)]

  if (length(grep("raw", colnames(CI)[1]))==1) {
    method="raw"
    P = CI$p.raw
    low = CI$lower.raw
    upp = CI$upper.raw
    ylab="unadjusted p-value"
  }
  if (length(grep("bh", colnames(CI)[1]))==1) {
    method="BH"
    P = CI$p.bh
    low = CI$lower.bh
    upp = CI$upper.bh
    ylab="FDR-adjusted p-value (BH)"
  }
  if (length(grep("by", colnames(CI)[1]))==1) {
    method="BY"
    P = CI$p.by
    low = CI$lower.by
    upp = CI$upper.by
    ylab="FDR-adjusted p-value (BY)"
  }

  plot(logFC, -log10(P), cex.axis=1.5, cex.lab=1.5, xlab="logFC", ylab=ylab, axes=FALSE, xlim=xlim, ylim=ylim, type="n")
  axis(1, cex.axis=1.5)
  axis(2, cex.axis=1.5, seq(ylim[1], ylim[2], length.out=7), round((1/10)^(seq(ylim[1], ylim[2], length.out=7)), 4))
  box()
  points(logFC[sel], -log10(P)[sel], pch=15, col=8, cex=0.5)
  for (j in sel) {
    points(c(low[j], upp[j]), c(-log10(P[j]), -log10(P[j])), type="l", col=8)
    points(c(low[j], low[j]), c(-log10(P[j])-0.01, -log10(P[j])+0.01), type="l", col=1)
    points(c(upp[j], upp[j]), c(-log10(P[j])-0.01, -log10(P[j])+0.01), type="l", col=1)
  }
  abline(v=0)
  abline(h=-log10(alpha), col=2, lwd=2)
}
