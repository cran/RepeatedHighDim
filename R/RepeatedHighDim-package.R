#' RepeatedHighDim Package
#'
#'@description
#'A comprehensive toolkit for repeated high-dimensional analysis.
#'
#'@details
#'The RepeatedHighDim-package is a collection of functions for the analysis of high-dimensional repeated measures data, e.g. from Omics experiments. It provides function for outlier detection, differential expression analysis, self-contained gene-set testing, and generation of correlated binary data.
#'
#'For more information and examples, please refer to the package documentation and the tutorial available at \url{https://software.klausjung-lab.de/}.
#'
#'
#' @section Functions:
#'
#'
#' This package includes the following functions:
#'
#'
#' \strong{B}:
#' \itemize{
#'   \item \code{\link{bag}}: Calculates the bag.
#' }
#'
#' \strong{D}:
#' \itemize{
#'   \item \code{\link{depmed}}: Calculates the depth median.
#' }
#'
#' \strong{F}:
#' \itemize{
#'   \item \code{\link{fc_ci}}: Calculates adjusted confidence intervals.
#'   \item \code{\link{fc_plot}}: Creates a volcano plot of adjusted confidence intervals.
#' }
#'
#' \strong{G}:
#' \itemize{
#'   \item \code{\link{GA_diagplot}}: Generates a diagnostic plot for comparing two correlation matrices.
#'   \item \code{\link{gem}}: Plots a gemstone to an interactive graphics device.
#'   \item \code{\link{GlobTestMissing}}: Detects global group effects.
#'   \item \code{\link{gridfun}}: Specifies a grid for calculating halfspace location depths.
#' }
#'
#' \strong{H}:
#' \itemize{
#'   \item \code{\link{hldepth}}: Calculates the halfspace location depth.
#' }
#'
#' \strong{I}:
#' \itemize{
#'   \item \code{\link{iter_matrix}}: Implements a genetic algorithm for generating correlated binary data.
#' }
#'
#' \strong{L}:
#' \itemize{
#'   \item \code{\link{loop}}: Calculates the fence and the loop.
#' }
#'
#' \strong{R}:
#' \itemize{
#'   \item \code{\link{RHighDim}}: Detects global group effects.
#'   \item \code{\link{rho_bounds}}: Calculates lower and upper bounds for pairwise correlations.
#'   \item \code{\link{rmvbinary_EP}}: Simulates correlated binary variables using the algorithm by Emrich and Piedmonte (1991).
#'   \item \code{\link{rmvbinary_QA}}: Simulates correlated binary variables using the algorithm by Qaqish (2003).
#' }
#'
#' \strong{S}:
#' \itemize{
#'   \item \code{\link{sequence_probs}}: Calculates probabilities for binary sequences.
#'   \item \code{\link{start_matrix}}: Sets up the start matrix.
#'   \item \code{\link{summary_RHD}}: Provides a summary of the RHighDim function.
#' }
#'
#' \strong{T}:
#' \itemize{
#'   \item \code{\link{TestStatSimple}}: Calculates the test statistic for RHighDim.
#'   \item \code{\link{TestStatSP}}: Calculates the test statistic for RHighDim.
#' }
#'
#'
#' @keywords bag gem outlier
#' @docType package
#' @name RepeatedHighDim
#' @aliases RepeatedHighDim-package
#'
#'
#'
#'
#'@seealso
#'For more information, please refer to the package's documentation and the tutorial: \url{https://software.klausjung-lab.de/}.
#'
#'
#'
#'@author
#'
#'
#' \strong{Maintainer}: Klaus Jung (\email{klaus.jung@tiho-hannover.de})
#'
#' \strong{Other contributors}:
#'
#' \itemize{
#'   \item Jochen Kruppa (\email{j.kruppa@hs-osnabrueck.de})
#'   \item Sergej Ruff (\email{Sergej.Ruff@tiho-hannover.de})
#' }
#'
#' If you have any questions, suggestions, or issues, please feel free to contact the maintainer, Klaus Jung (\email{klaus.jung@tiho-hannover.de}).
#'
#'
NULL

