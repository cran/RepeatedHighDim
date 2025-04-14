#' COVID-19 Markers Dataset
#'
#' Marker data from COVID-19 analysis using CLR transformation, 5 PCs, and 1000 features.
#'
#' @name covid_markers
#' @keywords internal
NULL

#' Robust COVID-19 Markers Dataset
#'
#' Robust marker data from COVID-19 analysis using CLR transformation, 5 PCs, and 1000 features.
#'
#' @name robust_covid_markers
#' @keywords internal
NULL

#' Robust COVID-19 Markers (02 Trim) Dataset
#'
#' Robust marker data with 0.2 trimming from COVID-19 analysis using CLR transformation.
#'
#' @name robust_covid_markers_02trim
#' @keywords internal
NULL

#' Robust COVID-19 Markers (03 Trim) Dataset
#'
#' Robust marker data with 0.3 trimming from COVID-19 analysis using CLR transformation.
#'
#' @name robust_covid_markers_03trim
#' @keywords internal
NULL

#' Robust COVID-19 Dataset
#'
#' Robust COVID-19 data using CLR transformation, 5 PCs, and 1000 features.
#'
#' @name robust_covid_data
#' @keywords internal
NULL



#' COVID-19 CLR Transformation Marker Effects (Internal)
#'
#' Internal dataset containing marker gene effects using CLR transformation
#' (5 PCs, 1000 features) for evaluating trimming effects in scTrimClust.
#' Used by \code{\link{scTC_trim_effect}}.
#'
#' @format A data frame with marker genes as rows and the following columns:
#' \describe{
#'   \item{CellAnnotation}{Additional cell annotation information}
#'   \item{X}{Row identifier}
#'   \item{avg_log2FC}{Average log2 fold-change}
#'   \item{gene}{Gene identifier}
#'   \item{p_val}{Raw p-value}
#'   \item{p_val_adj}{Adjusted p-value}
#'   \item{pct.1}{Percentage of cells expressing the gene in cluster}
#'   \item{pct.2}{Percentage of cells expressing the gene in other clusters}
#'   \item{cluster}{Cluster assignment}
#' }
#' @name scTC_eff_clr
#' @keywords internal
NULL

#' Robust COVID-19 CLR Transformation Effects (Internal)
#'
#' Internal dataset containing robust marker gene effects using CLR transformation
#' (5 PCs, 1000 features) for evaluating trimming effects in scTrimClust.
#' Used by \code{\link{scTC_trim_effect}}.
#'
#' @format A data frame with the same structure as \code{scTC_eff_clr}
#' @name scTC_eff_clr_robust
#' @keywords internal
NULL

#' COVID-19 LogNormalized Marker Effects (Internal)
#'
#' Internal dataset containing marker gene effects using LogNormalization
#' (5 PCs, 1000 features) for evaluating trimming effects in scTrimClust.
#' Used by \code{\link{scTC_trim_effect}}.
#'
#' @format A data frame with the same structure as \code{scTC_eff_clr}
#' @name scTC_eff_log
#' @keywords internal
NULL

#' Robust COVID-19 LogNormalized Effects (Internal)
#'
#' Internal dataset containing robust marker gene effects using LogNormalization
#' (5 PCs, 1000 features) for evaluating trimming effects in scTrimClust.
#' Used by \code{\link{scTC_trim_effect}}.
#'
#' @format A data frame with the same structure as \code{scTC_eff_clr}
#' @name scTC_eff_log_robust
#' @keywords internal
NULL



#' ProcessedSingle-Cell Data
#'
#' A pre-processed Seurat object containing synthetic single cell data.
#'
#' @format A Seurat object with the following characteristics:
#' \describe{
#'   \item{Assays}{RNA assay with 200 features}
#'   \item{Cells}{150 single-cell samples}
#'   \item{Variable features}{100 most variable genes}
#'   \item{Layers}{counts (raw), data (normalized), scale.data (scaled)}
#'   \item{Dimensional reductions}{PCA, t-SNE}
#'   \item{Normalization}{LogNormalize with scale factor 10,000}
#' }
#'
#' @details
#' The object contains synthetic data with 150 cells for 3 cell types.
#' Processing steps match the Seurat tutorial and include:
#' \itemize{
#'   \item Identification of 100 most variable features
#'   \item Data scaling and centering
#'   \item PCA dimensional reduction (20 principal components)
#'   \item t-SNE dimensional reduction
#' }
#'
#' @name seurat_obj
#' @keywords internal
NULL
