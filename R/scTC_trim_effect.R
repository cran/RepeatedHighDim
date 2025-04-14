#' @title scTC_trim_effect: Compare scTrimClust trimming  against default Seurat analysis
#'
#' @description
#' Visualizes the impact of scTrimClust's trimming by comparing gene sets between:
#' 1) Default Seurat analysis (no trimming)
#' 2) scTrimClust post-trimming results
#'
#'
#'
#' @param method_pairs A named list of method comparisons. Each element should be a list
#'        with two components (\emph{data1 (untrimmed)} and \emph{data2 (trimmed)}) containing data frames with:
#'        \itemize{
#'          \item \emph{cluster}: Cluster identifiers
#'          \item \emph{gene}: Gene identifiers
#'        }
#' @param method_colors Named vector of colors for method annotations. Names should match
#'        the names in \emph{method_pairs}.
#' @param set_colors Named vector of colors for set annotations (S1-S3). Default:
#'        \emph{c("S1:standard", "S2:intersect", "S3:trimmed")} with grey colors.
#' @param heatmap_color_palette Color mapping function for heatmap. Default:
#'        \emph{colorRamp2(seq(0, 100, 1), heat.colors(101, rev = TRUE))}.
#' @param column_title Main title for the heatmap columns.
#' @param row_names_side Side for row names ("left" or "right"). Default: "right".
#' @param legend_name Title for the heatmap legend. Default: "No. of markers".
#' @param row_names_gp Graphics parameters for row names. Default: 10.
#' @param column_title_gp Graphics parameters for column title. Default: 12.
#'
#' @return A \emph{Heatmap} object from the ComplexHeatmap package.
#'
#' @details
#' scTC_trim_effect creates a heatmap showing the percentage differences in gene sets between method pairs
#' across clusters.
#' The heatmap shows three components for each method comparison:
#' \itemize{
#'   \item Column 1-3: Unique to method1 (untrimmed), Intersection, Unique to method2 (trimmed)
#'   \item Rows represent (cell) clusters with counts from first method in parentheses
#'   \item Columns are split by method pairs
#' }
#'
#' @examples
#' \dontrun{
#'
#' method_pairs <- list(
#'   CLR = list(
#'     data1 = RepeatedHighDim:::scTC_eff_clr,
#'     data2 = RepeatedHighDim:::scTC_eff_clr_robust
#'   ),
#'   LogNorm = list(
#'     data1 = RepeatedHighDim:::scTC_eff_log,
#'     data2 = RepeatedHighDim:::scTC_eff_log_robust
#'   )
#' )
#'
#' method_colors <- setNames(grey.colors(2), c("CLR", "LogNorm"))
#'
#' scTC_trim_effect(
#'   method_pairs = method_pairs,
#'   method_colors = method_colors,
#'   column_title = "nPCs:5, nFeatures:1000"
#' )
#'
#' set_colors <- grey.colors(3)
#' names(set_colors) <- c("S1:standard", "S2:intersect", "S3:trimmed")
#'
#' scTC_trim_effect(
#'   method_pairs = method_pairs,
#'   method_colors = method_colors,
#'   set_colors = setNames(c("blue", "green", "red"), names(set_colors)),
#'   heatmap_color_palette = colorRamp2(c(0, 50, 100), c("white", "pink", "purple")),
#'   column_title = "Custom Color Example"
#' )
#' }
#'
#' @importFrom ComplexHeatmap Heatmap HeatmapAnnotation
#' @importFrom circlize colorRamp2
#' @importFrom stats setNames
#' @importFrom grDevices heat.colors
#' @export
scTC_trim_effect <- function(
    method_pairs,
    method_colors,
    set_colors = setNames(c("#4D4D4D", "#AEAEAE", "#E6E6E6"), c("S1:standard", "S2:intersect", "S3:trimmed")),
    heatmap_color_palette = colorRamp2(seq(0, 100, 1), heat.colors(101, rev = TRUE)),
    column_title = "", row_names_side = "right", legend_name = "No. of\nmarkers",row_names_gp = 10,
    column_title_gp = 12) {
  cells <- sort(unique(unlist(lapply(method_pairs, function(m) {
    unique(c(m$data1$cluster, m$data2$cluster))
  }))))

  ncells <- length(cells)



  P_list <- list()
  counts_matrix <- matrix(NA, nrow = ncells, ncol = length(method_pairs))
  colnames(counts_matrix) <- names(method_pairs)


  for (idx in seq_along(method_pairs)) {
    method_name <- names(method_pairs)[idx]
    m <- method_pairs[[method_name]]

    P_method <- matrix(NA, nrow = 3, ncol = ncells)
    method_counts <- numeric(ncells)

    for (i in 1:ncells) {
      cell <- cells[i]

      set1 <- m$data1$gene[m$data1$cluster == cell]
      set2 <- m$data2$gene[m$data2$cluster == cell]


      P_method[1, i] <- length(setdiff(set1, set2))
      P_method[2, i] <- length(intersect(set1, set2))
      P_method[3, i] <- length(setdiff(set2, set1))


      method_counts[i] <- length(set1)
    }


    P_method <- 100 * prop.table(P_method, 2)
    P_list[[method_name]] <- P_method
    counts_matrix[, idx] <- method_counts
  }


  P_combined <- do.call(rbind, P_list)
  P <- t(P_combined)


  cname <- paste0(cells, " (", apply(counts_matrix, 1, paste, collapse = ", "), ")")
  o <- order(rowMeans(counts_matrix), decreasing = TRUE)
  P <- P[o, ]
  rownames(P) <- cname[o]


  method_factor <- rep(names(method_pairs), each = 3)
  set_factor <- rep(names(set_colors), length(method_pairs))

  ha <- HeatmapAnnotation(
    Method = method_factor,
    Set = set_factor,
    col = list(Method = method_colors, Set = set_colors)
  )


  heatplot <- Heatmap(P,
    col = heatmap_color_palette,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    name = legend_name,
    top_annotation = ha,
    column_title = column_title,
    column_split = rep(names(method_pairs), each = 3),
    row_names_side = row_names_side,
    row_names_gp = grid::gpar(fontsize = row_names_gp),
    column_title_gp = grid::gpar(fontsize = column_title_gp)
  )

  return(heatplot)
}
