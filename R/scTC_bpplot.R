#' @title scTC_bpplot: Post-trim breakpoint heatmap for scTrimClust results
#'
#' @description Generates a heatmap showing the percentage overlap of marker genes between
#' the original (untrimmed) cluster markers and markers identified after trimming
#' at various percentages.
#'
#' @param ... Two or more data.frames/tibbles containing marker genes from 'FindAllMarkers'.
#' @param trim_percent_vector Numeric vector of trim percentages.
#' @param plot_title Character string for the heatmap title.
#' @param legend_title Character string for the legend title.
#' @param color Color palette for the heatmap.
#'
#' @return A ComplexHeatmap object
#'
#' @details scTC_bpplot compares marker genes between the original (untrimmed) clustering
#' and trimmed versions. For each cluster, it calculates what percentage of the
#' original markers are retained at each trim level. Clusters are ordered by the
#' number of markers in the original (untrimmed) results.
#'
#' At least two data.frames/tibbles containing marker genes from 'FindAllMarkers'
#' (from the 'Seurat' package) should be provided to \emph{...} as input.
#' The first data frame should be the original (untrimmed) results, followed by
#' trimmed results.
#'
#' \emph{trim_percent_vector} must be a numeric vector of trim percentages
#' corresponding to each input data frame (e.g., c(0,10,20,30,40) for untrimmed,
#' 10%, 20%, 30%, and 40% trimmed results). Must have same length as number of
#' input data.frames/tibbles.
#'
#'
#'
#' @examples
#' \dontrun{
#' scTC_bpplot(
#'   covid_markers = RepeatedHighDim:::covid_markers,
#'   robust_covid_markers = RepeatedHighDim:::robust_covid_markers,
#'   robust_covid_markers_02trim = RepeatedHighDim:::robust_covid_markers_02trim,
#'   robust_covid_markers_03trim = RepeatedHighDim:::robust_covid_markers_03trim,
#'   robust_covid_data = RepeatedHighDim:::robust_covid_data,
#'   trim_percent_vector = c(0, 10, 20, 30, 40),
#'   plot_title = "CLR, nPCs:5, nFeatures:1000",
#'   legend_title = "Percent markers of non-trimmed"
#' )
#' }
#'
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ComplexHeatmap Heatmap
#' @importFrom circlize colorRamp2
#' @export
scTC_bpplot <- function(...,
                        trim_percent_vector,
                        plot_title = "scTrimClust: Post-trim Breakpoint Heatmap",
                        legend_title = "Percent markes\nof non-trimmed",
                        color = brewer.pal(n = 11, name = "RdYlBu")) {
  de_list <- list(...)

  check_input_marker_matrices(de_list)
  check_trim_percent_vector(trim_percent_vector, length(de_list))

  if (length(color) != 11) {
    stop(
      "ERROR: 'color' must have exactly 11 values (for 0%, 10%, ..., 100%).\n",
      "You provided ", length(color), " colors: ", paste(color, collapse = ", "), "\n"
    )
  }




  de_list <- lapply(de_list, function(df) {
    df$cluster <- as.character(df$cluster)
    df
  })


  cells <- sort(unique(unlist(lapply(de_list, `[[`, "cluster"))))
  ncells <- length(cells)


  gene_sets <- lapply(cells, function(cl) {
    lapply(de_list, function(df) {
      genes <- df$gene[df$cluster == cl]
      if (length(genes) == 0) NA else genes
    })
  })


  nA <- vapply(gene_sets, function(x) {
    if (all(is.na(x[[1]]))) 0 else length(x[[1]])
  }, numeric(1))


  P <- vapply(gene_sets, function(sets) {
    base_set <- sets[[1]]
    if (all(is.na(base_set))) {
      return(rep(NA, length(sets)))
    }

    vapply(sets, function(s) {
      if (all(is.na(s))) {
        return(NA)
      }
      length(intersect(base_set, s)) / length(base_set)
    }, numeric(1))
  }, numeric(length(de_list)))


  P <- t(P) * 100



  colnames(P) <- paste0(trim_percent_vector, "%")
  rownames(P) <- paste0(cells, "(", nA, ")")



  P <- P[order(nA, decreasing = TRUE), , drop = FALSE]
  P[which(is.na(P))] <- 0




  colo <- colorRamp2(seq(0, 100, 10), color)
  hplot <- Heatmap(P,
    col = colo,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    name = legend_title,
    column_title = plot_title
  )

  return(hplot)
}


#' @noRd
check_input_marker_matrices <- function(de_list) {
  if (length(de_list) < 2) {
    stop("At least 2 data frames must be provided (original + at least one trimmed)")
  }


  required_cols <- c("gene", "cluster")
  for (i in seq_along(de_list)) {
    if (!all(required_cols %in% colnames(de_list[[i]]))) {
      stop(
        "Data frame ", i, " is missing required columns. ",
        "All inputs must contain: ", paste(required_cols, collapse = ", ")
      )
    }
  }

  return(TRUE)
}

#' @noRd
check_trim_percent_vector <- function(trim_percent_vector, n_datasets) {
  if (!is.numeric(trim_percent_vector)) {
    stop("'trim_percent_vector' must be numeric")
  }


  if (length(trim_percent_vector) != n_datasets) {
    stop(
      "'trim_percent_vector' length (", length(trim_percent_vector),
      ") must match number of input data frames (", n_datasets, ")"
    )
  }


  if (any(trim_percent_vector < 0 | trim_percent_vector > 100)) {
    stop("All trim percentages must be between 0 and 100")
  }


  if (trim_percent_vector[1] != 0) {
    stop("First value in trim_percent_vector must be 0 (untrimmed reference)")
  }



  return(TRUE)
}
