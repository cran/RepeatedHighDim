#' scTrimClust: Cluster visualization with alpha hull-based outlier detection
#'
#' @description
#' Visualizes cell clusters in low-dimensional space (t-SNE, UMAP, etc.) and identifies/removes
#' potential outliers based on their distance from cluster alpha hulls.
#'
#' @param object A Seurat object containing dimensionality reduction results.
#' @param dims Integer vector of length 2 specifying which dimensions to plot (e.g., c(1, 2)).
#' @param cells Vector of cells to include (NULL uses all cells).
#' @param cols Vector of colors for clusters.
#' @param pt.size Point size for cells.
#' @param reduction Name of dimensionality reduction to use (e.g., "umap", "tsne").
#' @param group.by Metadata column to group cells by (default: 'ident' uses cluster IDs).
#' @param split.by Metadata column to split plots by (creates multiple facets).
#' @param shape.by Metadata column to determine point shapes.
#' @param order Vector specifying order to plot cells (affects z-ordering).
#' @param shuffle Logical to randomly shuffle plotting order.
#' @param seed Random seed for reproducibility when shuffle=TRUE.
#' @param label Logical to add cluster labels.
#' @param label.size Size of cluster labels.
#' @param label.color Color of cluster labels.
#' @param label.box Logical to add background box to labels.
#' @param repel Logical to use ggrepel for label placement.
#' @param alpha Transparency level for points (0-1).
#' @param stroke.size Size of point borders.
#' @param cells.highlight Specific cells to highlight.
#' @param cols.highlight Color(s) for highlighted cells.
#' @param sizes.highlight Size(s) for highlighted cells.
#' @param na.value Color for NA values.
#' @param ncol Number of columns for faceted plots.
#' @param combine Logical to combine multiple plots into one.
#' @param raster Logical to rasterize points (for large datasets).
#' @param raster.dpi Resolution for rasterized points.
#' @param add.alpha.hull Logical to compute and plot alpha hulls.
#' @param hull.alpha Alpha parameter for hull calculation. Higher values produce smoother hulls
#'        that encompass more cells (default: 2).
#' @param hull.color Color of the alpha hull lines (default: Null = same color as cluster points).
#' @param hull.size Thickness of the alpha hull lines (default: 0.5).
#' @param outlier.quantile Quantile threshold (0-1) for outlier detection based on hull distance.
#'        Cells with distances below this quantile are considered outliers (default: 0.4).
#' @param remove.outliers Logical - whether to remove outliers from the returned Seurat object
#'        (default: FALSE).
#' @param outlier.alpha Transparency level for outlier points (0-1; default: 0.1).
#' @param outlier.color Single color to use for all outlier points. If NULL, uses cluster colors.
#' @param outlier.colors A named vector of colors to be assigned to outliers.If NULL, uses cluster colors.
#' @param outline.color Color for the outline of points. If NULL, no outline is added.
#' @param outline.size Thickness of the outline around points (default: 0.5).
#' @param outline.alpha Transparency of the outline around points (default: 1).
#' @param outline.outliers Logical whether to add outlines to outlier points (default: FALSE).
#'
#' @return A list containing:
#' \itemize{
#'   \item \emph{plot}: ggplot object of the visualization with hulls and highlighted outliers
#'   \item \emph{object}: Modified Seurat object with outliers removed (if remove.outliers=TRUE)
#'   \item \emph{outlier_coords}: Dataframe containing coordinates of outlier cells, their IDs and cluster assignments
#'   \item \emph{hull_info}: List containing alpha hull geometries (if add.alpha.hull=TRUE)
#' }
#'
#' @examples
#' \dontrun{
#'
#' scTrimClust(RepeatedHighDim:::seurat_obj,reduction = 'tsne',
#' group.by = 'CellType',
#' hull.alpha = 2,
#' remove.outliers = FALSE,
#' outlier.quantile = 0.2,
#' outlier.alpha = 0.3,
#' outlier.color = "red",
#' pt.size = 5,
#' outline.color = "black",
#' outline.outliers = TRUE)
#'
#' # second example with custom outlier col per cluster
#'
#' scTrimClust(RepeatedHighDim:::seurat_obj,reduction = 'tsne',
#' group.by = 'CellType',
#' hull.alpha = 2,
#' remove.outliers = FALSE,
#' outlier.quantile = 0.2,
#' outlier.alpha = 0.3,
#' outlier.colors = c('TypeA'="black",
#' 'TypeB'='violet','TypeC' ='pink'),
#' pt.size = 5,
#' outline.color = "black",
#' outline.outliers = TRUE)$plot
#'
#' }
#'
#' @importFrom Seurat DimPlot
#' @importFrom alphahull ahull
#' @importFrom stats quantile
#' @importFrom rlang .data is_integerish
#' @export
scTrimClust <- function(
    object,
    dims = c(1, 2),
    cells = NULL,
    cols = NULL,
    pt.size = NULL,
    reduction = NULL,
    group.by = NULL,
    split.by = NULL,
    shape.by = NULL,
    order = NULL,
    shuffle = FALSE,
    seed = 1,
    label = FALSE,
    label.size = 4,
    label.color = 'black',
    label.box = FALSE,
    repel = FALSE,
    alpha = 1,
    stroke.size = NULL,
    cells.highlight = NULL,
    cols.highlight = '#DE2D26',
    sizes.highlight = 1,
    na.value = 'grey50',
    ncol = NULL,
    combine = TRUE,
    raster = NULL,
    raster.dpi = c(512, 512),
    add.alpha.hull = TRUE,
    hull.alpha = 2,
    hull.color = NULL,
    hull.size = 0.5,
    outlier.quantile = 0.4,
    remove.outliers = FALSE,
    outlier.alpha = 0.1,
    outlier.color = NULL,
    outlier.colors = NULL,
    outline.color = NULL,
    outline.size = 0.5,
    outline.alpha = 1,
    outline.outliers = FALSE
) {
  if (!rlang::is_integerish(dims, n = 2L, finite = TRUE) || !all(dims > 0L)) {
    stop("'dims' must be a two-length integer vector")
  }

  reduction <- reduction %||% SeuratObject::DefaultDimReduc(object = object)
  cells <- cells %||% Seurat::Cells(
    x = object,
    assay = Seurat::DefaultAssay(object = object[[reduction]])
  )
  dims <- paste0(Seurat::Key(object = object[[reduction]]), dims)
  group.by <- group.by %||% 'ident'

  data <- Seurat::FetchData(
    object = object,
    vars = c(dims, group.by),
    cells = cells,
    clean = 'project'
  )

  group.by <- colnames(data)[3:ncol(data)]
  for (group in group.by) {
    if (!is.factor(data[[group]])) {
      data[[group]] <- factor(data[[group]])
    }
  }

  tsne_coords <- data[, dims]
  clusters <- data[[group.by]]
  dims_names <- dims

  plots <- lapply(
    X = group.by,
    FUN = function(x) {
      current_clusters <- data[[x]]
      unique_clusters <- levels(current_clusters)
      non_outlier <- logical(nrow(data))
      outlier_coords <- data.frame(x = numeric(), y = numeric())
      nonoutliers_coords <- data.frame(x = numeric(), y = numeric())
      all_data_names <- c()

      for (cluster in unique_clusters) {
        cluster_idx <- which(current_clusters == cluster)

        if (length(cluster_idx) < 3) {
          non_outlier[cluster_idx] <- TRUE
          next
        }
        x_coords <- tsne_coords[cluster_idx, 1]
        y_coords <- tsne_coords[cluster_idx, 2]

        ahull <- tryCatch(
          alphahull::ahull(x_coords, y_coords, alpha = hull.alpha),
          error = function(e) NULL
        )
        if (is.null(ahull)) {
          non_outlier[cluster_idx] <- TRUE
          next
        }

        arcs <- ahull$arcs
        if (nrow(arcs) == 0) {
          non_outlier[cluster_idx] <- TRUE
          next
        }

        hull_vertices <- unique(arcs[, "end1"])
        hull_points <- ahull$xahull[hull_vertices, ]
        if (nrow(hull_points) == 0) next

        distances <- sapply(1:length(x_coords), function(i) {
          min(sqrt((x_coords[i] - hull_points[,1])^2 + (y_coords[i] - hull_points[,2])^2))
        })
        Q <- stats::quantile(distances, probs = outlier.quantile, na.rm = TRUE)
        nonoutliers <- which(distances > Q)
        outliersindex <- which(!(seq_along(cluster_idx) %in% nonoutliers))
        non_outlier[cluster_idx] <- seq_along(cluster_idx) %in% outliersindex

        if (length(nonoutliers) > 0) {
          nonoutlier_cell_ids <- rownames(data)[cluster_idx[nonoutliers]]
          nonoutlier_cell_types <- data[cluster_idx[nonoutliers], x]

          nonoutliers_coords <- rbind(
            nonoutliers_coords,
            data.frame(
              x = x_coords[nonoutliers],
              y = y_coords[nonoutliers],
              CellID = nonoutlier_cell_ids,
              CellType = nonoutlier_cell_types,
              stringsAsFactors = FALSE
            )
          )
        }

        if (length(outliersindex) > 0) {
          outlier_cell_ids <- rownames(data)[cluster_idx[outliersindex]]
          outlier_cell_types <- data[cluster_idx[outliersindex], x]

          outlier_coords <- rbind(
            outlier_coords,
            data.frame(
              x = x_coords[outliersindex],
              y = y_coords[outliersindex],
              CellID = outlier_cell_ids,
              CellType = outlier_cell_types,
              stringsAsFactors = FALSE
            )
          )
        }

        if (remove.outliers) {
          cluster_mask <- data[[3]] %in% cluster
          clust_data <- data[cluster_mask, ]
          data_name <- rownames(clust_data[nonoutliers, ])
          all_data_names <- c(all_data_names, data_name)
        }
      }

      if (remove.outliers) {
        data <- data[rownames(data) %in% all_data_names, ]
        object <- object[, which(colnames(object) %in% all_data_names)]
      }

      plot_data <- data[, c(dims, x, split.by, shape.by)]
      plot_data$alpha_value <- ifelse(rownames(plot_data) %in% outlier_coords$CellID,
                                      outlier.alpha,
                                      alpha)


      plot <- Seurat::SingleDimPlot(
        data = plot_data,
        dims = dims,
        col.by = x,
        cols = cols,
        pt.size = pt.size,
        shape.by = shape.by,
        order = order,
        alpha = plot_data$alpha_value,
        stroke.size = stroke.size,
        label = FALSE,
        cells.highlight = cells.highlight,
        cols.highlight = cols.highlight,
        sizes.highlight = sizes.highlight,
        na.value = na.value,
        raster = raster,
        raster.dpi = raster.dpi
      )


      if (add.alpha.hull) {
        ahull_list <- list()
        d_ahull_coords <- data.frame(x1 = NULL, y1 = NULL, x2 = NULL, y2 = NULL, cluster = NULL)

        for (cluster in unique_clusters) {
          cluster_idx <- which(current_clusters == cluster)
          if (length(cluster_idx) < 3) next
          x_coords <- tsne_coords[cluster_idx, 1]
          y_coords <- tsne_coords[cluster_idx, 2]

          ahull <- tryCatch(
            alphahull::ahull(x_coords, y_coords, alpha = hull.alpha),
            error = function(e) NULL
          )
          ahull_list[[cluster]] <- ahull
          if (is.null(ahull) || nrow(ahull$arcs) == 0) next

          arcs <- ahull$arcs
          d_ahull <- data.frame(
            x1 = ahull$xahull[arcs[,7], 1],
            y1 = ahull$xahull[arcs[,7], 2],
            x2 = ahull$xahull[arcs[,8], 1],
            y2 = ahull$xahull[arcs[,8], 2],
            cluster = cluster,
            stringsAsFactors = FALSE
          )

          d_ahull_coords <- rbind(d_ahull_coords, d_ahull)
        }


        if (nrow(d_ahull_coords) > 0) {
          if (is.null(hull.color)) {

            if (is.null(cols)) {
              cols <- scales::hue_pal()(length(unique_clusters))
            }
            color_mapping <- setNames(cols, unique_clusters)

            plot <- plot +
              ggplot2::geom_segment(
                data = d_ahull_coords,
                ggplot2::aes(x = .data$x1, y = .data$y1,
                             xend = .data$x2, yend = .data$y2,
                             color = .data$cluster),
                size = hull.size,
                show.legend = FALSE
              ) +
              ggplot2::scale_color_manual(values = color_mapping)
          } else {

            plot <- plot +
              ggplot2::geom_segment(
                data = d_ahull_coords,
                ggplot2::aes(x = .data$x1, y = .data$y1,
                             xend = .data$x2, yend = .data$y2),
                color = hull.color,
                size = hull.size,
                show.legend = FALSE
              )
          }


          plot$layers <- c(plot$layers[[length(plot$layers)]], plot$layers[-length(plot$layers)])
        }
      }


      if (!is.null(outlier.color)) {
        outlier_data <- plot_data[rownames(plot_data) %in% outlier_coords$CellID, ]
        plot <- plot +
          ggplot2::geom_point(
            data = outlier_data,
            ggplot2::aes(x = .data[[dims[1]]], y = .data[[dims[2]]]),
            color = outlier.color,
            size = pt.size %||% 1,
            alpha = outlier.alpha
          )
      }

      if (!is.null(outline.color)) {
        if (outline.outliers) {
          outline_data <- plot_data
        } else {
          outline_data <- plot_data[!(rownames(plot_data) %in% outlier_coords$CellID), ]
        }

        outline_layer <- ggplot2::geom_point(
          data = outline_data,
          ggplot2::aes(x = .data[[dims[1]]], y = .data[[dims[2]]]),
          color = outline.color,
          size = pt.size %||% 1,
          stroke = outline.size,
          alpha = outline.alpha,
          shape = 1
        )
        plot$layers <- c(list(outline_layer), plot$layers)
      }

      final_results <- list(
        plot = plot,
        object = object,
        nonoutliers_coords = nonoutliers_coords,
        outlier_coords = outlier_coords
      )

      if (add.alpha.hull) final_results$d_ahull_coords <- d_ahull_coords
      if (add.alpha.hull) final_results$ahull_list <- ahull_list

      return(final_results)
    }
  )

  if (add.alpha.hull) d_ahull_coords <- plots[[1]]$d_ahull_coords
  if (add.alpha.hull) ahull_list <- plots[[1]]$ahull_list
  object <- plots[[1]]$object
  nonoutliers_coords <- plots[[1]]$nonoutliers_coords
  outlier_coords <- plots[[1]]$outlier_coords
  plots <- plots[[1]]$plot

  if (combine) {
    plots <- patchwork::wrap_plots(plots, ncol = ncol)
  }

  if(!is.null(outlier.colors)){


    outlier_coords$Colc <- outlier.colors[match(outlier_coords[,4], names(outlier.colors))]


    plots <- plots +
      ggplot2::geom_point(
        data = outlier_coords,
        ggplot2::aes(x = .data$x, y = .data$y),
        color = 'white',
        size = pt.size %||% 1,
        alpha = 1
      )+
      ggplot2::geom_point(
        data = outlier_coords,
        ggplot2::aes(x = .data$x, y = .data$y),
        color = outlier_coords$Colc,
        size = pt.size %||% 1,
        alpha = outlier.alpha
      )
  }

  exported_results <- list(
    plot = plots,
    object = object,
    nonoutliers_coords = nonoutliers_coords,
    outlier_coords = outlier_coords
  )

  if (add.alpha.hull) exported_results$d_ahull_coords <- d_ahull_coords
  if (add.alpha.hull) exported_results$ahull_list <- ahull_list

  return(exported_results)
}
