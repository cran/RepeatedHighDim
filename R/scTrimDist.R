#' ScTrimDist: Trim extreme cells based on kNN distance within cell types
#'
#' @description
#' Identifies and removes extreme (outlier) cells within each cell type or cluster
#' based on k-nearest neighbour (kNN) distances computed in the normalized
#' high-dimensional gene expression space. Cells located in sparsely populated
#' regions at the periphery of clusters are excluded prior to downstream analyses.
#'
#' @details
#' For each cell type (or cluster), a kNN search is performed using the normalized
#' gene expression matrix obtained from a standard Seurat preprocessing workflow.
#' For a given cell \eqn{i} in cluster \eqn{k}, the Euclidean distances
#' \eqn{D_{(j,i)}^k} to its \eqn{j = 1, \ldots, K} nearest neighbours are computed.
#'
#' The minimum distance
#' \deqn{
#' \min D_i^k = \min_{j = 1, \ldots, K} D_{(j,i)}^k
#' }
#' is used as a measure of local neighbourhood density. Cells with large minimum
#' distances are interpreted as extreme or non-representative cells.
#'
#' A fraction \eqn{\alpha} (specified via \code{keep_frac}) of the most extreme cells
#' is removed per cluster, defined as cells with
#' \deqn{
#' \min D_i^k > Q_{1 - \alpha}
#' }
#' where \eqn{Q_{1 - \alpha}} is the \eqn{(1 - \alpha)} quantile of the minimum
#' kNN distance distribution within the cluster.
#'
#' After trimming, the remaining cells are re-normalized and reprocessed using
#' standard Seurat workflows. Cell type annotations are assigned using a
#' **precomputed SingleR result** supplied by the user, and cluster-specific
#' marker genes are identified.
#'
#' @param seurat_obj A \code{Seurat} object containing single-cell expression data.
#' @param celltype_col Character scalar specifying the column in
#'   \code{seurat_obj@meta.data} defining cell types or clusters.
#' @param knn_k Integer specifying the number of nearest neighbours.
#' @param keep_frac Numeric in (0,1) specifying the fraction of most extreme cells
#'   to remove per cell type.
#' @param normalization_method Normalization method passed to
#'   \code{Seurat::NormalizeData}.
#' @param nfeatures Number of variable features selected.
#' @param assay Assay used for expression data extraction.
#' @param npcs Number of principal components used downstream.
#' @param resolution Clustering resolution for \code{FindClusters}.
#' @param log2FC_filter Minimum log2 fold-change threshold for marker filtering.
#'   If \code{NULL}, no filtering is applied.
#' @param pred A \code{SingleR} result object. Row names must correspond to cell
#'   barcodes; \code{pred$labels} is used for annotation.
#' @param verbose Logical indicating whether progress messages are printed.
#'
#' @return A named list containing:
#' \itemize{
#'   \item \code{plot_outliers}: ggplot showing t-SNE with outliers highlighted.
#'   \item \code{trimmed_object}: Seurat object after trimming and reprocessing.
#'   \item \code{all_markers}: Data frame of marker genes.
#'   \item \code{knn_res}: List of kNN results per cell type.
#' }
#'
#' @import Seurat
#' @importFrom FNN get.knn
#' @importFrom ggplot2 ggplot aes geom_point geom_histogram geom_vline labs theme_bw scale_color_manual scale_fill_manual
#' @importFrom dplyr filter
#' @importFrom scales hue_pal
#'
scTrimDist <- function(
    seurat_obj,
    celltype_col,
    knn_k = 30,
    keep_frac = 0.05,
    normalization_method = "LogNormalize",
    nfeatures = 2000,
    assay = "RNA",
    npcs = 20,
    resolution = 0.5,
    log2FC_filter = 1,
    pred,
    verbose = TRUE
) {


  stopifnot(inherits(seurat_obj, "Seurat"))
  stopifnot(celltype_col %in% colnames(seurat_obj@meta.data))
  stopifnot(is.numeric(keep_frac), keep_frac > 0, keep_frac < 1)
  stopifnot(!is.null(pred$labels))
  stopifnot(!is.null(rownames(pred)))

  meta <- seurat_obj@meta.data
  celltypes <- unique(meta[[celltype_col]])


  expr_mat <- t(as.matrix(
    GetAssayData(seurat_obj, assay = assay, layer = "data")
  ))


  if (!"tsne" %in% names(seurat_obj@reductions)) {
    stop("t-SNE reduction not found in seurat_obj.")
  }
  tsne_coords <- Embeddings(seurat_obj, reduction = "tsne")

  outlier_flag <- rep(FALSE, nrow(expr_mat))
  names(outlier_flag) <- rownames(expr_mat)
  knn_res_list <- vector("list", length(celltypes))
  names(knn_res_list) <- celltypes


  for (ct in celltypes) {

    cell_idx <- which(meta[[celltype_col]] == ct)

    if (length(cell_idx) <= 3) {
      if (verbose) {
        message("Skipping cell type ", ct, " (n = ", length(cell_idx), ")")
      }
      next
    }

    this_k <- min(knn_k, length(cell_idx) - 1)
    ct_expr <- expr_mat[cell_idx, , drop = FALSE]

    knn_res <- FNN::get.knn(ct_expr, k = this_k)
    min_dists <- apply(knn_res$nn.dist, 1, min, na.rm = TRUE)

    knn_res_list[[ct]] <- knn_res
    threshold <- quantile(min_dists, probs = 1 - keep_frac, na.rm = TRUE)

    outliers <- cell_idx[min_dists > threshold]
    outlier_flag[outliers] <- TRUE


    hist_plot <- ggplot(
      data.frame(dist = min_dists),
      aes(x = .data$dist)
    ) +
      geom_histogram(bins = 50, fill = "steelblue", color = "white") +
      geom_vline(xintercept = threshold, color = "red", linetype = "dashed") +
      labs(
        title = paste("Min kNN distance:", ct),
        x = "Min kNN distance",
        y = "Cell count"
      ) +
      theme_bw()




  }


  plot_df <- data.frame(
    tsne1 = tsne_coords[, 1],
    tsne2 = tsne_coords[, 2],
    celltype = meta[[celltype_col]],
    status = ifelse(outlier_flag, "Outlier", "Normal")
  )

  p_outliers <- ggplot(plot_df, aes(.data$tsne1, .data$tsne2)) +
    geom_point(
      data = subset(plot_df, .data$status == "Normal"),
      aes(color = .data$celltype),
      alpha = 0.3,
      size = 1
    ) +
    geom_point(
      data = subset(plot_df, .data$status == "Outlier"),
      aes(fill = .data$celltype),
      shape = 21,
      color = "black",
      size = 1.5
    ) +
    scale_color_manual(values = hue_pal()(length(celltypes))) +
    scale_fill_manual(values = hue_pal()(length(celltypes))) +
    labs(x = "t-SNE 1", y = "t-SNE 2") +
    theme_bw()

  if (verbose) {
    message("Removed ", sum(outlier_flag), " cells (",
            round(mean(outlier_flag) * 100, 2), "%).")
  }


  trimmed <- subset(seurat_obj, cells = names(outlier_flag)[!outlier_flag])

  trimmed <- NormalizeData(trimmed,
                           normalization.method = normalization_method,
                           verbose = verbose)
  trimmed <- FindVariableFeatures(trimmed,
                                  nfeatures = nfeatures,
                                  verbose = verbose)
  trimmed <- ScaleData(trimmed, verbose = verbose)
  trimmed <- RunPCA(trimmed, npcs = npcs, verbose = verbose)
  trimmed <- RunTSNE(trimmed, dims = 1:npcs, verbose = verbose)
  trimmed <- RunUMAP(trimmed, dims = 1:npcs, verbose = verbose)
  trimmed <- FindNeighbors(trimmed, dims = 1:npcs, verbose = verbose)
  trimmed <- FindClusters(trimmed,
                          resolution = resolution,
                          verbose = verbose)


  common_cells <- intersect(
    rownames(trimmed@meta.data),
    rownames(pred)
  )

  trimmed$CellAnnotation <- NA
  trimmed$CellAnnotation[common_cells] <-
    pred$labels[match(common_cells, rownames(pred))]

  Idents(trimmed) <- trimmed$CellAnnotation


  markers <- FindAllMarkers(trimmed,
                            only.pos = TRUE,
                            verbose = verbose)

  if (!is.null(log2FC_filter)) {
    markers <- dplyr::filter(
      markers,
      .data$avg_log2FC > log2FC_filter,
      .data$p_val_adj < 0.05
    )
  }


  return(list(
    plot_outliers = p_outliers,
    trimmed_object = trimmed,
    all_markers = markers,
    knn_res = knn_res_list
  ))
}
