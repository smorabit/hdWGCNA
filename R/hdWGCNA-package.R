#' @keywords internal
"_PACKAGE"

# Package-level imports.
#
# All functions from packages in Imports that are called without `pkg::` are listed here.
# Packages in Suggests are never imported: call them as `pkg::fn()` after
# `check_suggests()` (see STYLE.md).
#
# Some imports deliberately mask base/stats functions:
# - Matrix generics (rowSums, colSums, colMeans, t, diag, which, mean, as.matrix) so that
#   these work on sparse matrices
# - WGCNA::cor (faster, supports WGCNA's correlation options); WGCNA::bicor is imported
#   because pickSoftThreshold() looks up corFnc = 'bicor' by name from the caller
# - fastcluster::hclust (faster drop-in replacement for stats::hclust, as used by WGCNA)

#' @importFrom dplyr %>% all_of arrange desc distinct group_by mutate mutate_at n_distinct pull
#'   sample_n select slice_max slice_min summarise summarise_all top_n vars
#' @importFrom fastcluster hclust
#' @importFrom ggplot2 aes aes_string annotate arrow coord_equal coord_fixed coord_flip
#'   element_blank element_line element_rect element_text facet_wrap geom_bar geom_hline
#'   geom_point geom_rect geom_segment geom_smooth geom_text geom_tile geom_vline ggplot
#'   ggplot_build ggtitle guide_colorbar guides labs margin scale_color_gradient
#'   scale_color_gradient2 scale_color_manual scale_colour_manual scale_fill_gradient
#'   scale_fill_gradient2 scale_fill_manual scale_shape_manual scale_size scale_x_continuous
#'   scale_y_continuous scale_y_discrete scale_y_reverse sec_axis theme theme_void unit xlab
#'   xlim ylab ylim
#' @importFrom ggraph circle geom_edge_fan geom_edge_link geom_edge_loop geom_node_label
#'   geom_node_point ggraph scale_edge_colour_gradient2
#' @importFrom ggrepel geom_text_repel
#' @importFrom grDevices adjustcolor dev.off pdf
#' @importFrom igraph V V<- vcount
#' @importFrom magrittr %<>%
#' @importFrom Matrix as.matrix colMeans colSums diag mean rowSums t which
#' @importFrom methods as is
#' @importFrom patchwork plot_annotation plot_layout wrap_plots
#' @importFrom proxy dist
#' @importFrom purrr map map_int
#' @importFrom scales alpha
#' @importFrom Seurat DoHeatmap FindMarkers NoLegend RotatedAxis ScaleData
#' @importFrom SeuratObject Assays CreateSeuratObject DefaultAssay GetAssayData Idents Idents<-
#'   Images SetAssayData VariableFeatures
#' @importFrom stats median na.omit p.adjust quantile reorder setNames
#' @importFrom SummarizedExperiment assay assay<- assays colData colData<- SummarizedExperiment
#' @importFrom tidygraph centrality_degree
#' @importFrom utils combn head packageVersion setTxtProgressBar
#' @importFrom WGCNA bicor checkSets cor
NULL
