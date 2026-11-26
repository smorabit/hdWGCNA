#' Simulated single-cell dataset with planted co-expression modules
#'
#' A small simulated scRNA-seq dataset for examples and tests. Counts were generated with
#' 'splatter' and five gene co-expression modules were then planted into it, so the true
#' module memberships are known (see \code{\link{hdw_sim_modules}}). The generating
#' script is \code{data-raw/hdw_sim.R} in the package source repository.
#'
#' The dataset contains three cell types (A, B, C) from eight samples, four
#' "control" and four "disease". The planted modules are:
#' \describe{
#'   \item{M1}{active in all cell types}
#'   \item{M2}{specific to cell type A}
#'   \item{M3}{specific to cell type B}
#'   \item{M4}{specific to cell type A and up-regulated in disease samples (a true
#'     differential module)}
#'   \item{M5}{specific to cell type C and only active in control samples (not preserved
#'     in disease)}
#' }
#'
#' @format A Seurat object with 1,000 genes and 4,000 cells, containing raw counts in the
#'   "RNA" assay and the following cell metadata:
#' \describe{
#'   \item{cell_type}{factor, simulated cell type (A, B, C)}
#'   \item{Sample}{factor, sample of origin (S01-S08)}
#'   \item{condition}{factor, "control" (S01-S04) or "disease" (S05-S08)}
#' }
#' @seealso \code{\link{hdw_sim_modules}}
#' @source Simulated; see \code{data-raw/hdw_sim.R}.
"hdw_sim"

#' Ground-truth modules of the simulated dataset
#'
#' The genes planted into each co-expression module of \code{\link{hdw_sim}}.
#'
#' @format A data frame with one row per planted gene (250 rows) and columns:
#' \describe{
#'   \item{module}{planted module (M1-M5)}
#'   \item{gene_name}{gene}
#'   \item{loading}{strength of the gene's coupling to the module's latent factor (0.5-1);
#'     genes with higher loadings are the true hub genes}
#'   \item{cell_types}{comma-separated cell types in which the module is active}
#'   \item{conditions}{comma-separated conditions in which the module is active}
#'   \item{disease_shift}{shift of the module's latent activity in disease samples}
#' }
#' @seealso \code{\link{hdw_sim}}
#' @source Simulated; see \code{data-raw/hdw_sim.R}.
"hdw_sim_modules"

#' Small Seurat object used by the legacy unit tests
#'
#' A small single-cell dataset used by the original hdWGCNA unit tests. It will be
#' replaced by \code{\link{hdw_sim}}.
#'
#' @format A Seurat object.
#' @keywords internal
"test_seurat"
