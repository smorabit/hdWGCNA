# MotifScan

This function scans the promoter regions of protein-coding genes for
transcription factor (TF) motifs. It extracts promoter sequences using
an Ensembl database (`EnsDb`) and then searches these sequences for TF
binding motifs using position frequency matrices (PFMs). The Seurat
object is updated with a matrix of motif-gene matches, a list of target
genes for each TF, and additional motif information.

## Usage

``` r
MotifScan(seurat_obj, pfm, EnsDb, species_genome, wgcna_name = NULL)
```

## Arguments

- seurat_obj:

  A Seurat object that will be updated with motif scan results.

- pfm:

  A list of position frequency matrices (PFMs), such as those from the
  JASPAR2020 database.

- EnsDb:

  An Ensembl database object (e.g., `EnsDb.Hsapiens.v86` or
  `EnsDb.Mmusculus.v79`) containing gene and promoter annotations.

- species_genome:

  A character string specifying the genome version (e.g., "hg38" for
  human or "mm10" for mouse).

- wgcna_name:

  A character string specifying the name of the WGCNA experiment to
  associate with the motif data (optional). If NULL, the active WGCNA
  experiment in `seurat_obj@misc` will be used.

## Value

A modified Seurat object containing the results of the motif scan,
including:

- `seurat_obj@misc$motif_matrix`: A binary matrix indicating motif
  matches for each gene.

- `seurat_obj@misc$motif_info`: A data frame containing motif names,
  IDs, and the number of target genes.

- `seurat_obj@misc$motif_targets`: A list of genes targeted by each
  motif.

- `seurat_obj@misc$pfm`: The original PFMs used for the motif scan.

## Details

The `MotifScan` function performs the following steps:

- It extracts promoter regions (typically 2 kb upstream of the
  transcription start site) of protein-coding genes from the provided
  `EnsDb`.

- It uses the `motifmatchr` package to search these promoters for TF
  motifs, using the input PFMs.

- The function returns the Seurat object updated with several key pieces
  of information:

  - A motif-gene match matrix indicating the presence or absence of each
    motif in the promoter of each gene.

  - A list of target genes for each TF based on motif presence.

  - A summary of motifs, including the number of target genes for each
    TF.

- This data is stored in the `seurat_obj`'s metadata and can be used for
  downstream analysis, such as regulatory network inference.
