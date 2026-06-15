# NormalizeCounts

Perform per-sample normalization of a count assay stored in a
SummarizedExperiment and add the normalized matrix as a new assay.

## Usage

``` r
NormalizeCounts(
  se,
  method = c("CPM", "logCPM", "logNorm", "VST", "rlog"),
  assay_name = "counts",
  new_assay_name = NULL,
  pseudocount = 1,
  ...
)
```

## Arguments

- se:

  A SummarizedExperiment containing a raw counts assay.

- method:

  Character; one of "CPM", "logCPM", "logNorm", "VST", "rlog".

  - "CPM": counts per million.

  - "logCPM": log2(CPM + pseudocount).

  - "logNorm": Seurat-style log1p(counts / size_factor \* 1e4) where
    size_factor = colSums(counts) / median(colSums(counts)).

  - "VST", "rlog": variance-stabilizing transform or rlog via DESeq2.

- assay_name:

  Character scalar; name of the assay in `se` to normalize (default:
  "counts").

- new_assay_name:

  Character or NULL; name to assign the normalized assay. If NULL,
  defaults to the chosen `method`.

- pseudocount:

  Numeric scalar added to CPM before log2 in "logCPM" (default: 1).

- ...:

  Additional arguments forwarded to DESeq2::vst or DESeq2::rlog when
  `method` is "VST" or "rlog".

## Value

A SummarizedExperiment identical to `se` but with a new assay named
`new_assay_name` containing the normalized matrix.

## Details

The function checks that `se` is a SummarizedExperiment and that
`assay_name` exists and is a (possibly sparse) matrix. For "VST" and
"rlog", DESeq2 must be installed; the function converts the assay to a
dense matrix and constructs a DESeqDataSet with design ~ 1 before
applying the transform. Errors are raised for invalid inputs.

## See also

DESeq2::vst, DESeq2::rlog

## Examples

``` r
# Basic usage (assuming `se` is a SummarizedExperiment with a "counts" assay)
# se_norm <- NormalizeSE(se, method = "logCPM")
```
