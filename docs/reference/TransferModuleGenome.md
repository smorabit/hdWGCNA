# TransferModuleGenome

Takes a module table and a gene mapping table (like from biomart) with
gene names from two genomes in order to switch

## Usage

``` r
TransferModuleGenome(
  modules,
  gene_mapping,
  genome1_col = NULL,
  genome2_col = NULL
)
```

## Arguments

- modules:

  module table from GetModules function

- gene_mapping:

  table that has the gene names for both genomes

- genome1_col:

  the column in gene_mapping that has gene names for the genome
  currently present in modules

- genome2_col:

  the column in gene_mapping that has gene names for the genome to
  transfer to

## Examples

``` r
TransferModuleGenome
#> function (modules, gene_mapping, genome1_col = NULL, genome2_col = NULL) 
#> {
#>     if (is.null(genome1_col)) {
#>         genome1_col <- colnames(gene_mapping)[1]
#>     }
#>     if (is.null(genome2_col)) {
#>         genome2_col <- colnames(gene_mapping)[2]
#>     }
#>     if (!all(table(gene_mapping[[genome2_col]]) == 1)) {
#>         stop("There can only be one value for each feature in genome2_col in the gene_mapping table.")
#>     }
#>     gene_mapping <- gene_mapping[, c(genome1_col, genome2_col)]
#>     gene_list <- modules$gene_name
#>     gene_mapping <- gene_mapping[gene_mapping[, genome1_col] %in% 
#>         gene_list, ]
#>     modules <- subset(modules, gene_name %in% gene_mapping[, 
#>         genome1_col])
#>     gene_match <- match(gene_list, gene_mapping[, genome1_col])
#>     gene_mapping <- na.omit(gene_mapping[gene_match, ])
#>     print(head(gene_mapping))
#>     modules$gene_name <- gene_mapping[, genome2_col]
#>     rownames(modules) <- modules$gene_name
#>     modules
#> }
#> <bytecode: 0x55d89e4c9a30>
#> <environment: namespace:hdWGCNA>
```
