
# test that the length of the gene list is expected:
test_that("List of selected genes is the correct length", {

    data(test_seurat)
    genes <- rownames(test_seurat)[1:200]

    test_seurat <- SetupForWGCNA(
        test_seurat, 
        wgcna_name = 'test',
        features = genes
    )
    wgcna_genes <- GetWGCNAGenes(test_seurat, 'test')

    expect_equal(length(wgcna_genes), length(genes))

})

# test that the length of the gene list is expected:
test_that("Duplicated genes are collapsed together", {

    data(test_seurat)
    genes <- rownames(test_seurat)[1:200]
    genes <- c(genes, genes)

    test_seurat <- SetupForWGCNA(
        test_seurat, 
        wgcna_name = 'test',
        features = genes
    )
    wgcna_genes <- GetWGCNAGenes(test_seurat, 'test')

    expect_equal(length(wgcna_genes), length(unique(genes)))

})

# test that selecting genes that aren't in the seurat object throws an error 
test_that("Selecting features that are not present throws an error", {

    data(test_seurat)
    genes <- paste0("Fake", seq(1:100))

    expect_error(SetupForWGCNA(
        test_seurat, 
        wgcna_name = 'test',
        features = genes
    ))

})

# test that the variable features are properly selected:
test_that("Variable features are selected when requested", {

    data(test_seurat)

    test_seurat <- SetupForWGCNA(
        test_seurat,
        wgcna_name = 'test',
        gene_select = 'variable'
    )
    wgcna_genes <- GetWGCNAGenes(test_seurat, 'test')

    expect_equal(all(wgcna_genes %in% VariableFeatures(test_seurat)), TRUE)

})