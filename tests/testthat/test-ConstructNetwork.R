# check that the network information gets added to the seurat object
test_that("Network info gets put into the Seurat object", {

    data(test_seurat)

    test_seurat <- ConstructNetwork(
        test_seurat,
        soft_power=5,
        verbose=FALSE
    )


    # # remove this in a sec:
    # test_seurat <- ModuleEigengenes(test_seurat)
    # test_seurat <- ModuleConnectivity(test_seurat)
    # head(GetModules(test_seurat))

    net <- GetNetworkData(test_seurat)

    # remove the TOM directory
    unlink('TOM', recursive=TRUE)
    #file.remove('ConsensusTOM-block.1.rda')

    expect_equal(is.null(net), FALSE)
})

# check that the TOM gets made
test_that("The TOM gets made", {

    data(test_seurat)

    wgcna_name <- 'test'
    
    test_seurat <- ConstructNetwork(
        test_seurat,
        soft_power=5,
        verbose=FALSE
    )
  
    check <- file.exists(paste0('TOM/', wgcna_name, '_TOM.rda'))

    # remove the TOM directory
    unlink('TOM', recursive=TRUE)
    #file.remove('ConsensusTOM-block.1.rda')

    expect_equal(check, TRUE)
})

# check that SetTOM / GetTOM round-trip works with a manually supplied matrix
test_that("SetTOM stores TOM in Seurat object and GetTOM retrieves it", {

    data(test_seurat)

    # build a small named dummy matrix
    genes <- GetWGCNAGenes(test_seurat)
    n <- length(genes)
    dummy_tom <- matrix(runif(n * n), nrow = n, ncol = n)
    rownames(dummy_tom) <- genes
    colnames(dummy_tom) <- genes

    test_seurat <- SetTOM(test_seurat, dummy_tom)
    retrieved <- GetTOM(test_seurat)

    expect_equal(retrieved, dummy_tom)
})

# check that GetTOM returns the in-object TOM without needing a file on disk
test_that("GetTOM returns in-object TOM without a file on disk", {

    data(test_seurat)

    genes <- GetWGCNAGenes(test_seurat)
    n <- length(genes)
    dummy_tom <- matrix(runif(n * n), nrow = n, ncol = n)
    rownames(dummy_tom) <- genes
    colnames(dummy_tom) <- genes

    test_seurat <- SetTOM(test_seurat, dummy_tom)

    # GetTOM should return from the object without touching any file
    expect_equal(is.null(GetTOM(test_seurat)), FALSE)
})

# check that store_tom_in_seurat=TRUE stores a correctly-shaped TOM in the object
test_that("ConstructNetwork with store_tom_in_seurat=TRUE stores TOM in Seurat object", {

    data(test_seurat)

    test_seurat <- ConstructNetwork(
        test_seurat,
        soft_power = 5,
        store_tom_in_seurat = TRUE,
        verbose = FALSE
    )

    tom <- GetTOM(test_seurat)
    modules <- GetModules(test_seurat)

    # remove the TOM directory
    unlink('TOM', recursive=TRUE)

    # TOM should be a square matrix with genes matching the modules table
    expect_equal(is.matrix(tom), TRUE)
    expect_equal(nrow(tom), nrow(modules))
    expect_equal(ncol(tom), nrow(modules))
    expect_equal(rownames(tom), modules$gene_name)
    expect_equal(colnames(tom), modules$gene_name)
})

# check that store_tom_in_seurat=TRUE with saveConsensusTOMs=FALSE throws an error
test_that("store_tom_in_seurat=TRUE with saveConsensusTOMs=FALSE throws an error", {

    data(test_seurat)

    expect_error(
        ConstructNetwork(
            test_seurat,
            soft_power = 5,
            store_tom_in_seurat = TRUE,
            saveConsensusTOMs = FALSE,
            verbose = FALSE
        )
    )
})