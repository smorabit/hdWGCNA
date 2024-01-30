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