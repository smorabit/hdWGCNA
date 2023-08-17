# check that the network information gets added to the seurat object
test_that("Network info gets put into the Seurat object", {

    data(test_seurat)

    test_seurat <- SetupForWGCNA(
        test_seurat,
        wgcna_name = 'test',
        features = rownames(test_seurat)
    )

    test_seurat <- MetacellsByGroups(
        test_seurat,
        group.by = c('cell_type'),
        k = 5, 
        max_shared = 5,
        ident.group = 'cell_type',
        target_metacells=50,
        min_cells=10
    )
    test_seurat <- NormalizeMetacells(test_seurat, verbose=FALSE)
    test_seurat <- SetDatExpr(test_seurat, assay = 'RNA', verbose=FALSE)
    test_seurat <- ConstructNetwork(
        test_seurat,
        soft_power=5,
        verbose=FALSE
    )
  
    net <- GetNetworkData(test_seurat)

    # remove the TOM directory
    unlink('TOM', recursive=TRUE)

    expect_equal(is.null(net), FALSE)
})

# check that the TOM gets made
test_that("The TOM gets made", {

    data(test_seurat)

    wgcna_name <- 'test'
    test_seurat <- SetupForWGCNA(
        test_seurat,
        wgcna_name = wgcna_name,
        features = rownames(test_seurat)
    )

    test_seurat <- MetacellsByGroups(
        test_seurat,
        group.by = c('cell_type'),
        k = 5, 
        max_shared = 5,
        ident.group = 'cell_type',
        target_metacells=50,
        min_cells=10
    )
    test_seurat <- NormalizeMetacells(test_seurat, verbose=FALSE)
    test_seurat <- SetDatExpr(test_seurat, assay = 'RNA', verbose=FALSE)
    test_seurat <- ConstructNetwork(
        test_seurat,
        soft_power=5,
        verbose=FALSE
    )
  
    check <- file.exists(paste0('TOM/', wgcna_name, '_TOM.rda'))

    # remove the TOM directory
    unlink('TOM', recursive=TRUE)

    expect_equal(check, TRUE)
})