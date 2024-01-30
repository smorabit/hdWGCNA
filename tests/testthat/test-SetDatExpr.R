# test that the output format is expected
test_that("SetDatExpr output is a data.frame", {

    data(test_seurat)

    test_seurat <- SetDatExpr(test_seurat, assay = 'RNA')
    datExpr <- GetDatExpr(test_seurat)
    expect_equal(class(datExpr), "data.frame")

})

# check that the datExpr has features for columns
test_that("features are columns", {

    data(test_seurat)

    test_seurat <- SetDatExpr(test_seurat, assay = 'RNA')
    datExpr <- GetDatExpr(test_seurat)

    check <- all(colnames(datExpr) %in% rownames(test_seurat))

    expect_equal(check, TRUE)

})

# check that the return_seurat flag works
test_that("return_seurat flag works", {

    data(test_seurat)

    datExpr <- SetDatExpr(test_seurat, assay = 'RNA', return_seurat=FALSE)
    
    expect_equal(class(datExpr), "data.frame")
})

# check that the use_metacells flag works
test_that("use_metacells flag works", {

    data(test_seurat)

    test_seurat <- SetDatExpr(test_seurat, assay = 'RNA', use_metacells=FALSE)
    datExpr <- GetDatExpr(test_seurat)

    check <- all(rownames(datExpr) %in% colnames(test_seurat))

    expect_equal(check, TRUE)
})


# test that group.by subsetting works as intended
test_that("group.by subsetting works", {

    data(test_seurat)

    test_seurat <- SetupForWGCNA(
        test_seurat,
        wgcna_name = 'test',
        features = rownames(test_seurat)
    )

    test_seurat <- MetacellsByGroups(
        test_seurat,
        group.by = c('Sample'),
        k = 5, 
        max_shared = 5,
        ident.group = 'Sample',
        target_metacells=10,
        min_cells=10
    )
    test_seurat <- NormalizeMetacells(test_seurat, verbose=FALSE)

    # selected groups to subset
    selected_groups <- c("C1", 'C2', 'C3')
    test_seurat <- SetDatExpr(
        test_seurat,
        group.by='Sample', 
        group_name = selected_groups,
        assay = 'RNA'
    )
    datExpr <- GetDatExpr(test_seurat)

    # check that the groups are in the rownames of the datExpr
    samples <- do.call(rbind, strsplit(rownames(datExpr), '_'))[,1]
    check <- all(selected_groups %in% samples)

    expect_equal(check, TRUE)

})
