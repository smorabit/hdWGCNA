
# check that the soft power table has the expected values
test_that("Soft power results has expected power values", {

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
    test_seurat <- SetDatExpr(test_seurat, assay = 'RNA')

    selected_powers <- c(4,6,8)
    test_seurat <- TestSoftPowers(test_seurat, powers=selected_powers, corFnc='cor')

    pt <- GetPowerTable(test_seurat)
    check <- all(selected_powers %in% pt$Power )

    expect_equal(check, TRUE)
})

