
# check that we have MEs calculated for all the modules
test_that("MEs calculated for all modules", {

    data(test_seurat)

    # get module info
    modules <- GetModules(test_seurat)
    mods <- levels(modules$module)

    # compute MEs
    test_seurat <- ModuleEigengenes(test_seurat, verbose=FALSE)    
    MEs <- GetMEs(test_seurat, harmonized=FALSE)

    # check
    check <- all(mods %in% colnames(MEs) )
    expect_equal(check, TRUE)
})


# check that harmonization works
test_that("ME Harmonization works", {

    data(test_seurat)

    test_seurat <- ModuleEigengenes(test_seurat, group.by = 'specimenID', verbose=FALSE)
    hMEs <- GetMEs(test_seurat, harmonized=TRUE)
    MEs <- GetMEs(test_seurat, harmonized=FALSE)

    check <- !is.logical(all.equal(hMEs, MEs)) & all.equal(dim(hMEs), dim(MEs))
    expect_equal(check, TRUE)
})
