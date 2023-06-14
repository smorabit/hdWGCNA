
# test that the output format is Seurat 
test_that("Metacell output is Seurat format", {

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

    m_obj <- GetMetacellObject(test_seurat)

    expect_equal("Seurat" %in% class(m_obj), TRUE)

})

# test that the genes are the same in the metacell and the main seurat object 
test_that("Metacell features are same as Seurat object",{

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

    m_obj <- GetMetacellObject(test_seurat)

    expect_equal(rownames(m_obj), rownames(test_seurat))

})


# test that the meta-data is getting passed 
test_that("Meta-data gets passed to the metacell object",{

    data(test_seurat)

    test_seurat <- SetupForWGCNA(
        test_seurat,
        wgcna_name = 'test',
        features = rownames(test_seurat)
    )

    test_seurat <- MetacellsByGroups(
        test_seurat,
        group.by = c('cell_type', 'annotation'),
        k = 5, 
        max_shared = 5,
        ident.group = 'cell_type',
        target_metacells=50,
        min_cells=10
    )

    m_obj <- GetMetacellObject(test_seurat)

    expect_equal(all(c('cell_type', 'annotation') %in% colnames(m_obj@meta.data)), TRUE)

})

# test that the idents are expected
test_that("Metacell object Idents are expected", {
    data(test_seurat)

    test_seurat <- SetupForWGCNA(
        test_seurat,
        wgcna_name = 'test',
        features = rownames(test_seurat)
    )

    test_seurat <- MetacellsByGroups(
        test_seurat,
        group.by = c('Sample', 'cell_type', 'annotation'),
        k = 5, 
        max_shared = 5,
        ident.group = 'Sample',
        target_metacells=10,
        min_cells=10
    )

    m_obj <- GetMetacellObject(test_seurat)

    check <- all(as.character(Idents(m_obj)) %in% as.character(test_seurat$Sample))
    expect_equal(check, TRUE)

})

# test that normalizing the metacells works as intended
test_that("Normalizing the metacell object works",{
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
    m_obj <- GetMetacellObject(test_seurat)
    expr <- Seurat::GetAssayData(test_seurat, slot='data')

    check <- all.equal(ceiling(expr), expr) == TRUE
    expect_equal(check, FALSE)

})
