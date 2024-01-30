
# check that the soft power table has the expected values
test_that("Soft power results has expected power values", {

    data(test_seurat)

    selected_powers <- c(4,6,8)
    test_seurat <- TestSoftPowers(test_seurat, powers=selected_powers, corFnc='cor')

    pt <- GetPowerTable(test_seurat)
    check <- all(selected_powers %in% pt$Power )

    expect_equal(check, TRUE)
})

