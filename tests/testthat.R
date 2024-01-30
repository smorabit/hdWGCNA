library(testthat)
library(hdWGCNA)

#test_check("hdWGCNA")
devtools::test()

# remove files
file.remove('tests/testthat/ConsensusTOM-block.1.rda')