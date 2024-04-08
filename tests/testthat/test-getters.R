library(testthat)
library(SingleCellAlleleExperiment)
library(scaeData)

example_data_5k <- scaeData::scaeDataGet(dataset="pbmc_5k")
lookup <- utils::read.csv(system.file("extdata",
                                      "pbmc_5k_lookup_table.csv",
                                      package="scaeData"))

scae <- read_allele_counts(example_data_5k$dir,
                           sample_names="example_data_wta",
                           filter_mode="yes",
                           lookup_file=lookup,
                           barcode_file=example_data_5k$barcodes,
                           gene_file=example_data_5k$features,
                           matrix_file=example_data_5k$matrix,
                           verbose=TRUE)

rn_c_scae <- rownames(counts(scae))

scae_ni_genes <- get_nigenes(scae)
scae_alleles <- scae_subset_alleles(scae)
scae_a_genes <- get_agenes(scae)
scae_functional <- scae_subset_functional(scae)

# non_immune genes layer
test_that("non-immune genes getter", {
  expect_equal(scae_ni_genes, scae[grepl("ENSG", rn_c_scae, fixed=TRUE),])
})

# alleles layer
test_that("alleles getter", {
  expect_equal(scae_alleles, scae[grepl("*", rn_c_scae, fixed=TRUE),])
})

# immune gene layer
test_that("immune genes getter", {
  expect_equal(scae_a_genes, scae[grepl("HLA-", rn_c_scae, fixed=TRUE),])
})

# functional class layer
test_that("functional class getter", {
  expect_equal(scae_functional, scae[grepl("class", rn_c_scae, fixed=TRUE),])
})

test_that("test wrapper getter", {
  #nonimmune
  expect_equal(scae_subset(scae, "nonimmune"), scae_ni_genes)
  #allele
  expect_equal(scae_subset(scae, "alleles"), scae_alleles)
  #immune
  expect_equal(scae_subset(scae, "immune_genes"), scae_a_genes)
  #functional
  expect_equal(scae_subset(scae, "functional_groups"), scae_functional)

  msg <- "Invalid layer specified, Choose from `nonimmune`, `alleles`, `immune_genes`, `functional_groups`"
  expect_message(scae_subset(scae, "wrong_layer"), regexp=msg)

})
