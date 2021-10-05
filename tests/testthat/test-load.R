context('load_data')
library(rnaseqtools)

load("test_data.rda")

# named vector. names are colnames that should be altered
# elements are what they should be changed to
colnames <- c('chr' = 'Chr', '#Chr' = 'Chr', 'start' = 'Start',
              'end' = 'End', 'strand' = 'Strand',
              'ID' = 'GeneID', 'e95 Ensembl Gene ID' = 'GeneID',
              'Gene ID' = 'GeneID', 'adjpval' = 'adjp', 'padj' = 'adjp',
              'Adjusted p value' = 'adjp', 'Gene name' = 'Name',
              'Gene' = 'GeneID' )
# create a test data frame with all the colnames that need changing
num_rows <- 5
test_df <- data.frame(
  matrix(seq_len(num_rows*length(colnames)), nrow = num_rows)
)
names(test_df) <- names(colnames)

test_that("standardise colnames",{
    expect_equal(names(standardise_colnames(test_df)),
                 unname(colnames))
})

expected_coltypes <- list(
  chr = readr::col_factor(),
  start = readr::col_integer(),
  end = readr::col_integer(),
  strand = readr::col_factor(),
  `e95 Ensembl Gene ID` = readr::col_character(),
  `Adjusted p value` = readr::col_double(),
  `Gene name` = readr::col_character(),
  `sample-1 count` = readr::col_integer(),
  `sample-2 count` = readr::col_integer(),
  `sample-3 count` = readr::col_integer(),
  `sample-4 count` = readr::col_integer(),
  `sample-5 count` = readr::col_integer(),
  `sample-6 count` = readr::col_integer(),
  `sample-1 normalised count` = readr::col_double(),
  `sample-2 normalised count` = readr::col_double(),
  `sample-3 normalised count` = readr::col_double(),
  `sample-4 normalised count` = readr::col_double(),
  `sample-5 normalised count` = readr::col_double(),
  `sample-6 normalised count` = readr::col_double()
)

test_that("coltypes",{
  expect_equal(set_col_types('test_data.tsv'), expected_coltypes)
})

test_that("colnames", {
  # colnames are currently not the same when tmp data is loaded back in
  # because colnames get changed
  expect_true(any(names(test_all_data) != names(load_rnaseq_data('test_data.tsv'))))
})
# correct test_all_data names
names(test_all_data) <- c('Chr', 'Start', 'End', 'Strand', 'GeneID', 'adjp',
                          'Name', paste0('sample-', seq_len(6), ' count'),
                          paste0('sample-', seq_len(6), ' normalised count'))

test_that("load RNAseq data",{
  expect_equal(as.data.frame(load_rnaseq_data('test_data.tsv')),
               as.data.frame(test_all_data), check.attributes = FALSE)
})

# test loading samples file
samples <- load_rnaseq_samples('test_samples.tsv')
test_that("load samples", {
  expect_equal(samples$sample, samples_data$sample)
  expect_equal(samples$condition, samples_data$condition)
  expect_equal(samples$sex, samples_data$sex)
})

## DeTCT data
detct_colnames <- c('chr' = 'Chr', '#Chr' = 'Chr', 'Region start' = 'RegionStart',
              'Region end' = 'RegionEnd', "3' end position" = '3PrimeEndPosition',
              "3' end strand" = '3PrimeEndStrand', "3' end read count" = '3PrimeEndReadCount',
              'p value' = 'pval', 'Adjusted p value' = 'adjp',
              "Distance to 3' end" = 'DistanceTo3PrimeEnd',
              'e98 Ensembl Gene ID' = 'GeneID',
              "Gene type" = 'GeneType',
              'e98 Ensembl Transcript ID' = 'TranscriptID',
              "Transcript type" = 'TranscriptType',
              'Gene name' = 'Name', 'Gene description' = 'Description')

# create a test data frame with all the colnames that need changing
num_rows <- 5
test_df <- data.frame(
  matrix(seq_len(num_rows*length(detct_colnames)), nrow = num_rows)
)
names(test_df) <- names(detct_colnames)

test_that("standardise colnames",{
  expect_equal(names(standardise_colnames(test_df)),
               unname(detct_colnames))
})

# check loading DeTCT data
detct_data <- load_rnaseq_data('test_detct_data.tsv')
test_that("load DeTCT data", {
  expect_equal(detct_data, test_detct_data, check.attributes = FALSE)
})

test_that("check col names", {
  expect_equal(colnames(detct_data),
               colnames(standardise_colnames(test_detct_data)),
               check.attributes = FALSE)
})

test_that("check col types", {
  expect_equal(purrr::map_chr(detct_data, class),
               purrr::map_chr(test_detct_data, class), check.attributes = FALSE)
})

teardown({
  unlink('test_data.tmp')
})
