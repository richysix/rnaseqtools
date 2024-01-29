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

# test loading a samples file where the first column isn't labelled
samples_txt <- load_rnaseq_samples('test_samples.txt')
test_that("load samples.txt", {
  expect_equal(colnames(samples_txt), c('sample', 'condition', 'sex'))
  expect_equal(class(samples_txt$condition), 'factor')
  expect_equal(class(samples_txt$sex), 'character')
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

# check loading DeTCT data with load_detct data
# this function adds a RegionID column
detct_data <- load_detct_data('test_detct_data.tsv')
test_that("load DeTCT data", {
  expect_equal(colnames(detct_data)[1], "RegionID")
})

# check check_samples_match_counts
test_that("check_samples_match_counts works", {
  expect_true(check_samples_match_counts(counts, samples_data))
  # remove one sample from samples
  samples_tmp <- samples_data[ samples_data$sample != "sample-1", ]
  expect_warning(check_samples_match_counts(counts, samples_tmp),
                 class = "missing_from_samples")
  # remove a sample from counts
  counts_tmp <- counts[ , colnames(counts) != "sample-6"]
  expect_warning(check_samples_match_counts(counts_tmp, samples_data),
                 class = "missing_from_counts")
  expect_warning(check_samples_match_counts(counts_tmp, samples_tmp),
                 class = "missing_from_both_samples_and_counts")
})

test_that("check_samples works", {
  expect_true(check_samples(test_all_data, samples_data))
  # check sample missing from sample data
  # remove one sample from samples
  samples_tmp <- samples_data[ samples_data$sample != "sample-1", ]
  expect_warning(check_samples(test_all_data, samples_tmp), "One or more extra samples in the counts data")
  # check sample missing from counts
  all_data_tmp <- test_all_data[ , colnames(test_all_data) != "sample-6 count"]
  expect_warning(check_samples_match_counts(all_data_tmp, samples_data), "One or more samples are missing from the counts data")
  all_data_tmp <- test_all_data[ , colnames(test_all_data) != "sample-6 normalised count"]
  expect_warning(check_samples_match_counts(all_data_tmp, samples_data), "One or more samples are missing from the counts data")
})

test_that("normalise_counts works", {
  test_data_counts <- test_all_data[ , !grepl("normalised", colnames(test_all_data)) ]
  expect_equal(normalise_counts(test_data_counts, samples_data), test_all_data)
})

teardown({
  unlink('test_data.tmp')
})
