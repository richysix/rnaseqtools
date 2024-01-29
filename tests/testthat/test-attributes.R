context('attributes')
library(rnaseqtools)
library(dplyr)

load("test_data.rda")

# this to test
# get_counts: get counts, norm counts, test missing samples
# get counts from test_all_data
test_counts_only <- get_counts(data = test_all_data, samples = samples_data)
test_that("get_counts",{
  expect_equal(colnames(test_counts_only),
               colnames(counts))
  expect_equal(dim(test_counts_only),
               dim(counts))
  expect_identical(test_counts_only, counts)
  expect_equal(get_counts(test_all_data), counts)
  expect_equal(get_counts(test_all_data, normalised = TRUE), norm_counts)
})

# create new samples objects for testing
sample_subset <- filter(samples_data, sex == "F")
sample_missing_from_counts <- rbind(samples_data,
                                    tibble(
                                      sample = "sample-7",
                                      condition = "het",
                                      sex = "M"
                                    ))

test_that("get_counts_errors_and_warnings",{
  expect_warning(get_counts(test_all_data, samples = sample_subset),
                 class = "missing_from_samples")
  expect_warning(get_counts(test_all_data, samples = sample_missing_from_counts),
                 class = "missing_from_counts")
})

# get_gene_metadata: check columns
metadata <- get_gene_metadata(test_all_data)
test_that("get_gene_metadata",{
  expect_equal(colnames(metadata), colnames(test_all_data)[ !grepl("count", colnames(test_all_data)) ])
  expect_equal(dim(metadata), dim(test_all_data[ , !grepl("count", colnames(test_all_data)) ]))
  expect_equal(metadata, test_all_data[ , !grepl("count", colnames(test_all_data)) ])
})

# subset_to_samples
test_that("subset_to_samples", {
  expect_warning(subset_to_samples(test_all_data, sample_subset),
                 class = "missing_from_samples")
  expect_warning(subset_to_samples(test_all_data, sample_missing_from_counts),
                 class = "missing_from_counts")
  expect_equal(dim(suppressWarnings(subset_to_samples(test_all_data, sample_subset))), c(100,13))
  expect_equal(suppressWarnings(subset_to_samples(test_all_data, sample_subset)),
               test_all_data[ , c(1:7,9,11,13,15,17,19)])
  test_all_data_no_counts <- test_all_data[ , !grepl("count", colnames(test_all_data)) |
                                              grepl("normalised count", colnames(test_all_data)) ]
  expect_warning(subset_to_samples(test_all_data_no_counts, samples_data),
                 class = "no_counts")
  test_all_data_no_norm_counts <- test_all_data[ , !grepl("normalised count", colnames(test_all_data)) ]
  expect_warning(subset_to_samples(test_all_data_no_norm_counts, samples_data),
                 class = "no_norm_counts")
})
