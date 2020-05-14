context('load_data')
library(rnaseqtools)

# setup test data file to load
num_rows <- 100
set.seed(208)
starts <- sample(1:10000, num_rows)
test_all_data <- data.frame(
  'chr' = sample(1:25, num_rows, replace = TRUE),
  'start' = starts,
  'end' = as.integer(starts + 100),
  'strand' = sample(c('1', '-1'), num_rows, replace = TRUE),
  'e95 Ensembl Gene ID' = paste0('ZFG', seq_len(num_rows)),
  'Adjusted p value' = runif(num_rows),
  'Gene name' = paste0('gene-', seq_len(num_rows))
)
test_all_data$chr <- factor(test_all_data$chr)
sample_num <- 6
counts_list <- vector("list", length = sample_num)
for(num in seq_len(sample_num)) {
  sample_name <- paste0('sample-', num, ' count')
  counts <- as.integer(floor(runif(num_rows)*100))
  test_all_data[[sample_name]] <- counts
  counts_list[[num]] <- counts
}
counts <- as.data.frame(do.call('cbind', counts_list))
names(counts) <- paste0('sample-', seq_len(sample_num))

norm_counts_list <- vector("list", length = sample_num)
for(num in seq_len(sample_num)) {
  sample_name <- paste0('sample-', num, ' normalised count')
  norm_counts <- runif(num_rows)*100
  test_all_data[[sample_name]] <- norm_counts
  norm_counts_list[[num]] <- norm_counts
}
norm_counts <- as.data.frame(do.call('cbind', norm_counts_list))
names(norm_counts) <- paste0('sample-', seq_len(sample_num))

names(test_all_data) <- gsub("\\.", " ", names(test_all_data))
write.table(test_all_data, file = "test_data.tmp", quote = FALSE,
            row.names = FALSE, col.names = TRUE, sep = "\t")

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

test_df <- data.frame(
  'Chr' = as.integer(sample(seq_len(num_rows), num_rows)),
  'Start' = as.numeric(sample(1:10000, num_rows)),
  'End' = as.numeric(sample(1:10000, num_rows)),
  'Strand' = as.integer(sample(c('1', '-1'), num_rows, replace = TRUE))
)
test_df <- standardise_coltypes(test_df)
test_that("coltypes",{
  expect_equal(class(test_df$Chr), 'factor')
  expect_equal(class(test_df$Start), 'integer')
  expect_equal(class(test_df$End), 'integer')
  expect_equal(class(test_df$Strand), 'factor')
})

test_that("load RNAseq data",{
  # colnames are currently not the same when tmp data is loaded back in
  # because colnames get changed
  expect_true(any(names(test_all_data) != names(load_rnaseq_data('test_data.tmp'))))
  # correct test_all_data names
  names(test_all_data) <- c('Chr', 'Start', 'End', 'Strand', 'GeneID', 'adjp',
                            'Name', paste0('sample-', seq_len(6), ' count'),
                            paste0('sample-', seq_len(6), ' normalised count'))

  expect_equal(load_rnaseq_data('test_data.tmp'),
               test_all_data)
})

test_that("get counts", {
  expect_equal(get_counts(test_all_data), counts)
  expect_equal(get_counts(test_all_data, normalised = TRUE), norm_counts)
})

teardown({
  unlink('test_data.tmp')
})
