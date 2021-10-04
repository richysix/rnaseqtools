library(tibble)
library(readr)
library(rprojroot)

root_path <- find_root(is_rstudio_project)

# setup test data file to load
num_rows <- 100
set.seed(208)
starts <- sample(1:10000, num_rows)
test_all_data <- tibble(
  'chr' = sample(1:25, num_rows, replace = TRUE),
  'start' = starts,
  'end' = as.integer(starts + 100),
  'strand' = sample(c('1', '-1'), num_rows, replace = TRUE),
  'e95 Ensembl Gene ID' = paste0('ZFG', seq_len(num_rows)),
  'Adjusted p value' = runif(num_rows),
  'Gene name' = paste0('gene-', seq_len(num_rows))
)
test_all_data$chr <- factor(test_all_data$chr, levels = unique(test_all_data$chr))
test_all_data$strand <- factor(test_all_data$strand)
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
counts <- as_tibble(counts)

norm_counts_list <- vector("list", length = sample_num)
for(num in seq_len(sample_num)) {
  sample_name <- paste0('sample-', num, ' normalised count')
  norm_counts <- runif(num_rows)*100
  test_all_data[[sample_name]] <- norm_counts
  norm_counts_list[[num]] <- norm_counts
}
norm_counts <- as.data.frame(do.call('cbind', norm_counts_list))
names(norm_counts) <- paste0('sample-', seq_len(sample_num))
norm_counts <- as_tibble(norm_counts)

# write a file for testing load_rnaseq_data
write.table(test_all_data, quote = FALSE,
            row.names = FALSE, col.names = TRUE, sep = "\t",
            file = file.path(root_path, "tests", "testthat", "test_data.tsv") )

# create a samples test file
samples_data <- tibble(
  sample = factor(paste0('sample-', seq_len(sample_num))),
  condition = factor(rep(c('wt', 'mut'), each = 3),
                     levels = c('wt', 'mut')),
  sex = rep(c('M', 'F'), 3)
)
# write a file for testing load_rnaseq_samples
write_tsv(samples_data, file = file.path(root_path, "tests", "testthat", "test_samples.tsv") )

# save test data as R objects
save(test_all_data, counts, norm_counts, samples_data,
     file = file.path(root_path, "tests", "testthat", "test_data.rda"))
