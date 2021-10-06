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
names(counts) <- paste0('sample-', seq_len(sample_num), ' count')
counts <- as_tibble(counts)

norm_counts_list <- vector("list", length = sample_num)
for(num in seq_len(sample_num)) {
  sample_name <- paste0('sample-', num, ' normalised count')
  norm_counts <- runif(num_rows)*100
  test_all_data[[sample_name]] <- norm_counts
  norm_counts_list[[num]] <- norm_counts
}
norm_counts <- as.data.frame(do.call('cbind', norm_counts_list))
names(norm_counts) <- paste0('sample-', seq_len(sample_num), ' normalised count')
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

samples_txt_data <- data.frame(
  row.names = factor(paste0('sample-', seq_len(sample_num))),
  condition = factor(rep(c('wt', 'mut'), each = 3),
                     levels = c('wt', 'mut')),
  sex = rep(c('M', 'F'), 3)
)
# write a file for testing load_rnaseq_samples
write.table(samples_txt_data, quote = FALSE,
            row.names = TRUE, col.names = NA, sep = "\t",
            file = file.path(root_path, "tests", "testthat", "test_samples.txt") )

# create test DeTCT data object
num_rows <- 100
set.seed(802)
starts <- sample(1:10000, num_rows)
test_metadata <- tibble(
  'Chr' = sample(1:25, num_rows, replace = TRUE),
  'Region start' = starts,
  'Region end' = as.integer(starts + 100),
  "3' end position" = as.character(starts + 200),
  "3' end strand" = sample(c('1', '-1'), num_rows, replace = TRUE),
  "3' end read count" = as.character(sample(1:1000, num_rows)),
  "p value" = runif(num_rows),
  'Adjusted p value' = runif(num_rows),
  "Distance to 3' end" = as.character(sample(1:1000, num_rows)),
  'e98 Ensembl Gene ID' = paste0('ZFG', seq_len(num_rows)),
  'Gene type' = "protein_coding",
  "e98 Ensembl Transcript ID" = paste0('ZFT', seq_len(num_rows)),
  "Transcript type" = "protein_coding:appris::",
  'Gene name' = paste0('gene-', seq_len(num_rows)),
  'Gene description' = paste0('gene-', seq_len(num_rows), ' extra stuff')
)
test_detct_data <- as_tibble(cbind(test_metadata, counts, norm_counts))
test_detct_data$Chr <- factor(test_detct_data$Chr, levels = unique(test_detct_data$Chr))
test_detct_data$`3' end strand` <- factor(test_detct_data$`3' end strand`)

write_tsv(test_detct_data, file = file.path(root_path, "tests", "testthat", "test_detct_data.tsv") )

# remove suffixes from counts and norm_counts
colnames(counts) <- sub('.count', "", colnames(counts))
colnames(norm_counts) <- sub('.normalised.count', "", colnames(norm_counts))

# save test data as R objects
save(test_all_data, test_detct_data, counts, norm_counts, samples_data,
     file = file.path(root_path, "tests", "testthat", "test_data.rda"))
