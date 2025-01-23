library(tibble)
library(magrittr)
library(readr)
library(rprojroot)
library(DESeq2)

root_path <- find_root(is_rstudio_project)

# create a samples test file
sample_num <- 6
samples_data <- tibble(
  sample = factor(paste0('sample-', seq_len(sample_num))),
  condition = factor(rep(c('wt', 'mut'), each = 3),
                     levels = c('wt', 'mut')),
  sex = rep(c('M', 'F'), 3)
)
# write a file for testing load_rnaseq_samples
write_tsv(samples_data, file = file.path(root_path, "tests", "testthat", "test_samples.tsv") )

# setup test data file to load
num_rows <- 100
set.seed(208)
starts <- sample(1:10000, num_rows)
test_all_data <- tibble(
  'chr' = sample(1:25, num_rows, replace = TRUE),
  'start' = starts,
  'end' = as.integer(starts + sample(100:1000, num_rows)),
  'strand' = sample(c('1', '-1'), num_rows, replace = TRUE),
  'e95 Ensembl Gene ID' = paste0('ZFG', seq_len(num_rows)),
  'Adjusted p value' = runif(num_rows),
  'Gene name' = paste0('gene-', seq_len(num_rows))
)
test_all_data$chr <- factor(test_all_data$chr, levels = unique(test_all_data$chr))
test_all_data$strand <- factor(test_all_data$strand)
counts_list <- vector("list", length = sample_num)
for (num in seq_len(sample_num)) {
  sample_name <- paste0('sample-', num, ' count')
  counts <- as.integer(floor(runif(num_rows)*100))
  test_all_data[[sample_name]] <- counts
  counts_list[[num]] <- counts
}
counts <- as.data.frame(do.call('cbind', counts_list))
names(counts) <- paste0('sample-', seq_len(sample_num))
counts <- as_tibble(counts)

DESeqData <- DESeq2::DESeqDataSetFromMatrix(counts, samples_data, design = ~ condition)
DESeqData <- DESeq2::estimateSizeFactors(DESeqData)
rowData(DESeqData) <- dplyr::select(test_all_data, chr:`Gene name`)
# get normalised counts
norm_counts <- DESeq2::counts(DESeqData, normalized = TRUE) |>
  as_tibble()
# add back "normalised" to colnames
colnames(norm_counts) <- sub("$", " normalised count", colnames(norm_counts))
# add count back to unnormalised counts
colnames(counts) <- sub("$", " count", colnames(counts))

test_all_data <- tibble(
  test_all_data,
  norm_counts
)

# write a file for testing load_rnaseq_data
write.table(test_all_data, quote = FALSE,
            row.names = FALSE, col.names = TRUE, sep = "\t",
            file = file.path(root_path, "tests", "testthat", "test_data.tsv") )

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

# make TPM
# add tx lengths to deseq object
TxLengths <- (test_all_data$end - test_all_data$start)
assays(DESeqData, withDimnames = FALSE)[["avgTxLength"]] <-
  matrix(rep(TxLengths, ncol(DESeqData)), ncol = ncol(DESeqData))

rnaseq_fpkm <- fpkm(DESeqData)
rownames(rnaseq_fpkm) <- rowData(DESeqData)$`e95 Ensembl Gene ID`

# calculate tpm from fpkm
fpkm_to_tpm <- function(fpkm) {
  total_fpkm_by_sample <- colSums(fpkm)
  do.call(cbind, lapply(1:ncol(fpkm), function(i) {
    exp( log(fpkm[,i]) - log(total_fpkm_by_sample[i]) + log(1e6) )
  }))
}
tpm <- fpkm_to_tpm(rnaseq_fpkm) |>
  as_tibble(.name_repair = "minimal") |>
  set_colnames(colnames(rnaseq_fpkm))
tpm_data <- tibble(
  dplyr::select(test_all_data, chr:`Gene name`),
  tpm |> set_colnames(sub("$", " tpm", colnames(tpm)))
)
write_tsv(tpm_data, file = file.path(root_path, "tests", "testthat", "test_tpm_data.tsv"))

# save test data as R objects
save(
  test_all_data,
  test_detct_data,
  counts,
  norm_counts,
  samples_data,
  tpm,
  tpm_data,
  file = file.path(root_path, "tests", "testthat", "test_data.rda")
)
