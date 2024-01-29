#' Load RNA-seq count data in from file
#'
#' \code{load_rnaseq_data} opens a file containing count data
#' and returns a tibble with column names adjusted to
#' take account of different column name formats.
#'
#' @param data_file char, Name of the file to open
#' @param ... Arguments passed to \code{\link[readr]{read_tsv}}
#'
#' @return tibble
#'
#' @examples
#' data <- load_rnaseq_data(data_file = 'all.tsv')
#'
#' @export
load_rnaseq_data <- function(data_file, ...) {
  # set col types
  coltypes <- set_col_types(data_file)

  # Read data
  data <- readr::read_tsv(data_file, col_types = do.call(readr::cols, coltypes), ...)
  data <- standardise_colnames(data)

  return(data)
}

#' Load DeTCT count data in from file
#'
#' \code{load_detct_data} opens a file containing count data
#' and returns a tibble with column names adjusted to
#' take account of different column name formats.
#' If it doesn't already exist, it creates a column called RegionID
#' made up of Chr:RegionStart:RegionEnd:3PrimeEndPosition:3PrimeEndStrand:GeneID
#'
#' @param data_file char, Name of the file to open
#' @param ... Arguments passed to \code{\link[readr]{read_tsv}}
#'
#' @return tibble
#'
#' @examples
#' data <- load_detct_data(data_file = 'all.tsv')
#'
#' @export
load_detct_data <- function(data_file, ...) {
  # set col types
  coltypes <- set_col_types(data_file)

  # Read data
  data <- readr::read_tsv(data_file, col_types = do.call(readr::cols, coltypes), ...)
  data <- standardise_colnames(data)

  if(!("RegionID" %in% colnames(data))) {
    data <- data %>%
      dplyr::mutate(., RegionID = paste(Chr, RegionStart, RegionEnd, `3PrimeEndPosition`,
                                        `3PrimeEndStrand`, GeneID, sep = ":")) %>%
      dplyr::relocate(., RegionID)
  }

  return(data)
}

#' Set the column types based on the column names
#'
#' \code{set_col_types} loads the first line of
#' the data and sets the column type based on the column names.
#' Chr and Strand are factor and Start and End are integer.
#' normalised count columns are double (float) and
#' count columns are integer.
#' All other columns are set to the type guessed by read_tsv
#'
#' @param data_file Name of file to open
#'
#' @return data.frame
#'
#' @examples
#' data <- standardise_colnames('test_data.tsv')
#'
set_col_types <- function(data_file){
  types_for_cols = c(
    "Chr" = "f",
    "Start" = "i",
    "End" = "i",
    "Strand" = "f",
    "RegionStart" = "i",
    "RegionEnd" = "i",
    "RegionID" = "c",
    "3PrimeEndPosition" = "c",
    "3PrimeEndStrand" = "f",
    "3PrimeEndReadCount" = "c",
    "DistanceTo3PrimeEnd" = "c",
    "pval" = "d",
    "adjp" = "d"
  )

  # get columns
  header <- readr::read_tsv(data_file, n_max = 1000,
                            col_types = readr::cols())
  standard_colnames <- colnames(standardise_colnames(header))

  coltype_for_column_name <- function(i, standard_colnames, types_for_cols, header) {
    standard_colname <- standard_colnames[i]
    original_colname <- colnames(header)[i]
    # check if the colname is in types_for_cols
    if (standard_colname %in% names(types_for_cols)) {
      if (types_for_cols[[standard_colname]] == "f") {
        return(readr::col_factor())
      } else if (types_for_cols[[standard_colname]] == "i") {
        return(readr::col_integer())
      } else if (types_for_cols[[standard_colname]] == "c") {
        return(readr::col_character())
      } else if (types_for_cols[[standard_colname]] == "d") {
        return(readr::col_double())
      }
    } else if (grepl("count$", standard_colname)) {
      if (grepl("normalised", standard_colname)) {
        return(readr::col_double())
      } else {
        # count columns should be integer
        return(readr::col_integer())
      }
    } else {
      # everything else should be as parsed
      return(readr::spec(header)$cols[[original_colname]])
    }
  }
  # set column types
  column_types <- purrr::map(seq_len(ncol(header)), coltype_for_column_name,
                             standard_colnames, types_for_cols, header)
  names(column_types) <- colnames(header)
  return(column_types)
}

#' Standardise count file column names
#'
#' \code{standardise_colnames} takes a dataframe
#' loaded from a RNA-seq all file and creates standard
#' column names.
#'
#' @param data data.frame, data to standardise
#'
#' @return data.frame
#'
#' @examples
#' data <- standardise_colnames(data)
#'
#' @export
standardise_colnames <- function(data) {
  # Support different column names
  names(data)[names(data) == 'chr']               <- 'Chr'
  names(data)[names(data) == '#Chr']              <- 'Chr'
  names(data)[names(data) == 'start']             <- 'Start'
  names(data)[names(data) == 'end']               <- 'End'
  names(data)[names(data) == 'strand']            <- 'Strand'
  names(data)[names(data) == 'ID']                <- 'GeneID'
  names(data)[ grepl("e[0-9]+ Ensembl Gene ID",
                     names(data)) ]               <- 'GeneID'
  names(data)[names(data) == 'Gene ID']           <- 'GeneID'
  names(data)[names(data) == 'Gene']              <- 'GeneID'
  names(data)[names(data) == 'adjpval']           <- 'adjp'
  names(data)[names(data) == 'padj']              <- 'adjp'
  names(data)[names(data) == 'Adjusted p value']  <- 'adjp'
  names(data)[names(data) == 'Gene name']         <- 'Name'

  # detct names
  names(data)[names(data) == 'p value']                 <- 'pval'
  names(data)[names(data) == 'Gene description']        <- 'Description'
  names(data)[names(data) == 'Region start']            <- 'RegionStart'
  names(data)[names(data) == 'Region end']              <- 'RegionEnd'
  names(data)[names(data) == "3' end position"]         <- '3PrimeEndPosition'
  names(data)[names(data) == "3' end strand"]           <- '3PrimeEndStrand'
  names(data)[names(data) == "3' end read count"]       <- '3PrimeEndReadCount'
  names(data)[names(data) == "Distance to 3' end"]      <- 'DistanceTo3PrimeEnd'
  names(data)[names(data) == "Gene type"]               <- 'GeneType'
  names(data)[ grepl("e[0-9]+ Ensembl Transcript ID",
                     names(data)) ]                     <- 'TranscriptID'
  names(data)[names(data) == "Transcript type"]         <- 'TranscriptType'

  return(data)
}

#' Load RNA-seq samples file
#'
#' \code{load_rnaseq_samples} opens a file containing sample
#' data and returns a tibble
#' It expects the file to contain columns named sample and
#' condition which it will make factors.
#'
#' @param samples_file Name of the file to open
#'
#' @return tibble
#'
#' @examples
#' data <- load_rnaseq_samples(samples_file = 'test_samples.tsv')
#'
#' @export
load_rnaseq_samples <- function(samples_file) {
  # check header
  samples <-
    tryCatch(readr::read_tsv(samples_file,
                             col_types = readr::cols(
                               sample = readr::col_factor(),
                               condition = readr::col_factor()
                             )
    ),
    warning = function(w){
      # print(w)
      # print(w$message)
      if(grepl("The following named parsers don't match the column names", w$message)) {
        # print("X1 match")
        header <- suppressWarnings(readr::read_tsv(samples_file, n_max = 2,
                                                   col_types = readr::cols()))
        cols_list <- readr::spec(header)$cols
        names(cols_list)[ names(cols_list) == "...1" ] <- "sample"
        cols_list[['sample']] <- readr::col_factor()
        if('condition' %in% names(cols_list)) {
          cols_list[['condition']] <- readr::col_factor()
        }
        return(readr::read_tsv(samples_file, skip = 1,
                               col_names = names(cols_list),
                               col_types = do.call(readr::cols, cols_list)))
      } else {
        rlang::warn(message = w$message, class = "sample_load")
      }
    },
    error = function(e){ stop(e) }
    # message = function(m){ message(m) }
    )

  return(samples)
}

#' Adjust RNAseq count data for library size
#'
#' `normalise_counts()` takes RNAseq count data and produces normalised counts
#' that are adjusted for library size (sequencing depth). The normalised count
#' columns are added to the supplied data frame with column names that are
#' the original column names with " normalised count" appended.
#'
#' @param rnaseq_data data frame of rnaseq data. Must contain columns labelled
#' "{sample} count".
#' @param samples data frame of sample metadata. Must included a `condition` column.
#'
#' @return tibble containing the original data and the normalised count columns added
#' @export
#'
#' @examples
#'
#' rnaseq_data <- test_all_data[ , !grepl("normalised", colnames(test_all_data)) ]
#' rnaseq_data <- normalise_counts(rnaseq_data, samples_data)
#'
normalise_counts <- function(rnaseq_data, samples){
  # get counts
  counts <- get_counts(rnaseq_data)
  # check samples match
  check_samples_match_counts(counts, samples)
  # make DESeq2 object
  DESeqData <- DESeq2::DESeqDataSetFromMatrix(counts, samples, design = ~ condition)
  DESeqData <- DESeq2::estimateSizeFactors(DESeqData)
  # get normalised counts
  norm_counts <- DESeq2::counts(DESeqData, normalized = TRUE)
  # add back "normalised" to colnames
  colnames(norm_counts) <- sub("$", " normalised count", colnames(norm_counts))
  # add to rnaseq data
  return(tibble::tibble(rnaseq_data, tibble::as_tibble(norm_counts)))
}


#' See if sample and count data.frames have matching samples
#'
#' @description
#' `check_samples()` takes sample and data data.frames and checks that
#'  the samples match between the two. It checks both `counts` and
#'  `normalised counts`. If there are samples
#'  in the samples data.frame that don't exist in the corresponding counts
#'  data.frame, an error is raised. If there are samples in the counts data.frame that
#'  aren't in the samples data.frame a warning is raised. If the samples
#'  match exactly, the function returns TRUE invisibly.
#'
#'  `check_samples_match_counts()` takes sample and count data.frames
#'  and checks that the samples match between the two. This does the actual checking
#'  that the samples match. `check_samples()` gets the counts and normalised counts
#'  form the data data.frame and passes them on to `check_samples_match_counts()`
#'
#' @param rnaseq_data data.frame of gene metadata, counts and normalised counts
#' @param count_data data.frame of counts only
#' @param samples data.frame of sample info
#'
#' @return TRUE invisibly if the samples match
#'
#' @export
#'
#' @examples
#'
#' check_samples_match_counts(test_all_data, samples_data)
#'
#' check_samples_match_counts(counts, samples_data)
#'
#' check_samples_match_counts(norm_counts, samples_data)
#'
check_samples <- function(rnaseq_data, samples) {
  # get counts
  # don't use the samples df
  counts <- get_counts(rnaseq_data)
  check_samples_match_counts(counts, samples)
  # same for normalised counts
  norm_counts <- get_counts(rnaseq_data, normalised = TRUE)
  check_samples_match_counts(norm_counts, samples)
  return(invisible(TRUE))
}

#' @rdname check_samples
#' @export
check_samples_match_counts <- function(count_data, samples) {
  samples_in_counts <- colnames(count_data)
  extra_samples <- setdiff(samples_in_counts, samples$sample)
  missing_samples <- setdiff(samples$sample, samples_in_counts)

  if(length(extra_samples) > 0) {
    if(length(missing_samples) > 0) {
      rlang::warn(class = "missing_from_both_samples_and_counts",
                  message = paste(
                    paste("One or more samples are missing from the counts data,",
                          paste0(missing_samples, collapse = ", ")),
                    paste("One or more extra samples in the counts data,",
                          paste0(extra_samples, collapse = ", ")),
                    sep = "\n")
      )
    } else {
      rlang::warn(class = "missing_from_samples",
                  message = paste("One or more extra samples in the counts data,",
                                  paste0(extra_samples, collapse = ", "))
      )
    }
  } else if(length(missing_samples) > 0) {
    rlang::warn(class = "missing_from_counts",
                message = paste("One or more samples are missing from the counts data,",
                                paste0(missing_samples, collapse = ", "))
    )
  }

  return(invisible(TRUE))
}
