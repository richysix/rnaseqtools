#' Get count data for an RNA-seq dataset
#'
#' \code{get_counts} returns a data.frame of counts,
#' either raw or normalised (default)
#'
#' @param data df rna-seq data to get counts from
#' @param samples df a samples df to subset the data to
#' @param normalised logical indicating whether to return
#' normalised counts or raw
#'
#' @return data.frame
#' If there is no count data in the supplied data frame
#' (i.e.) no column names containing "count" then NULL is returned
#'
#' @examples
#' norm_counts <- get_counts(data)
#' counts <- get_counts(data, normalised = FALSE)
#' counts_subset <- get_counts(data, samples)
#'
#' @export
get_counts <- function(data, samples = NULL, normalised = FALSE) {
  if (normalised) {
    count_data <- get_cols(data, samples, ".normalised.counts?$")
  } else {
    count_data <- get_cols(data, samples, ".counts?$")
  }
  return(count_data)
}

## TO DO. Create get_tpm function
#' Get TPM data for an RNA-seq dataset
#'
#' \code{get_tpm} returns a data.frame of tpms,
#'
#' @param data df rna-seq data to get tpms from
#' @param samples df a samples df to subset the data to
#'
#' @return data.frame
#' If there is no tpm data in the supplied data frame
#' (i.e.) no column names containing "tpm" then NULL is returned
#'
#' @examples
#' tpms <- get_tpms(data)
#' tpms_subset <- get_tpms(data, samples)
#'
#' @export
get_tpms <- function(data, samples = NULL) {
  return(get_cols(data, samples, ".tpms?$"))
}

# alias for get_tpms
#' @export
get_tpm <- function(data, samples = NULL) {
  return(get_tpms(data, samples))
}


#' Get columns from an RNA-seq dataset
#'
#' \code{get_cols} is a helper function to get either counts,
#' normalised counts or tpm from an RNA-seq data frame
#'
#' @param data df rna-seq data to get columns from
#' @param samples df a samples df to subset the data to
#' @param cols Str to use as a search term
#'
#' @return data.frame
#' If there no matching columns in the supplied data frame
#' then NULL is returned
#'
#' @examples
#' counts <- get_cols(data, samples, ".counts?$")
#' norm_counts <- get_cols(data, samples, ".normalised.counts?$")
#' tpm <- get_cols(data, samples, ".tpms?$")
#'
get_cols <- function(data, samples = NULL, cols = ".counts?$") {
  selected_cols <- data[ , grepl(cols, names(data))]
  if (grepl("count", cols) & !grepl("normalised", cols)) {
    selected_cols <- selected_cols[ , !grepl("normalised", colnames(selected_cols))]
  }
  names(selected_cols) <- gsub(cols, "", names(selected_cols))
  if (ncol(selected_cols) == 0) {
    return(NULL)
  }

  # Subset and reorder count data
  if (!is.null(samples)) {
    check_samples_match_counts(selected_cols, samples)
    available_samples <- intersect(samples$sample, colnames(selected_cols))
    selected_cols <- dplyr::select(selected_cols, dplyr::all_of(available_samples))
  }

  return(selected_cols)
}

#' Get non count data for an RNA-seq dataset
#'
#' \code{get_gene_metadata} returns a data.frame of metadata,
#' associated with the genes of an RNA-seq dataset
#'
#' @param data df rna-seq data to get gene metadata from
#'
#' @return data.frame
#'
#' @examples
#' non_count_data <- get_gene_metadata(data)
#'
#' @export
get_gene_metadata <- function(data) {
  non_count_data <- data[ , !grepl("count", colnames(data)) & !grepl("tpm", colnames(data)) ]
  return(non_count_data)
}

#' Subset an RNA-seq dataset using a samples file
#'
#' \code{subset_to_samples} returns a data.frame of metadata,
#' counts and normalised counts subset to the supplied sample
#' data frame
#'
#' @param data df rna-seq data to get counts from
#' @param samples df a samples df to subset the data to
#' @param counts logical indicating whether to include raw counts in the results
#' @param normalised_counts logical indicating whether to include normalised counts in the results
#' @param tpm logical indicating whether to include tpm values in the results
#'
#' @return data.frame
#'
#' @examples
#' non_count_data <- subset_to_samples(data, samples)
#'
#' @export
subset_to_samples <-
  function(data, samples, counts = TRUE, normalised_counts = TRUE, tpm = FALSE) {
  subset_metadata <- get_gene_metadata(data)

  if (counts) {
    # check if the samples match the counts
    sample_counts <- tryCatch(
      get_counts(data, samples, normalised = FALSE),
      warning = function(w) {
        all_counts <- get_counts(data, normalised = FALSE)
        # subset to intersection of samples
        available_samples <- intersect(samples$sample, colnames(all_counts))
        rlang::warn(class = class(w),
                    message = paste(w$message, "Only samples in both were returned"))
        all_counts[ , available_samples ]
      }
    )
    if (is.null(sample_counts)) {
      rlang::warn(class = "no_counts",
                  "The rnaseq data.frame has no count columns")
    } else {
      # add back 'count' to end of column names
      colnames(sample_counts) <- paste(colnames(sample_counts), 'count')
    }
  } else {
    sample_counts = NULL
  }

  if (normalised_counts) {
    # subset to samples
    sample_norm_counts <- tryCatch(
      get_counts(data, samples, normalised = TRUE),
      warning = function(w) {
        all_norm_counts <- get_counts(data, normalised = TRUE)
        # subset to intersection of samples
        available_samples <- intersect(samples$sample, colnames(all_norm_counts))
        rlang::warn(class = class(w),
                    message = paste(w$message, "Only samples in both were returned"))
        all_norm_counts[ , available_samples ]
      }
    )
    if (is.null(sample_norm_counts)) {
      normalised_counts = FALSE
      rlang::warn(class = "no_norm_counts",
                  "The rnaseq data.frame has no normalised count columns")
    } else {
      # add back 'normalised count' to end of column names
      colnames(sample_norm_counts) <- paste(colnames(sample_norm_counts), 'normalised count')
    }
  } else {
    sample_norm_counts = NULL
  }

  if (tpm) {
    # subset to samples
    sample_tpm <- tryCatch(
      get_tpms(data, samples),
      warning = function(w) {
        all_tpm <- get_tpms(data)
        # subset to intersection of samples
        available_samples <- intersect(samples$sample, colnames(all_tpm))
        rlang::warn(class = class(w),
                    message = paste(w$message, "Only samples in both were returned"))
        all_tpm[ , available_samples ]
      }
    )
    if (is.null(sample_tpm)) {
      rlang::warn(class = "no_tpm",
                  "The rnaseq data.frame has no normalised count columns")
    } else {
      # add back 'tpm' to end of column names
      colnames(sample_tpm) <- paste(colnames(sample_tpm), 'tpm')
    }
  } else {
    sample_tpm = NULL
  }

  subset_data <- list(
    subset_metadata,
    sample_counts,
    sample_norm_counts,
    sample_tpm
  )
  return(purrr::list_cbind(subset_data))
}
