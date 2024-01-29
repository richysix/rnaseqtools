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
    count_data <- data[,grepl(".normalised.counts?$", names(data))]
    names(count_data) <- gsub(".normalised.counts?$", "", names(count_data))
  } else {
    count_data <- data[,grepl(".counts?$", names(data)) &
                         !grepl("normalised", names(data))]
    names(count_data) <- gsub(".counts?$", "", names(count_data))
  }
  if (ncol(count_data) == 0) {
    return(NULL)
  }

  # Subset and reorder count data
  if (!is.null(samples)) {
    check_samples_match_counts(count_data, samples)
    available_samples <- intersect(samples$sample, colnames(count_data))
    count_data <- dplyr::select(count_data, dplyr::one_of(available_samples))
  }

  return(count_data)
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
  non_count_data <- data[ , !grepl("count", colnames(data)) ]
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
#'
#' @return data.frame
#'
#' @examples
#' non_count_data <- subset_to_samples(data, samples)
#'
#' @export
subset_to_samples <- function(data, samples, counts = TRUE, normalised_counts = TRUE) {
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
    # add back 'count' to end of column names
    colnames(sample_counts) <- paste(colnames(sample_counts), 'count')
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
    # add back 'normalised count' to end of column names
    colnames(sample_norm_counts) <- paste(colnames(sample_norm_counts), 'normalised count')
  }

  if (counts) {
    if (normalised_counts) {
      subset_data <- tibble::as_tibble(cbind(subset_metadata, sample_counts, sample_norm_counts))
    } else {
      subset_data <- tibble::as_tibble(cbind(subset_metadata, sample_counts))
    }
  } else {
    if (normalised_counts) {
      subset_data <- tibble::as_tibble(cbind(subset_metadata, sample_norm_counts))
    } else {
      subset_data <- tibble::as_tibble(subset_metadata)
    }
  }
  return(subset_data)
}
