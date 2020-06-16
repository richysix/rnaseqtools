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

  # Subset and reorder count data
  # TO ADD: check all samples exist in counts
  # and if any are in counts that aren't in samples
  if (!is.null(samples)) {
    count_data <- select(count_data, one_of(as.character(samples$sample)))
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
#' associated with the genes of an RNA-seq dataset
#'
#' @param data df rna-seq data to get counts from
#' @param samples df a samples df to subset the data to
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
    # get counts
    counts <- get_counts(data, samples, normalised = FALSE)
    # add back 'count' to end of column names
    colnames(counts) <- paste(colnames(counts), 'count')
  }

  if (normalised_counts) {
    # get normalised counts and add normalised count back to col names
    norm_counts <- get_counts(data, samples, normalised = TRUE)
    colnames(norm_counts) <- paste(colnames(norm_counts), 'normalised count')
  }

  if (counts) {
    if (normalised_counts) {
      subset_data <- as_tibble(cbind(subset_metadata, counts, norm_counts))
    } else {
      subset_data <- as_tibble(cbind(subset_metadata, counts))
    }
  } else {
    if (normalised_counts) {
      subset_data <- as_tibble(cbind(subset_metadata, norm_counts))
    } else {
      subset_data <- as_tibble(subset_metadata)
    }
  }
  return(subset_data)
}
