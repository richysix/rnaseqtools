#' Get count data for an RNA-seq dataset
#'
#' \code{get_counts} returns a data.frame of counts,
#' either raw or normalised (default)
#'
#' @param data df rna-seq data to get counts from
#' @param samples df a samples df to subset the data to
#' @param normalised logical indicatig whether to return
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
    count_data <- count_data[, samples$sample ]
  }

  return(count_data)
}
