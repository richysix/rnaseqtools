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
  # data <- standardise_coltypes(data)

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
    "chr" = "f",
    "#Chr" = "f",
    "Chr" = "f",
    "start" = "i",
    "Start" = "i",
    "end" = "i",
    "End" = "i",
    "strand" = "f",
    "Strand" = "f"
  )

  # get columns
  header <- readr::read_tsv(data_file, n_max = 1000)
  # set column types
  column_types <- vector('list', length = ncol(header))
  names(column_types) <- colnames(header)
  for(colname in colnames(header)) {
    # check if the colname is in types_for_cols
    if (colname %in% names(types_for_cols)) {
      if (types_for_cols[[colname]] == "f") {
        column_types[[colname]] <- readr::col_factor()
      } else if (types_for_cols[[colname]] == "i") {
        column_types[[colname]] <- readr::col_integer()
      }
    } else if (grepl("count$", colname)) {
      if (grepl("normalised", colname)) {
        column_types[[colname]] <- readr::col_double()
      } else {
        # count columns should be integer
        column_types[[colname]] <- readr::col_integer()
      }
    } else {
      # everything else should be as parsed
      column_types[[colname]] <- readr::spec(header)$cols[[colname]]
    }
  }
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
  # Read data
  samples <-
    readr::read_tsv(samples_file,
                    col_types = readr::cols(
                      sample = readr::col_factor(),
                      condition = readr::col_factor()
                    )
    )

  return(samples)
}
