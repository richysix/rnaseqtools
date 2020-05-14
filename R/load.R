#' Load RNA-seq count data in from file
#'
#' \code{load_rnaseq_data} opens a file containing count data
#' and returns a data.frame with column names adjusted to
#' take account of different column name formats.
#'
#' @param data_file char, Name of the file to open
#'
#' @return data.frame
#'
#' @examples
#' data <- load_rnaseq_data(data_file = 'all.tsv')
#'
#' @export
load_rnaseq_data <- function(data_file) {
  # Read data
  data <- read.delim(data_file, header=TRUE, check.names=FALSE)
  data <- standardise_colnames(data)
  data <- standardise_coltypes(data)

  return(data)
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

#' Standardise column classes
#'
#' \code{standardise_coltypes} takes a dataframe
#' loaded from a RNA-seq all file and creates standard
#' column types. For example, Chr should be a factor even
#' though some times all the values will be integers
#'
#' @param data data.frame, data to standardise
#'
#' @return data.frame
#'
#' @examples
#' data <- standardise_coltypes(data)
#'
#' @export
standardise_coltypes <- function(data) {
  # make Chr a factor
  data$Chr <- factor(data$Chr)
  # make sure start and end are integer
  data$Start <- as.integer(data$Start)
  data$End <- as.integer(data$End)
  # make Strand a factor
  data$Strand <- factor(data$Strand)

  return(data)
}

