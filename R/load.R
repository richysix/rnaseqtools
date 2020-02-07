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
#' data <- load_rnaseq_data(datafile = 'all.tsv')
#'
#' @export
load_rnaseq_data <- function(data_file) {
  # Read data
  data <- read.delim(data_file, header=TRUE, check.names=FALSE)

  # Support different column names
  names(data)[names(data) == 'chr']               <- 'Chr'
  names(data)[names(data) == '#Chr']              <- 'Chr'
  names(data)[names(data) == 'start']             <- 'Start'
  names(data)[names(data) == 'end']               <- 'End'
  names(data)[names(data) == 'strand']            <- 'Strand'
  names(data)[names(data) == 'ID']                             <- 'GeneID'
  names(data)[ grepl("e[0-9]+ Ensembl Gene ID", names(data)) ] <- 'GeneID'
  names(data)[names(data) == 'Gene ID']                        <- 'GeneID'
  names(data)[names(data) == 'adjpval']           <- 'adjp'
  names(data)[names(data) == 'padj']              <- 'adjp'
  names(data)[names(data) == 'Adjusted p value']  <- 'adjp'
  names(data)[names(data) == 'Gene name']         <- 'Name'

  return(data)
}
