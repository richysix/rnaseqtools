% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load.R
\name{load_detct_data}
\alias{load_detct_data}
\title{Load DeTCT count data in from file}
\usage{
load_detct_data(data_file, ...)
}
\arguments{
\item{data_file}{char, Name of the file to open}

\item{...}{Arguments passed to \code{\link[readr]{read_tsv}}}
}
\value{
tibble
}
\description{
\code{load_detct_data} opens a file containing count data
and returns a tibble with column names adjusted to
take account of different column name formats.
If it doesn't already exist, it creates a column called RegionID
made up of Chr:RegionStart:RegionEnd:3PrimeEndPosition:3PrimeEndStrand:GeneID
}
\examples{
data <- load_detct_data(data_file = 'all.tsv')

}
