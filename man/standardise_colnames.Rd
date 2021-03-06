% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load.R
\name{standardise_colnames}
\alias{standardise_colnames}
\title{Standardise count file column names}
\usage{
standardise_colnames(data)
}
\arguments{
\item{data}{data.frame, data to standardise}
}
\value{
data.frame
}
\description{
\code{standardise_colnames} takes a dataframe
loaded from a RNA-seq all file and creates standard
column names.
}
\examples{
data <- standardise_colnames(data)

}
