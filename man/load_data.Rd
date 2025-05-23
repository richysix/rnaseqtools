% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load.R
\name{load_data}
\alias{load_data}
\title{Load data}
\usage{
load_data(data_file, ...)
}
\arguments{
\item{data_file}{Name of file to load}

\item{...}{Other arguments passed on to the [readr] functions}
}
\value{
data.frame The loaded data
}
\description{
\code{load_data} is the underlying function used to load data
It automatically detects the delimiter based on the filename and uses the
appropriate [readr] function.
Possibilities are [readr::read_csv()], [readr::read_tsv()] or [readr::readr_delim()]
[readr::read_delim()] will require the delimiter passing in as an argument to
[load_rnaseq_data()] or [load_detct_data()]
}
