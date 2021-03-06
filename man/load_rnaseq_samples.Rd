% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load.R
\name{load_rnaseq_samples}
\alias{load_rnaseq_samples}
\title{Load RNA-seq samples file}
\usage{
load_rnaseq_samples(samples_file)
}
\arguments{
\item{samples_file}{Name of the file to open}
}
\value{
tibble
}
\description{
\code{load_rnaseq_samples} opens a file containing sample
data and returns a tibble
It expects the file to contain columns named sample and
condition which it will make factors.
}
\examples{
data <- load_rnaseq_samples(samples_file = 'test_samples.tsv')

}
