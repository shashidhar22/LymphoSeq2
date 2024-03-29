% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getLorenz.R
\name{getLorenz}
\alias{getLorenz}
\title{Calculate Lorenz curve}
\usage{
getLorenz(sample_table)
}
\arguments{
\item{sample_table}{A tibble for a single repertoire_id generated using the
LymphoSeq function \code{\link[=readImmunoSeq]{readImmunoSeq()}} or \code{\link[=productiveSeq]{productiveSeq()}}.
The table should contain the "duplicate_frequency" column.}
}
\value{
Returns a Lorenz curve tibble.
}
\description{
Calculate a Lorenz curve derived from the frequency of the amino acid
sequences.
}
