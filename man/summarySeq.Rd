% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clonality.R
\name{summarySeq}
\alias{summarySeq}
\title{Get summary statistics for each repertoire_id in the analysis}
\usage{
summarySeq(study_table)
}
\arguments{
\item{study_table}{A tibble consisting of antigen receptor sequencing
imported by the LymphoSeq2 function \code{\link[=readImmunoSeq]{readImmunoSeq()}}. "junction_aa",
"duplicate_count", and "duplicate_frequency" are required columns. Note that
clonality is usually calculated from productive junction sequences.
Therefore, it is not recommended to run this function using a productive
sequence list aggregated by amino acids.}
}
\value{
Tibble summarizing the sequence information for each repertoire_id
}
\description{
Get summary statistics for each repertoire_id in the analysis
}
