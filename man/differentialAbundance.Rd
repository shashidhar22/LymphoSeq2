% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/differentialAbundance.R
\name{differentialAbundance}
\alias{differentialAbundance}
\title{Differential abundance analysis}
\usage{
differentialAbundance(
  study_table,
  repertoire_ids = NULL,
  abundance = "duplicate_count",
  type = "junction_aa",
  q = 1,
  zero = 1,
  parallel = FALSE
)
}
\arguments{
\item{study_table}{A tibble consisting of antigen receptor sequences
imported by the LymphoSeq2 function \code{\link[=readImmunoSeq]{readImmunoSeq()}}.}

\item{repertoire_ids}{A character vector of two repertoire_ids in study_table
to be compared. If \code{NULL} (the default), the first two repertoire_ids from
study_table will be used.}

\item{abundance}{The input value for the Fisher exact test. "duplicate_count"
is the default value and is also the recommended value.}

\item{type}{A character vector indicating whether "junction_aa" (the default)
or "junction" sequences should be used.  If "junction_aa" is specified, then
run \code{\link[=productiveSeq]{productiveSeq()}} first.}

\item{q}{A numeric value between 0.0 and 1.0 indicating the threshold Holms
adjusted P value (also known as the false discovery rate or q value) to
subset the results. Any sequences with a q value greater than this value
will not be shown.}

\item{zero}{A numeric value to set all zero values to when calculating the
log2 transformed fold change between samples 1 and 2. This does not apply to
the p and q value calculations.}

\item{parallel}{A Boolean value
\itemize{
\item \code{TRUE} : Enable parallel processing
\item \code{FALSE} (the default): Disable parallel processing
}}
}
\value{
Returns a data frame with columns corresponding to the frequency of
the abundance measure in samples 1 and 2, the P value, Q value
(Holms adjusted P value, also known as the false discovery rate), and log2
transformed fold change.
}
\description{
Use a Fisher exact test to calculate differential abundance of each sequence
in two samples and reports the log2 transformed fold change, P value and
adjusted P value.
}
\examples{
file_path <- system.file("extdata", "TCRB_sequencing",
 package = "LymphoSeq2")
study_table <- LymphoSeq2::readImmunoSeq(path = file_path, threads = 1)
study_table <- LymphoSeq2::topSeqs(study_table, top = 100)
amino_table <- LymphoSeq2::productiveSeq(
  study_table = study_table,
  aggregate = "junction_aa"
)
LymphoSeq2::differentialAbundance(
  study_table = amino_table,
  repertoire_ids = c("TRB_Unsorted_949", "TRB_Unsorted_1320"),
  type = "junction_aa", q = 0.01, zero = 0.001
)
}
