% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/exportFasta.R
\name{exportFasta}
\alias{exportFasta}
\title{Export sequences in fasta format}
\usage{
exportFasta(
  study_table,
  type = "junction",
  names = c("rank", "junction_aa", "duplicate_count")
)
}
\arguments{
\item{study_table}{A tibble consisting of antigen receptor sequences
imported by the LymphoSeq function \code{\link[=readImmunoSeq]{readImmunoSeq()}}.}

\item{type}{A character vector indicating whether "junction_aa" or "junction"
sequences should be exported.  If "junction_aa" is specified, then run
\code{\link[=productiveSeq]{productiveSeq()}} first.}

\item{names}{A character vector of one or more column names to name the
sequences.If "rank" is specified, then the rank order of the sequences by
frequency is used.}
}
\value{
Exports fasta files to the working directory.
}
\description{
Export junction or amino acid sequences in fasta format.
}
\examples{
file_path <- system.file("extdata", "TCRB_sequencing",
 package = "LymphoSeq2")
# Export raw data
study_table <- LymphoSeq2::readImmunoSeq(path = file_path, threads = 1)
study_table <- LymphoSeq2::topSeqs(study_table, top = 100)
LymphoSeq2::exportFasta(
  study_table = study_table, type = "junction",
  names = c("junction_aa", "duplicate_count")
)
# Export only productive junction amino acid sequences
amino_table <- LymphoSeq2::productiveSeq(
  study_table = study_table,
  aggregate = "junction_aa"
)
LymphoSeq2::exportFasta(
  study_table = amino_table, type = "junction_aa",
  names = "duplicate_frequency"
)
}
