% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plotTrack.R
\name{plotTrack}
\alias{plotTrack}
\title{Clone tracking alluvial plot}
\usage{
plotTrack(
  clone_table,
  alist = NULL,
  apal = NULL,
  breaks = 2,
  alphas = 0,
  breaks_pal = c("#7fc97f", "#beaed4")
)
}
\arguments{
\item{clone_table}{A tibble of productive amino acid sequences
generated by LymphoSeq2 function \code{\link[=cloneTrack]{cloneTrack()}}}

\item{alist}{An optional list of amino acid, if a list is provided only those
sequences will be highlighted}

\item{apal}{An optional list of palette colors used for the amino acids to be
highlighted}

\item{breaks}{Add an additional band behind the alluvial plot to highlight
group, specifies after which bar the break should appear}

\item{alphas}{Set alpha for the band}

\item{breaks_pal}{Specify palette for the bands}
}
\value{
An alluvial diagram tracking particular clone across samples.
}
\description{
Creates alluvial tracking amino acid frequencies across multiple samples
}
\details{
The plot is made using the package ggplot2 and can be reformatted
using ggplot2 functions.  See examples below.
}
\examples{
file_path <- system.file("extdata", "TCRB_sequencing",
 package = "LymphoSeq2")
study_table <- LymphoSeq2::readImmunoSeq(path = file_path, threads = 1)
study_table <- LymphoSeq2::topSeqs(study_table, top = 100)
amino_table <- LymphoSeq2::productiveSeq(study_table,
  aggregate = "junction_aa"
)
clone_table <- LymphoSeq2::cloneTrack(
  study_table = amino_table,
  sample_list = c(
    "TRB_CD8_949",
    "TRB_CD8_CMV_369"
  )
)
LymphoSeq2::plotTrack(clone_table)
}
