% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/topFreq.R
\name{topFreq}
\alias{topFreq}
\title{Top frequencies}
\usage{
topFreq(productive_table, percent = 0.1)
}
\arguments{
\item{percent}{The minimum \% frequency that the sequence appears in any of 
the listed samples.}

\item{productive_aa}{A tibble of productive amino acid sequences 
imported using the function LymphoSeq function productiveSeq where the 
aggregate parameter was set to "aminoAcid".}
}
\value{
A data frame of amino acid sequences and the number of samples that 
the sequence appears in along with the minimum, maximum, and mean frequency 
across all samples.  
For T cell receptor beta sequences, additionally reported is the 
\% prevalence that the sequence appears in 55 healthy donor blood samples.  
Also provided is the antigen specificity of that sequence if known by 
comparing it to a database of previously reported sequences in the 
literature.  The prevalenceTRB and publishedTRB databases are located in a 
separate package called LymphoSeqDB that should be loaded automatically.
}
\description{
Creates a data frame of the top productive amino acid sequences that have a 
specified minimum frequency threshold and reports the number of samples that 
the sequence appears in along with the minimum, maximum, and mean frequency 
across all samples.  For T cell receptor beta sequences, the \% prevalence 
and antigen specificity of that sequence is also provided.
}
\examples{
file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")

study_table <- readImmunoSeq(path = file.path)

productive_aaa <- productiveSeq(study_table = study_table, aggregate = "aminoAcid")

top_freq <- topFreq(productive_table = productive_aa, percent = 0.1)
}
\seealso{
Refer to the LymphoSeqDB package for details regarding the 
prevalenceTRB and publishedTRB database.
}
