#' Read ImmunoSeq metadata files
#'
#' Imports tab-separated value (.tsv) files containing metadata information
#' and converts them into the MiAIRR format for metadata. Additional information
#' is stored in its original format
#' 
#' @details May import tab-delimited files containing antigen receptor 
#' sequencing from with the following header set.  
#' 
#' | MiAIRR field ID | Type1           | Type2                       | Type3               | Type4                          | Type5                               | Type6             |
#' | --------------- | --------------- | --------------------------- | ------------------- | ------------------------------ | ----------------------------------- | ----------------- |
#' | junction        | nucleotide      | nucleotide                  | nucleotide          | nucleotide.CDR3.in.lowercase . | nucleotide( CDR3   in   lowercase ) | nSeqCDR3          |
#' | junction_aa     | aminoAcid       | aminoAcid                   | aminoAcid           | aminoAcid.CDR3.in.lowercase .  | aminoAcid( CDR3   in   lowercase )  | aaSeqCDR3         |
#' | duplicate_count | count ( reads ) | count ( templates / reads ) | count ( templates ) | cloneCount                     | cloneCount                          | cloneCount        |
#' | v_call          | vGeneName       | vGeneName                   | vGeneName           | vGene                          | vGene                               | allVHitsWithScore |
#' | d_call          | dGeneName       | dGeneName                   | dGeneName           | dGene                          | dGene                               | allDHitsWithScore |
#' | j_call          | jGeneName       | jGeneName                   | jGeneName           | jGene                          | jGene                               | allJHitsWithScore |
#'
#' @md
#' @name readImmunoSeq
NULL
#' @describeIn readImmunoSeq Read a set of files 
#' @param path Path to the directory containing tab-delimited files.  Only 
#' files with the extension .tsv are imported.  The names of the data frames are 
#' the same as names of the files.
#' 
#' @return Returns a tibble with MiAIRR headers and repertoire_id
#'
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' study_table <- readImmunoSeq(path = file.path, 
#'                              recursive = FALSE)
#' @export
#' @import tidyverse dplyr
#' @rdname readImmunoSeq