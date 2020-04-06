#' Unique sequences
#' 
#' Aggregates all productive sequences within a list of data frames by count.
#' 
#' @param productive_table A tibble of productive amino acid sequences 
#' imported using the function LymphoSeq function productiveSeq where the 
#' aggregate parameter was set to "aminoAcid". 
#' @return A data frame of unique amino acid sequences from the list of 
#' data frames aggregated by count.
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' study_table <- readImmunoSeq(path = file.path)
#' 
#' productive_table <- productiveSeq(study_table = study_table, aggregate = "aminoAcid")
#' 
#' unique_seqs <- uniqueSeqs(productive_table = productive_table, unique_type = "aminoAcid")
#' @export
#' @importFrom plyr llply ldply
#' @importFrom stats aggregate
#' @import tidyverse
uniqueSeqs <- function(productive_table = productive_table, unique_type = "aminoAcid") {
    # Add checks to see if the tibble is a prudctive table
    unique_seq <- tibble()
    if (unique_type == "nucleotide") {
        unique_seq <- productive_table %>% group_by(nucleotide) %>% summarize(count = sum(count)) %>% arrange(desc(count))
    } else if (unique_type == "aminoAcid") {
        unique_seq <- productive_table %>% group_by(aminoAcid) %>% summarize(count = sum(count)) %>% arrange(desc(count))
    }
    return(unique_seq)
}