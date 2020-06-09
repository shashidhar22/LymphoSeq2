#' Unique sequences
#' 
#' Aggregates all productive sequences within a list of data frames by duplicate_count.
#' 
#' @param productive_table A tibble of productive amino acid sequences 
#' imported using the function LymphoSeq function productiveSeq where the 
#' aggregate parameter was set to "junction_aa". 
#' @return A data frame of unique amino acid sequences from the list of 
#' data frames aggregated by duplicate_count.
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' study_table <- readImmunoSeq(path = file.path)
#' 
#' productive_table <- productiveSeq(study_table = study_table, aggregate = "junction_aa")
#' 
#' unique_seqs <- uniqueSeqs(productive_table = productive_table, unique_type = "junction_aa")
#' @export
#' @importFrom plyr llply ldply
#' @importFrom stats aggregate
#' @import tidyverse
uniqueSeqs <- function(productive_table = productive_table, unique_type = "junction_aa") {
    # Add checks to see if the tibble is a prudctive table
    unique_seq <- tibble::tibble()
    if (unique_type == "junction") {
        unique_seq <- productive_table %>% 
                      dplyr::group_by(junction) %>% 
                      dplyr::summarize(duplicate_count = sum(duplicate_count)) %>% 
                      dplyr::arrange(desc(duplicate_count))
    } else if (unique_type == "junction_aa") {
        unique_seq <- productive_table %>% 
                      dplyr::group_by(junction_aa) %>% 
                      dplyr::summarize(duplicate_count = sum(duplicate_count)) %>% 
                      dplyr::arrange(desc(duplicate_count))
    }
    return(unique_seq)
}