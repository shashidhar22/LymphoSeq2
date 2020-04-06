#' Remove sequence
#' 
#' Removes an amino acid sequence and associated data from all instances within 
#' a list of data frames and then recomputes the frequencyCount.
#' 
#' @param file.list A list of data frames imported using the LymphoSeq function 
#' readImmunoSeq.  "aminoAcid", "count", and "frequencyCount" are required columns.
#' @param sequence A character vector of one or more amino acid sequences to 
#' remove from the list of data frames.
#' @return Returns a list of data frames like the one imported except all rows 
#' with the specified amino acid sequence are removed.  The frequencyCount is 
#' recalculated.
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' file.list <- readImmunoSeq(path = file.path)
#' 
#' searchSeq(list = file.list, sequence = "CASSDLIGNGKLFF")
#' 
#' cleansed <- removeSeq(file.list = file.list, sequence = "CASSDLIGNGKLFF")
#' 
#' searchSeq(list = cleansed, sequence = "CASSDLIGNGKLFF")
#' @export
#' @import  tidyverse
removeSeq <- function(file.list, sequence) {
    study_table <- study_table %>% 
                   filter(!aminoAcid %in% sequence) %>% 
                   group_by(sample) %>%
                   mutate(frequencyCount = `count` / sum(`count`)) %>% 
                   ungroup()
    return(study_table)
}