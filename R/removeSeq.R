#' Remove sequence
#' 
#' Removes an amino acid sequence and associated data from all instances within 
#' a list of data frames and then recomputes the frequencyCount.
#' 
#' @param study_table A list of data frames imported using the LymphoSeq function 
#' readImmunoSeq.  "aminoAcid", "count", and "frequencyCount" are required columns.
#' @param sequence A character vector of one or more amino acid sequences to 
#' remove from the list of data frames.
#' @return Returns a list of data frames like the one imported except all rows 
#' with the specified amino acid sequence are removed.  The frequencyCount is 
#' recalculated.
#' @examples
#' file_path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq2")
#' stable <- readImmunoSeq(path = file_path)
#' searchSeq(stable, sequence = "CASSDLIGNGKLFF")
#' cleansed <- removeSeq(stable, sequence = "CASSDLIGNGKLFF")
#' searchSeq(cleansed, sequence = "CASSDLIGNGKLFF")
#' @export
#' @import  tidyverse
removeSeq <- function(study_table, sequence) {
    study_table <- study_table %>% 
                   dplyr::filter(!junction_aa %in% sequence) %>% 
                   dplyr::group_by(repertoire_id) %>%
                   dplyr::mutate(duplicate_frequency = `duplicate_count` / sum(`duplicate_count`)) %>% 
                   dplyr::ungroup()
    return(study_table)
}
