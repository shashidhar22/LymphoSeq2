#' Top sequences
#' 
#' Creates a tibble of a selected number of top productive sequences 
#' from the study table.
#' 
#' @param productive_table A tibble of productive sequences generated 
#' by the LymphoSeq function productiveSeq.  "duplicate_frequency" and "junction_aa" 
#' are a required columns.
#' @param top The number of top productive sequences in each data frame to subset 
#' by their frequencies.
#' @return Returns a tibble of a selected number of top productive sequences 
#' from a list of data frames.
#' @seealso \code{\link{chordDiagramVDJ}}
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' study_table <- readImmunoSeq(path = file.path)
#' 
#' productive_table <- productiveSeq(study_table = study_table, aggregate = "junction_aa")
#' 
#' top_seqs <- topSeqs(productive_table = productive_table, top = 1)
#' @export
#' @import tidyverse
topSeqs <- function(productive_table, top = 1) {
    top_seqs <- productive_table %>% 
                dplyr::group_by(rerpertoire_id) %>% 
                arrange(desc(duplicate_frequency)) %>% 
                top_n(top, wt=duplicate_frequency)
    return(top_seqs)
} 