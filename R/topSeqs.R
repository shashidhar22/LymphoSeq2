#' Top sequences
#' 
#' Creates a tibble of a selected number of top productive sequences 
#' from the study table.
#' 
#' @param productive_table A tibble of productive sequences generated 
#' by the LymphoSeq2 function productiveSeq.  "duplicate_frequency" and "junction_aa" 
#' are a required columns.
#' @param top The number of top productive sequences in each data frame to subset 
#' by their frequencies.
#' @return Returns a tibble of a selected number of top productive sequences 
#' from a list of data frames.
#' @seealso \code{\link{chordDiagramVDJ}}
#' @examples
#' file_path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq2
#' stable <- readImmunoSeq(path = file_path)
#' atable <- productiveSeq(study_table = stable, aggregate = "junction_aa")
#' top_seqs <- topSeqs(productive_table = atable, top = 1)
#' @export
#' @import tidyverse dtplyr
topSeqs <- function(productive_table, top = 1) {
    top_seqs <- productive_table %>%
                dtplyr::lazy_dt() %>%
                dplyr::group_by(repertoire_id) %>%
                dplyr::arrange(desc(duplicate_frequency)) %>%
                dplyr::slice_head(n = top) %>%
                dplyr::as_tibble()
    return(top_seqs)
}
