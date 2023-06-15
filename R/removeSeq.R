#' Remove sequence
#'
#' Removes an amino acid sequence and associated data from all instances within
#' study table and then recomputes the duplicate_frequency.
#'
#' @param study_table A tibble imported using the LymphoSeq2 function
#' [readImmunoSeq()]. "junction_aa", "duplicate_count", and
#' "duplicate_frequency" are required columns.
#' @param sequence A character vector of one or more amino acid sequences to
#' remove from the study table
#' @return Returns a tibble like the one imported except all rows
#' with the specified amino acid sequence are removed.  The
#' "duplicate_frequency" is recalculated.
#' @examples
#' library(magrittr)
#' file_path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq2")
#' study_table <- LymphoSeq2::readImmunoSeq(path = file_path)
#' Lstudy_table <- LymphoSeq2::topSeqs(study_table, top = 100)
#' LymphoSeq2::searchSeq(study_table, sequence = "CASSDLIGNGKLFF")
#' cleaned_table <- LymphoSeq2::removeSeq(study_table, sequence = "CASSDLIGNGKLFF")
#' LymphoSeq2::searchSeq(cleaned_table, sequence = "CASSDLIGNGKLFF")
#' @export
#' @import  magrittr
removeSeq <- function(study_table, sequence) {
  study_table <- study_table %>%
    dplyr::filter(!junction_aa %in% sequence) %>%
    dplyr::group_by(repertoire_id) %>%
    dplyr::mutate(duplicate_frequency = `duplicate_count` / sum(`duplicate_count`)) %>%
    dplyr::ungroup()
  return(study_table)
}
