#' Top sequences
#'
#' Creates a tibble of a selected number of top productive sequences
#' from the study table.
#'
#' @param productive_table A tibble of productive sequences generated
#' by the LymphoSeq2 function [productiveSeq()]. "duplicate_frequency" and
#' "junction_aa" are a required columns.
#' @param top The number of top productive sequences in each data frame to
#'  subset by their frequencies.
#' @return Returns a tibble of a selected number of top productive sequences
#' from a list of data frames.
#' @seealso [LymphoSeq2::chordDiagramVDJ()]
#' @examples
#' file_path <- system.file("extdata", "TCRB_sequencing",
#'  package = "LymphoSeq2")
#' study_table <- LymphoSeq2::readImmunoSeq(path = file_path, threads = 1)
#' amino_table <- LymphoSeq2::productiveSeq(study_table = study_table,
#'  aggregate = "junction_aa")
#' top_seqs <- LymphoSeq2::topSeqs(productive_table = amino_table, top = 1)
#' @export
topSeqs <- function(productive_table, top = 1) {
  top_seqs <- productive_table |>
    dplyr::group_by(repertoire_id) |>
    dplyr::arrange(desc(duplicate_frequency)) |>
    dplyr::slice_head(n = top) |>
    dplyr::as_tibble()
  return(top_seqs)
}
