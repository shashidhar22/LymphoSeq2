#' Unique sequences
#'
#' Aggregates all productive sequences within a list of data frames by duplicate_count.
#'
#' @param productive_table A tibble of productive amino acid sequences
#' imported using the function LymphoSeq2 function [productiveSeq()] where the
#' aggregate parameter was set to "junction_aa".
#' @param unique_type Use `"junction_aa"` (the default) to aggregate by amino acid sequences.
#' Use `"junction"` to aggregate by nucleotide sequences.
#' @return A data frame of unique amino acid sequences from the list of
#' data frames aggregated by "duplicate_count"
#' @examples
#' file_path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq2")
#' study_table <- LymphoSeq2::readImmunoSeq(path = file_path)
#' amino_table <- LymphoSeq2::productiveSeq(study_table = study_table, aggregate = "junction_aa")
#' unique_seqs <- LymphoSeq2::uniqueSeqs(productive_table = amino_table, unique_type = "junction_aa")
#' @export
#' @import magrittr
uniqueSeqs <- function(productive_table = productive_table, unique_type = "junction_aa") {
  # Add checks to see if the tibble is a productive table
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
