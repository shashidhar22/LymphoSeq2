#' Export sequences in fasta format
#'
#' Export junction or amino acid sequences in fasta format.
#'
#' @param study_table A tibble consisting of antigen receptor sequences
#' imported by the LymphoSeq function [readImmunoSeq()].
#' @param type A character vector indicating whether "junction_aa" or "junction" 
#' sequences should be exported.  If "junction_aa" is specified, then run 
#' [LymphoSeq2::productiveSeq()] first.
#' @param names A character vector of one or more column names to name the 
#' sequences.If "rank" is specified, then the rank order of the sequences by 
#' frequency is used.
#' @return Exports fasta files to the working directory.
#' @examples
#' file_path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq2")
#' # Export raw data
#' study_table <- LymphoSeq2::readImmunoSeq(path = file_path)
#' LymphoSeq2::exportFasta(study_table = study_table, type = "junction", 
#'   names = c("junction_aa", "duplicate_count"))
#' # Export only productive junction amino acid sequences
#' amino_table <- LymphoSeq2::productiveSeq(study_table = study_table, 
#'   aggregate = "junction_aa")
#' LymphoSeq2::exportFasta(study_table = amino_table, type = "junction_aa", 
#'   names = "duplicate_frequency")
#' @export
#' @import magrittr
exportFasta <- function(study_table, type = "junction",
                        names = c("rank", "junction_aa", "duplicate_count")) {
  if (type == "junction") {
    study_table <- study_table %>%
      dplyr::arrange(repertoire_id, dplyr::desc(duplicate_frequency)) %>%
      tibble::rowid_to_column() %>%
      dplyr::mutate(sequences = junction) %>%
      tidyr::unite(fasta_name, names)
  } else if (type == "junction_aa") {
    study_table <- LymphoSeq2::productiveSeq(study_table)
    study_table <- study_table %>%
      dplyr::arrange(repertoire_id, dplyr::desc(duplicate_frequency)) %>%
      tibble::rowid_to_column() %>%
      dplyr::mutate(sequences = junction_aa) %>%
      tidyr::unite(fasta_name, names)
  }
  study_table %>%
    dplyr::group_by(repertoire_id) %>%
    dplyr::group_split() %>%
    purrr::map(~ writeFasta(.x, type))
  message(paste("Fasta files saved to", getwd()))
}
#' Write FASTA file
#' 
#' @inheritParams exportFasta
writeFasta <- function(study_table, type) {
  repertoire_id <- study_table$repertoire_id[1]
  if (type == "junction") {
    fasta <- Biostrings::DNAStringSet(study_table$sequences)
  } else if (type == "junction_aa") {
    fasta <- Biostrings::AAStringSet(study_table$sequences)
  }
  names(fasta) <- study_table$fasta_name
  Biostrings::writeXStringSet(fasta, paste(repertoire_id, "fasta", sep = "."))
}
