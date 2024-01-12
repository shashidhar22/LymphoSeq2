#' Common sequences in two or more repertoire_ids
#'
#' Creates a data frame of the common sequences in two or more repertoire_ids,
#' reporting their frequencies in each.
#'
#' @param study_table A list of productive amino acid sequences generated
#' by the LymphoSeq2 function [productiveSeq()] where aggregate = "junction_aa".
#' @param repertoire_ids A character vector of two or more repertoire_id names in
#' study_table.
#' @return Returns a data frame of the common sequences between two or more files
#' displaying their frequencies in each.
#' @seealso [LymphoSeq2::productiveSeq()], [LymphoSeq2::commonSeqsVenn()],
#' [LymphoSeq2::commonSeqsPlot()], [LymphoSeq2::commonSeqsBar()]
#' @examples
#' file_path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq2")
#' study_table <- LymphoSeq2::readImmunoSeq(path = file_path, threads = 1)
#' study_table <- LymphoSeq2::topSeqs(study_table, top = 100)
#' amino_table <- LymphoSeq2::productiveSeq(
#'   study_table = study_table,
#'   aggregate = "junction_aa"
#' )
#' LymphoSeq2::commonSeqs(
#'   repertoire_ids = c("TRB_Unsorted_0", "TRB_Unsorted_32"),
#'   study_table = amino_table
#' )
#' @export
#' @import magrittr
commonSeqs <- function(study_table, repertoire_ids = NULL) {
  if (base::is.null(repertoire_ids)) {
    repertoire_ids <- study_table %>%
      dplyr::pull(repertoire_id) %>%
      base::unique()
  }
  common_seqs <- study_table %>%
    LymphoSeq2::cloneTrack(sample_list = repertoire_ids) %>%
    dplyr::filter(seen > 1) %>%
    LymphoSeq2::seqMatrix()
  return(common_seqs)
}
