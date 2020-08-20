#' Common sequences in two or more repertoire_ids
#' 
#' Creates a data frame of the common sequences in two or more repertoire_ids, reporting 
#' their frequencies in each.
#' 
#' @param study_table A list of productive amino acid sequences generated
#' by the LymphoSeq function productiveSeq where aggregate = "junction_aa".
#' @param repertoire_ids A character vector of two or more repertoire_id names in 
#' productive.aa.
#' @return Returns a data frame of the common sequences between two or more files 
#' displaying their frequencies in each.
#' @seealso \code{\link{commonSeqsVenn}}
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' study_table <- readImmunoSeq(path = file.path)
#' 
#' productive_aa <- productiveSeq(study_table = study_table, aggregate = "junction_aa")
#' 
#' commonSeqs(repertoire_ids = c("TRB_Unsorted_0", "TRB_Unsorted_32"), 
#'    productive_aa = productive_aa)
#' @export
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
