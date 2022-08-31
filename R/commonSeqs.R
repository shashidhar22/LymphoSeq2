#' Common sequences in two or more repertoire_ids
#' 
#' Creates a data frame of the common sequences in two or more repertoire_ids, reporting 
#' their frequencies in each.
#' 
#' @param study_table A list of productive amino acid sequences generated
#' by the LymphoSeq2 function productiveSeq where aggregate = "junction_aa".
#' @param repertoire_ids A character vector of two or more repertoire_id names in 
#' study_table.
#' @return Returns a data frame of the common sequences between two or more files 
#' displaying their frequencies in each.
#' @seealso \code{\link{productiveSeq}} \code{\link{commonSeqsVenn}} \code{\link{commonSeqsPlot}} \code{\link{commonSeqsBar}}
#' @examples
#' file_path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq2")
#' 
#' stable <- readImmunoSeq(path = file_path)
#' 
#' atable <- productiveSeq(study_table = stable, aggregate = "junction_aa")
#' 
#' commonSeqs(repertoire_ids = c("TRB_Unsorted_0", "TRB_Unsorted_32"),
#'    study_table = atable)
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
