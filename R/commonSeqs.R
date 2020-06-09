#' Common sequences in two or more repertoire_ids
#' 
#' Creates a data frame of the common sequences in two or more repertoire_ids, reporting 
#' their frequencies in each.
#' 
#' @param repertoire_ids A character vector of two or more repertoire_id names in 
#' productive.aa.
#' @param productive_aa A list of productive amino acid sequences generated
#' by the LymphoSeq function productiveSeq where aggregate = "junction_aa".
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
#' @importFrom plyr llply
#' @import ggplot2
#' @import tidyverse
commonSeqs <- function(repertoire_ids, productive_aa) {
    
    common_seqs <- productive_aa %>% 
                   filter(repertoire_id %in% repertoire_ids) %>%
                   add_column(seen = 1) %>%
                   pivot_wider(id_cols=junction_aa, 
                               names_from=repertoire_id, 
                               values_from=duplicate_count, 
                               values_fill=list(duplicate_count=0)) %>%
                   filter_at(vars(-junction_aa), all_vars(. != 0)) 
    return(common_seqs)
}