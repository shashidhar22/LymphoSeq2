#' Common sequences in two or more samples
#' 
#' Creates a data frame of the common sequences in two or more samples, reporting 
#' their frequencies in each.
#' 
#' @param samples A character vector of two or more sample names in 
#' productive.aa.
#' @param productive_aa A list of productive amino acid sequences generated
#' by the LymphoSeq function productiveSeq where aggregate = "aminoAcid".
#' @return Returns a data frame of the common sequences between two or more files 
#' displaying their frequencies in each.
#' @seealso \code{\link{commonSeqsVenn}}
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' study_table <- readImmunoSeq(path = file.path)
#' 
#' productive_aa <- productiveSeq(study_table = study_table, aggregate = "aminoAcid")
#' 
#' commonSeqs(samples = c("TRB_Unsorted_0", "TRB_Unsorted_32"), 
#'    productive_aa = productive_aa)
#' @export
#' @importFrom plyr llply
#' @import ggplot2
#' @import tidyverse
commonSeqs <- function(samples, productive_aa) {
    
    common_seqs <- productive_aa %>% 
                   filter(sample %in% samples) %>%
                   add_column(seen = 1) %>%
                   pivot_wider(id_cols=aminoAcid, names_from=sample, values_from=frequencyCount, values_fill=list(frequencyCount=0)) %>%
                   filter_at(vars(-aminoAcid), all_vars(. != 0)) 
    return(common_seqs)
}