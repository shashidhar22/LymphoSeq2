#' Differential abundance analysis
#' 
#' Use a Fisher exact test to calculate differential abdunance of each sequence in 
#' two samples and reports the log2 transformed fold change, P value and adjusted P value.
#' 
#' @param study_table A tibble consisting of antigen receptor sequences 
#' imported by the LymphoSeq function readImmunoSeq.
#' @param repertoire_ids A character vector of two repertoire_ids in study_table 
#' to be compared. If NULL, the first two repertoire_ids from study_table will be used.
#' @param abundance The input value for the Fisher exact test.  "duplicate_count"
#' is recommend but "duplicate_count" may also be used.
#' @param type A character vector indicating whether "junction_aa" or "junction" sequences
#' should be used.  If "junction_aa" is specified, then run productiveSeqs first.
#' @param q A numeric value between 0.0 and 1.0 indicating the threshold Holms adjusted 
#' P value (also knowns as the false discovery rate or q value) to subset the results with.  
#' Any sequences with a q value greater than this value will not be shown.
#' @param zero A numeric value to set all zero values to when calculating the log2
#' transformed fold change between samples 1 and 2.  This does not apply to the 
#' p and q value calculations.
#' @param parallel A boolean indicating wheter parallel processing should be used or not.
#' @return Returns a data frame with columns corresponding to the frequency of the abudance
#' measure in samples 1 and 2, the P value, Q value (Holms adjusted P value, also knowns as 
#' the false discovery rate), and log2 transformed fold change.
#' @examples
#' file_path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq2")
#' 
#' stable <- readImmunoSeq(path = file_path)
#' 
#' atable <- productiveSeq(study_table = stable, aggregate = "junction_aa")
#' 
#' differentialAbundance(study_table = atable, repertoire_ids = c("TRB_Unsorted_949", "TRB_Unsorted_1320"), 
#'                       type = "junction_aa", q = 0.01, zero = 0.001)
#' @export
#' @import tidyverse
differentialAbundance <- function(study_table, repertoire_ids = NULL, abundance = "duplicate_count", 
                                  type = "junction_aa", q = 1, zero = 1, parallel = FALSE) {
    if (base::is.null(repertoire_ids)) {
        repertoire_ids <- study_table %>% 
                          dplyr::pull(repertoire_id) 
        repertoire_ids <- repertoire_ids[1:3]
    }
    fisher_table <- study_table %>%
                    dplyr::filter(repertoire_id %in% repertoire_ids) %>%
                    LymphoSeq2::seqMatrix(by = "duplicate_count") %>%
                    dplyr::mutate(not_x = base::sum(!!base::as.name(repertoire_ids[1])) - !!base::as.name(repertoire_ids[1]),
                                  not_y = base::sum(!!base::as.name(repertoire_ids[2])) - !!base::as.name(repertoire_ids[2])) %>%
              
                    dplyr::rowwise() %>%
                    dplyr::mutate(fisher = list(LymphoSeq2::fisherFunction(!!base::as.name(repertoire_ids[1]), 
                                                                       !!base::as.name(repertoire_ids[2]), 
                                                                       not_x, 
                                                                       not_y)))  %>%
                    tidyr::unnest_wider(fisher) %>%
                    dplyr::select(junction_aa,
                                 !!base::as.name(repertoire_ids[1]),
                                 !!base::as.name(repertoire_ids[2]),
                                 p, q, l2fc)
    
    
    return(fisher_table)
}

#' 
#' @export
fisherFunction <- function(x, y, not_x, not_y) {
    matrix <- matrix(c(x, y, not_x, not_y), nrow = 2)
    fisher <- stats::fisher.test(matrix, workspace = 2e6)
    q <- stats::p.adjust(fisher$p, method = "holm")
    l2fc <- base::log2(x/y)
    return(list( "p" = fisher$p, "q" = q, "l2fc" = l2fc))
}
