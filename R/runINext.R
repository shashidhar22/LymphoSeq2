#' Run iNEXT on repertoire_ids
#' 
#' Given a repertoire_id table, for each generate rarefaction curves to estimate 
#' repertoire diversity. The method used to generate the rarefaction curve
#' is derived from Chao et al., (2014) using the iNEXT library
#'
#' @param sample_table A tibble consisting antigen receptor sequencing 
#' data imported by the LymphoSeq function readImmunoSeq. "junction_aa", "duplicate_count", and 
#' "frequencyCount" are required columns.
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' study_table <- readImmunoSeq(path = file.path)
#' 
#' productive_aa <- productiveSeq(study_table = study_table, 
#'   aggregate = "junction_aa", prevalence = TRUE)
#'
#' sample_table <- productive_aa %>% filter(repertoire_id == "TRB_Unsorted_1320")
#'
#' @export
#' @import tidyverse
#' @import iNEXT
#' @import purrr
runINext <- function(sample_table, progress, color="repertoire_id") {
    progress$tick()$print()
    repertoire_id <- sample_table$repertoire_id[1]
    rarefaction_tables <- sample_table %>% 
                          dplyr::pivot_wider(names_from=repertoire_id, 
                                      id_cols=junction_aa, 
                                      values_from=duplicate_count, 
                                      values_fill=list(duplicate_count=0)) %>% 
                          dplyr::select(-junction_aa) %>% as.matrix()
    rarefaction_tables <- iNEXT::iNEXT(rarefaction_tables, 
                                       q=0, 
                                       datatype="abundance", 
                                       endpoint=100000, 
                                       se=TRUE, 
                                       conf=0.95, 
                                       nboot=10)
    rarefaction_tables <- rarefaction_tables$iNextEst[[repertoire_id]] %>% 
                          dplyr::as_tibble() %>% 
                          dplyr::mutate(repertoire_id = repertoire_id)
    return(rarefaction_tables)
}
