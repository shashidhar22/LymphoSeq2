#' Run iNEXT on samples
#' 
#' Given a sample table, for each generate rarefaction curves to estimate 
#' repertoire diversity. The method used to generate the rarefaction curve
#' is derived from Chao et al., (2014) using the iNEXT library
#'
#' @param sample_table A tibble consisting antigen receptor sequencing 
#' data imported by the LymphoSeq function readImmunoSeq. "aminoAcid", "count", and 
#' "frequencyCount" are required columns.
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' study_table <- readImmunoSeq(path = file.path)
#' 
#' productive_aa <- productiveSeq(study_table = study_table, 
#'   aggregate = "aminoAcid", prevalence = TRUE)
#'
#' sample_table <- productive_aa %>% filter(sample == "TRB_Unsorted_1320")
#'
#' @export
#' @import tidyverse
#' @import iNEXT
#' @import purrr
runINext <- function(sample_table, progress) {
    progress$tick()$print()
    sample = sample_table$sample[1]
    rarefaction_tables <- sample_table %>% pivot_wider(names_from=sample, id_cols=aminoAcid, values_from=count, values_fill=list(count=0)) %>% select(-aminoAcid) %>% as.matrix()
    rarefaction_tables <- iNEXT(rarefaction_tables, q=0, datatype="abundance", endpoint=100000, se=TRUE, conf=0.95, nboot=10)
    rarefaction_tables <- rarefaction_tables$iNextEst[[sample]] %>% as_tibble() %>% add_column(sample = sample)
    return(rarefaction_tables)
}
