#' Productive sequences
#' 
#' Remove unproductive CDR3 sequences from a list of data frames.
#' 
#' @param study_table A tibble consisting antigen receptor sequencing 
#' data imported by the LymphoSeq function readImmunoSeq. "aminoAcid", "count", and 
#' "frequencyCount" are required columns.
#' @param aggregate Indicates whether the values of "count", "frequencyCount", 
#' and "esimatedNumberGenomes" should be aggregated by amino acid or nucleotide 
#' sequence.  Acceptable values are "aminoAcid" or "nucleotide".  If "aminoAcid" 
#' is selected, then resulting data frame will have columns corresponding to 
#' aminoAcid, count, frequencyCount, and estimatedNumberGenomes.  
#' If "nucleotide" is selected then all columns in the 
#' original list will be present in the outputted list.  The difference in 
#' output is due to the fact that the same amino acid CDR3 sequence may be 
#' encoded by multiple unique nucleotide sequences with differing V, D, and J 
#' genes.
#' @param prevalence A Boolean value indicating if a new column should be added 
#' to each of the data frames giving the prevalence of each CDR3 amino acid 
#' sequence in 55 healthy donor peripheral blood samples.  TRUE means the column 
#' is added and FALSE means it is not.  Values range from 0 to 100\% where 
#' 100\% means the sequence appeared in the blood of all 55 individuals.  The 
#' prevalenceTRB database is located in a separate package called LymphoSeqDB 
#' that should be loaded automatically.
#' @return Returns a list of data frames of productive amino acid sequences with
#' recomputed values for "count", "frequencyCount", and 
#' "esimatedNumberGenomes".  A productive sequences is defined as a sequences 
#' that is in frame and does not have an early stop codon.
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' study_table <- readImmunoSeq(path = file.path)
#' 
#' productive.nt <- productiveSeq(study_table = study_table, 
#'    aggregate = "nucleotide", prevalence = FALSE)
#' 
#' productive.aa <- productiveSeq(study_table = study_table, 
#'   aggregate = "aminoAcid", prevalence = TRUE)
#' @seealso Refer to the LymphoSeqDB package for details regarding the 
#' prevalenceTRB database.
#' @export
#' @import tidyverse
#' @import LymphoSeqDB
productiveSeq <- function(study_table, aggregate = "aminoAcid", prevalence = FALSE) {
    if (aggregate == "nucleotide" & prevalence) {
        stop("In order to add prevalence to your list of data frames, aggregate must be equal 'aminoAcid'.", call. = FALSE)
    }
    if (aggregate == "nucleotide") {
        study_table <- study_table %>% filter(`function` == "in-frame") %>% group_by( sample, nucleotide) %>% arrange(desc(`count`)) %>%    
            summarize(aminoAcid = first(aminoAcid), `count` = sum(`count`), vGeneName = first(vGeneName), `function` = first(`function`),
            jGeneName = first(jGeneName), dGeneName = first(dGeneName), vFamilyName = first(vFamilyName), dFamilyName = first(dFamilyName),
            jFamilyName = first(jFamilyName)) %>% mutate(frequencyCount = `count` / sum(`count`)) %>% ungroup()
    }  else if (aggregate == "aminoAcid") {
        study_table <- study_table  %>% filter(`function` == "in-frame") %>% group_by(sample, aminoAcid) %>% arrange(desc(`count`)) %>% 
            summarize(`count` = sum(`count`), vGeneName = first(vGeneName), `function` = first(`function`),
            jGeneName = first(jGeneName), dGeneName = first(dGeneName), vFamilyName = first(vFamilyName), dFamilyName = first(dFamilyName),
            jFamilyName = first(jFamilyName)) %>% mutate(frequencyCount = `count` / sum(`count`)) %>% ungroup()
    }
    if (prevalence) {
        study_table <- left_join(study_table, LymphoSeqDB::prevalenceTRB, by="aminoAcid") %>% mutate(prevalence = replace_na(prevalence, 0))

    }
    return(study_table)
}