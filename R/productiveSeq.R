#' Productive sequences
#' 
#' Remove unproductive CDR3 sequences from a list of data frames.
#' 
#' @param study_table A tibble consisting antigen receptor sequencing 
#' data imported by the LymphoSeq function readImmunoSeq. "junction_aa", "duplicate_count", and 
#' "duplicate_frequency" are required columns.
#' @param aggregate Indicates whether the values of "duplicate_count", "duplicate_frequency", 
#' and "esimatedNumberGenomes" should be aggregated by amino acid or junction 
#' sequence.  Acceptable values are "junction_aa" or "junction".  If "junction_aa" 
#' is selected, then resulting data frame will have columns corresponding to 
#' junction_aa, duplicate_count, and duplicate_frequency.  
#' If "junction" is selected then all columns in the 
#' original list will be present in the outputted list.  The difference in 
#' output is due to the fact that the same amino acid CDR3 sequence may be 
#' encoded by multiple unique junction sequences with differing V, D, and J 
#' genes.
#' @param prevalence A Boolean value indicating if a new column should be added 
#' to each of the data frames giving the prevalence of each CDR3 amino acid 
#' sequence in 55 healthy donor peripheral blood samples.  TRUE means the column 
#' is added and FALSE means it is not.  Values range from 0 to 100\% where 
#' 100\% means the sequence appeared in the blood of all 55 individuals.  The 
#' prevalenceTRB database is located in a separate package called LymphoSeqDB 
#' that should be loaded automatically.
#' @return Returns a list of data frames of productive amino acid sequences with
#' recomputed values for "duplicate_count", "duplicate_frequency", and 
#' "esimatedNumberGenomes".  A productive sequences is defined as a sequences 
#' that is in frame and does not have an early stop codon.
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' study_table <- readImmunoSeq(path = file.path)
#' 
#' productive.nt <- productiveSeq(study_table = study_table, 
#'    aggregate = "junction", prevalence = FALSE)
#' 
#' productive.aa <- productiveSeq(study_table = study_table, 
#'   aggregate = "junction_aa", prevalence = TRUE)
#' @seealso Refer to the LymphoSeqDB package for details regarding the 
#' prevalenceTRB database.
#' @export
#' @import tidyverse
#' @import LymphoSeqDB
productiveSeq <- function(study_table, aggregate = "junction_aa", prevalence = FALSE) {
    if (aggregate == "junction" & prevalence) {
        stop("In order to add prevalence to your list of data frames, aggregate must be equal 'junction_aa'.", 
             call. = FALSE)
    }
    if (aggregate == "junction") {
        study_table <- study_table %>% 
                       dplyr::filter(reading_frame == "in-frame") %>% 
                       dplyr::group_by(repertoire_id, junction) %>% 
                       dplyr::arrange(desc(duplicate_count)) %>%    
                       dplyr::summarize(junction_aa = first(junction_aa), 
                                        duplicate_count = sum(duplicate_count), 
                                        v_call = first(v_call), 
                                        reading_frame = first(reading_frame),
                                        j_call = first(j_call), 
                                        d_call = first(d_call), 
                                        v_family = first(v_family), 
                                        d_family = first(d_family),
                                        j_family = first(j_family)) %>% 
                       dplyr::mutate(duplicate_frequency = duplicate_count / sum(duplicate_count)) %>% 
                       dplyr::ungroup()
    }  else if (aggregate == "junction_aa") {
        study_table <- study_table  %>% 
                       dplyr::filter(reading_frame == "in-frame") %>% 
                       dplyr::group_by(repertoire_id, junction_aa) %>% 
                       dplyr::arrange(desc(duplicate_count)) %>% 
                       dplyr::summarize(duplicate_count = sum(duplicate_count), 
                                        v_call = first(v_call), 
                                        reading_frame = first(reading_frame),
                                        j_call = first(j_call), 
                                        d_call = first(d_call), 
                                        v_family = first(v_family), 
                                        d_family = first(d_family),
                                        j_family = first(j_family)) %>% 
                       dplyr::mutate(duplicate_frequency = duplicate_count / sum(duplicate_count)) %>% 
                       dplyr::ungroup()
    }
    if (prevalence) {
        study_table <- dplyr::left_join(study_table, LymphoSeqDB::prevalenceTRB, by="junction_aa") %>% 
                       dplyr::mutate(prevalence = replace_na(prevalence, 0))

    }
    return(study_table)
}