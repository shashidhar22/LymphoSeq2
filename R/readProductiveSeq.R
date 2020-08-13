#' Productive sequences
#' 
#' \code{productiveSeq} Remove unproductive CDR3 sequences from a list of data frames.
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
#' file_path <- base::system.file("extdata", "TCRB_study", package = "LymphoSeq2")
#' 
#' study_table <- LymphoSeq2::readImmunoSeq(path = file_path)
#' 
#' productive_aa <- LymphoSeq2::productiveSeq(study_table = study_table, 
#'                                            aggregate = "junction_aa", 
#'                                            prevalence = TRUE)
#'                                            
#' @export
#' @import tidyverse
#' @import progress
productiveSeq <- function(study_table, aggregate = "junction_aa", prevalence = FALSE) {
    if (aggregate == "junction" & prevalence) {
        stop("In order to add prevalence to your list of data frames, aggregate must be equal 'junction_aa'.", 
             call. = FALSE)
    }
    nsample <- study_table %>% 
               dplyr::pull(repertoire_id) %>% 
               base::unique() %>% 
               base::length()
    progress <- dplyr::progress_estimated(nsample)
    agg_table <- study_table %>% 
                 dplyr::group_by(repertoire_id) %>%
                 dplyr::group_split() %>%
                 purrr::map(~ aggreateSeq(.x, aggregate, prevalence, progress)) %>% 
                 dplyr::bind_rows()
    progress$stop()
    return(agg_table)
}

#' Group productive sequences by repertoire
#' 
#' @param progress Dplyr progress bar
#' @inheritParams productiveSeq
aggreateSeq <- function(study_table, aggregate, prevalence, progress) {
    progress$tick()$print()
    if (aggregate == "junction") {
        study_table <- study_table %>% 
                       dplyr::filter(reading_frame == "in-frame") %>% 
                       dplyr::group_by(repertoire_id, junction) %>% 
                       dplyr::arrange(desc(duplicate_count)) %>%    
                       dplyr::summarize(junction_aa = dplyr::first(junction_aa), 
                                        duplicate_count = base::sum(duplicate_count), 
                                        v_call = dplyr::first(v_call), 
                                        reading_frame = dplyr::first(reading_frame),
                                        j_call = dplyr::first(j_call), 
                                        d_call = dplyr::first(d_call), 
                                        v_family = dplyr::first(v_family), 
                                        d_family = dplyr::first(d_family),
                                        j_family = dplyr::first(j_family)) %>%
                       dplyr::ungroup() %>%
                       dplyr::mutate(duplicate_frequency = duplicate_count / base::sum(duplicate_count))
    }  else if (aggregate == "junction_aa") {
        study_table <- study_table  %>% 
                       dplyr::filter(reading_frame == "in-frame") %>% 
                       dplyr::group_by(repertoire_id, junction_aa) %>% 
                       dplyr::arrange(desc(duplicate_count)) %>% 
                       dplyr::summarize(duplicate_count = base::sum(duplicate_count), 
                                        v_call = dplyr::first(v_call), 
                                        reading_frame = dplyr::first(reading_frame),
                                        j_call = dplyr::first(j_call), 
                                        d_call = dplyr::first(d_call), 
                                        v_family = dplyr::first(v_family), 
                                        d_family = dplyr::first(d_family),
                                        j_family = dplyr::first(j_family)) %>%
                       dplyr::ungroup() %>%
                       dplyr::mutate(duplicate_frequency = duplicate_count / base::sum(duplicate_count))
    }
    if (prevalence) {
        prev_table <- LymphoSeq2::prevalenceTRB
        study_table <- dplyr::left_join(study_table, prev_table, by=c("junction_aa" = "aminoAcid")) %>% 
                       dplyr::mutate(prevalence = tidyr::replace_na(prevalence, 0))
        
    }
    return(study_table)
    
    
}
