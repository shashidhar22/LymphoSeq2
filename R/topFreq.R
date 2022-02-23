#' Top frequencies
#' 
#' Creates a data frame of the top productive amino acid sequences that have a 
#' specified minimum frequency threshold and reports the number of samples that 
#' the sequence appears in along with the minimum, maximum, and mean frequency 
#' across all samples.  For T cell receptor beta sequences, the \% prevalence 
#' and antigen specificity of that sequence is also provided.
#' 
#' @param productive_table A tibble of productive amino acid sequences 
#' imported using the function LymphoSeq function productiveSeq where the 
#' aggregate parameter was set to "junction_aa".  
#' @param frequency The minimum frequency that the sequence appears in any of 
#' the listed samples.
#' @return A data frame of amino acid sequences and the number of samples that 
#' the sequence appears in along with the minimum, maximum, and mean frequency 
#' across all samples.  
#' For T cell receptor beta sequences, additionally reported is the 
#' \% prevalence that the sequence appears in 55 healthy donor blood samples.  
#' Also provided is the antigen specificity of that sequence if known by 
#' comparing it to a database of previously reported sequences in the 
#' literature. 
#' @examples
#' file_path <- system.file("extdata", "TCRB_sequencing", 
#'  package = "LymphoSeq2")
#' stable <- readImmunoSeq(path = file_path)
#' atable <- productiveSeq(study_table = stable, aggregate = "junction_aa")
#' top_freq <- topFreq(productive_table = atable, frequency = 0.1)
#' @export
topFreq <- function(productive_table, frequency = 0.1) {
    productive_table <- productive_table
    top_freq <- productive_table %>%
                dplyr::filter(duplicate_frequency >= frequency) %>%
                dplyr::group_by(junction_aa) %>%
                dplyr::summarise(minFrequency = min(duplicate_frequency),
                                 maxFrequency = max(duplicate_frequency),
                                 meanFrequency = mean(duplicate_frequency),
                                 numberSamples = length(duplicate_frequency > 0)) %>%
                dplyr::arrange(desc(numberSamples), desc(meanFrequency))
    
    top_freq <- dplyr::left_join(top_freq, LymphoSeq2::prevalenceTRB, by=c("junction_aa" = "aminoAcid")) %>% 
                dplyr::mutate(prevalence = tidyr::replace_na(0))
    antigen_table <- LymphoSeq2::publishedTRB %>% 
        dplyr::as_tibble() %>% 
        dplyr::select(aminoAcid, antigen)
    top_freq <- dplyr::left_join(top_freq, antigen_table, by=c("junction_aa" = "aminoAcid"))  
    return(top_freq)
}