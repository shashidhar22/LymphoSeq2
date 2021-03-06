#' Top frequencies
#' 
#' Creates a data frame of the top productive amino acid sequences that have a 
#' specified minimum frequency threshold and reports the number of samples that 
#' the sequence appears in along with the minimum, maximum, and mean frequency 
#' across all samples.  For T cell receptor beta sequences, the \% prevalence 
#' and antigen specificity of that sequence is also provided.
#' 
#' @param productive_aa A tibble of productive amino acid sequences 
#' imported using the function LymphoSeq function productiveSeq where the 
#' aggregate parameter was set to "junction_aa".  
#' @param percent The minimum \% frequency that the sequence appears in any of 
#' the listed samples.
#' @return A data frame of amino acid sequences and the number of samples that 
#' the sequence appears in along with the minimum, maximum, and mean frequency 
#' across all samples.  
#' For T cell receptor beta sequences, additionally reported is the 
#' \% prevalence that the sequence appears in 55 healthy donor blood samples.  
#' Also provided is the antigen specificity of that sequence if known by 
#' comparing it to a database of previously reported sequences in the 
#' literature.  The prevalenceTRB and publishedTRB databases are located in a 
#' separate package called LymphoSeqDB that should be loaded automatically.
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' study_table <- readImmunoSeq(path = file.path)
#' 
#' productive_aaa <- productiveSeq(study_table = study_table, aggregate = "junction_aa")
#' 
#' top_freq <- topFreq(productive_table = productive_aa, percent = 0.1)
#' @seealso Refer to the LymphoSeqDB package for details regarding the 
#' prevalenceTRB and publishedTRB database.
#' @export
#' @import LymphoSeqDB
#' @import tidyverse
#' @importFrom dplyr group_by summarise
topFreq <- function(productive_table, percent = 0.1) {
    top_aa <- productive_table %>%  
              dplyr::arrange(desc(duplicate_frequency)) %>% 
              dplyr::top_frac(percent, wt=duplicate_frequency) %>% 
              dplyr::as_vector()
    top_freq <- productive_table %>% 
                dplyr::filter(junction_aa %in% top_aa) %>%
                dplyr::group_by(junction_aa) %>%
                dplyr::summarise(minFrequency = min(duplicate_frequency),
                                 maxFrequency = max(duplicate_frequency),
                                 meanFrequency = mean(duplicate_frequency),
                                 numberSamples = length(duplicate_frequency > 0)) %>%
                dplyr::arrange(desc(numberSamples), desc(meanFrequency))
    

    top_freq <- dplyr::left_join(top_freq, LymphoSeqDB::prevalenceTRB, by="junction_aa") %>% 
                dplyr::mutate(prevalence = dplyr::replace_na(0))
    antigen_table <- LymphoSeqDB::publishedTRB %>% 
                     dplyr::as_tibble() %>% 
                     dplyr::select(junction_aa, antigen)
    top_freq <- dplyr::left_join(top_freq, antigen_table, by="junction_aa")  
    return(top_freq)
}