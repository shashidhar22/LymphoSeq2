#' Clonality (summary statistics)
#'
#' `clonality()` Creates a tibble giving the total number of sequences, number
#' of unique productive sequences, number of genomes, entropy, clonality, Gini
#' coefficient, TCR/BCR convergence, and the frequency of the top productive 
#' sequences for any given sample.
#'
#' @param study_table A tibble consisting of antigen receptor sequencing
#' imported by the LymphoSeq2 function [readImmunoSeq()]. "junction_aa",
#' "duplicate_count", and "duplicate_frequency" are required columns. Note that
#' clonality is usually calculated from productive junction sequences.
#' Therefore, it is not recommended to run this function using a productive
#' sequence list aggregated by amino acids.
#' @param rarefy A Boolean value
#'    * `TRUE` : Rarefied diversity metrics are calculated by sampling down each
#'    repertoire in the input table down to the repertoire with the smallest 
#'    number of sequences and calculating the diversity metrics on the sampled
#'    data. The process is repeated for the number of iterations specified by 
#'    the user and the diversity metrics are averaged over the number of 
#'    iterations. Default 100 (the default)
#'    * `FALSE` (the default): Diversity metrics will be calculated considering 
#'    the raw repertoire data for each of the samples 
#' @param iterations Number of iterations to run the sampled clonality metrics.
#' @param min_count The minimum depth to which each repertoire in the study must
#' be sampled to. Default 1000 (the default)
#' @return Returns a tibble giving the total number of sequences, number of
#' unique productive sequences, number of genomes, clonality, Gini coefficient,
#' Simpson index, inverse Simpson index, and the frequency of the top
#' productive sequence for each repertoire_id.
#' @details Clonality is derived from the Shannon entropy, which is calculated
#' from the frequencies of all productive sequences divided by the logarithm of
#' the total number of unique productive sequences.  This normalized entropy
#' value is then inverted (1 - normalized entropy) to produce the clonality
#' metric.
#'
#' The Gini coefficient is an alternative metric used to calculate repertoire
#' diversity and is derived from the Lorenz curve.  The Lorenz curve is drawn
#' such that x-axis represents the cumulative percentage of unique sequences and
#' the y-axis represents the cumulative percentage of reads.  A line passing
#' through the origin with a slope of 1 reflects equal frequencies of all clones.
#' The Gini coefficient is the ratio of the area between the line of equality
#' and the observed Lorenz curve over the total area under the line of equality.
#' Both Gini coefficient and clonality are reported on a scale from 0 to 1 where
#' 0 indicates all sequences have the same frequency and 1 indicates the
#' repertoire is dominated by a single sequence.
#'
#' TCR/BCR convergence is defined as the average number of productive CDR3
#' nucleotide sequences that form the same productive CDR3 amino acid sequence.
#'
#' Sequencing depth and amount of input available can often confound diversity
#' metrics. For example, a peripheral blood sample can appear to be more clonal
#' than a tumor sample when it is not sequenced to adequate depth. To overcome
#' this we can sample down the sample down all repertoires to the depth of the
#' sample with the least number of sequences and then calculate the diversity
#' metrics. Repeating this process multiple times and averaging the diversity
#' metrics can give a more accurate representation of sample diversity and
#' enable comparison of repertoire samples from different experiments and
#' different tissue of origin
#' @examples
#' file_path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq2")
#' study_table <- LymphoSeq2::readImmunoSeq(path = file_path)
#' study_table <- LymphoSeq2::topSeqs(study_table, top = 100)
#' raw_clonality <- LymphoSeq2::clonality(study_table)
#' sampled_clonality <- LymphoSeq2::clonality(study_table,
#'   rarefy = TRUE,
#'   iterations = 100,
#'   min_count = 100
#' )
#' @seealso [LymphoSeq2::lorenzCurve()]
#' @export
#' @import magrittr
clonality <- function(study_table, rarefy = FALSE, iterations = 100, min_count = 1000) {
    if (rarefy) {
        low_count <- study_table %>% 
            dplyr::group_by(repertoire_id) %>%
            dplyr::summarize(total = sum(duplicate_count)) %>% 
            dplyr::filter(total < min_count) %>%
            dplyr::pull(repertoire_id)
        if (length(low_count) >= 1) {
            warning(stringr::str_c("Dropping the following samples since they have less than ",
                    min_count, "sequences \n", stringr::str_c(low_count, sep =""), 
                    sep = ""))
        }
        study_table <- study_table %>%
            dplyr::filter(!(repertoire_id %in% low_count)) %>%
            dplyr::group_by(repertoire_id) %>%
            dplyr::group_split() %>%
            purrr::map(~iterativeSummary(.x, iterations, min_count)) %>%
            dplyr::bind_rows()      
    } else {
        study_table <- study_table %>%
            dplyr::group_by(repertoire_id) %>%
            dplyr::group_split() %>%
            purrr::map(summarySeq) %>%
            dplyr::bind_rows()
    }
  
  return(study_table)
}

#' Get summary statistics for each repertoire_id in the analysis
#'
#' @inheritParams clonality
#' @return Tibble summarizing the sequence information for each repertoire_id
#'
#' @export
#' @import magrittr 

summarySeq <- function(study_table) {
  productive <- LymphoSeq2::productiveSeq(study_table, aggregate = "junction")
  frequency <- productive %>%
    dplyr::pull(duplicate_frequency)
  counts <- productive %>%
    dplyr::pull(duplicate_count)
  entropy <- -base::sum(frequency * base::log2(frequency), na.rm = TRUE)
  clonality <- 1 - base::round(entropy / base::log2(base::nrow(productive)), digits = 6)
  convergence <- productive %>% 
      LymphoSeq2::topSeqs(top = 100) %>%
      dplyr::group_by(junction_aa) %>% 
      dplyr::summarise(convergence = length(unique(junction))) %>%
      dplyr::pull(convergence) %>% 
      mean()
  study_summary <- tibble::tibble(
    repertoire_id = study_table$repertoire_id[1],
    total_sequences = base::nrow(study_table),
    unique_productive_sequences = base::nrow(productive),
    total_count = base::sum(study_table$duplicate_count),
    clonality = clonality,
    gini_coefficient = ineq::Gini(productive$duplicate_frequency),
    top_productive_sequence = base::max((productive$duplicate_frequency) * 100),
    convergence = convergence
  )
  return(study_summary)
}

#' Uncount, sample and iteratively calculate repertoire summary statistics
#'
#' @inheritParams clonality
#' @return Tibble summarizing the sequence information for each repertoire_id
#' normalized by depth of sequencing
#'
#' @export
iterativeSummary <- function(study_table, iterations, min_count = 1000) {
    uncount_table <- study_table %>% 
        tidyr::uncount(weights = duplicate_count) 

    summary_table <- purrr::map(1:iterations, 
                                \(i)sampledSummary(uncount_table, min_count)) %>%
        dplyr::bind_rows() %>% 
        dplyr::group_by(repertoire_id) %>%
        dplyr::summarise_all(mean)
    return(summary_table)
}


#' Sample repertoire
#'
#' @inheritParams clonality
#' @return Tibble summarizing the sequence information for each repertoire_id
#' randomly sampled to the minimum count 
#'
#' @export
sampledSummary <- function(study_table, min_count) {
    study_table <- tibble::as_tibble(study_table) %>% 
        dplyr::sample_n(min_count) %>%
        dplyr::select(-duplicate_frequency) %>% 
        dplyr::group_by_all() %>% 
        dplyr::summarise(duplicate_count = dplyr::n()) %>%
        dplyr::mutate(duplicate_frequency = duplicate_count/sum(duplicate_count)) %>%
        dplyr::ungroup() 
    summary_table <- summarySeq(study_table) 
    return(summary_table)
}
    
        
