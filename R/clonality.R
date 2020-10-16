#' Clonality
#' 
#' Creates a tibble giving the total number of sequences, number of unique 
#' productive sequences, number of genomes, entropy, clonality, Gini 
#' coefficient, and the frequency (\%) of the top productive sequences form a repertoire_id tibble.
#' 
#' @param study_table A tibble consisting of antigen receptor 
#' sequencing imported by the LymphoSeq function readImmunoSeq. "junction_aa", "duplicate_count", 
#' and "duplicate_frequency" are required columns. Note that clonality is usually calculated from 
#' productive junction sequences. Therefore, it is not recommended to run this function using a 
#' productive sequence list aggregated by amino acids.
#' @return Returns a tibble giving the total number of sequences, number of 
#' unique productive sequences, number of genomes, clonality, Gini coefficient, 
#' and the frequency (\%) of the top productive sequence, simpson index,
#' inverse simpson index, hill diversity index, chao diversity index, and kemp diversity index
#' for each repertoire_id.
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
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq2")
#' 
#' study_table <- readImmunoSeq(path = file.path)
#' 
#' clonality(study_table = study_table)
#' @seealso \code{\link{lorenzCurve}}
#' @export
#' @importFrom ineq Gini
#' @import breakaway plyr reshape jsonlite httr vegan 
clonality <- function(study_table) {
    study_table <- study_table %>% 
                   dplyr::group_by(repertoire_id) %>% 
                   dplyr::group_split() %>% 
                   purrr::map(summarySeq) %>% 
                   dplyr::bind_rows()
    return(study_table)
}

#' Get summary statistics for each repertoire_id in the analysis
#'
#' @param sample_table immune repertoire tibble for each a repertoire_id
#'
#' @return Tibble summarizing the sequence information for each repertoire_id
#'
#' @export
#' @import tidyverse breakaway vegan

summarySeq <- function(study_table) {
    productive <- LymphoSeq2::productiveSeq(study_table, aggregate="junction")
    frequency <- productive %>% 
                 dplyr::pull(duplicate_frequency)
    counts <- productive %>% 
              dplyr::pull(duplicate_count)
    entropy <- -base::sum(frequency * base::log2(frequency), na.rm=TRUE)
    clonality <- 1 - base::round(entropy/base::log2(base::nrow(productive)), digits = 6)
    simpson_index <- vegan::diversity(frequency, index = "simpson")
    inverse_simpson <- vegan::diversity(frequency, index = "invsimpson")
    chao_estimate <- breakaway::chao1(counts)$estimate
    kemp_estimate <- breakaway::kemp(counts)$estimate
    hill_estimate <- breakaway::true_hill(frequency, q = 0)
    breakaway <- breakaway::breakaway(counts)$estimate
    model <- base::paste(breakaway::breakaway(counts)$model, "breakaway", sep = "_")
    study_summary <- tibble::tibble(repertoire_id = study_table$repertoire_id[1], 
                                    total_sequences = base::nrow(study_table), 
                                    unique_productive_sequences = base::nrow(productive),
                            total_count = base::sum(study_table$duplicate_count), 
                            clonality = clonality, 
                            simpson_index = simpson_index,
                            inverse_simpson = inverse_simpson,
                            gini_coefficient = ineq::Gini(productive$duplicate_frequency), 
                            top_productive_sequence = base::max((productive$duplicate_frequency) * 100),
                            chao_estimate = chao_estimate,
                            kemp_estimate = kemp_estimate,
                            hill_estimate = hill_estimate)
    return(study_summary)
}
