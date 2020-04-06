#' Clonality
#' 
#' Creates a tibble giving the total number of sequences, number of unique 
#' productive sequences, number of genomes, entropy, clonality, Gini 
#' coefficient, and the frequency (\%) of the top productive sequences form a sample tibble.
#' 
#' @param study_table A tibble consisting of antigen receptor 
#' sequencing imported by the LymphoSeq function readImmunoSeq. "aminoAcid", "count", 
#' and "frequencyCount" are required columns.  "estimatedNumberGenomes" is optional.  
#' Note that clonality is usually calculated from productive nucleotide sequences.  
#' Therefore, it is not recommended to run this function using a productive sequence
#' list aggregated by amino acids.
#' @return Returns a tibble giving the total number of sequences, number of 
#' unique productive sequences, number of genomes, clonality, Gini coefficient, 
#' and the frequency (\%) of the top productive sequence in each sample.
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
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' study_table <- readImmunoSeq(path = file.path)
#' 
#' clonality(study_table = study_table)
#' @seealso \code{\link{lorenzCurve}}
#' @export
#' @importFrom ineq Gini
clonality <- function(study_table) {
    study_table <- study_table %>% group_by(sample) %>% group_split() %>% map(summarySeq) %>% reduce(rbind)
    return(study_table)
}

#' Get summary statistics for each sample in the analysis
#'
#' @param sample_table immune repertoire tibble for each a sample
#'
#' @return Tibble summarizing the sequence infomration for each sample
#'
#' @export
#' @import tidyverse

summarySeq <- function(study_table) {
    productive <- productiveSeq(study_table, aggregate="nucleotide")
    frequency <- productive$frequencyCount
    entropy <- -sum(frequency * log2(frequency), na.rm=TRUE)
    clonality <- 1 - round(entropy/log2(nrow(productive)), digits = 6)
    study_summary <- tibble(sample = study_table$sample[1], totalSequences = nrow(study_table), uniqueProductiveSequences = nrow(productive),
                            totalCount = sum(study_table$count), clonality = clonality, 
                            giniCoefficient = ineq::Gini(productive$frequencyCount), topProductiveSequence = max((productive$frequencyCount) * 100),
                            estimatedNumberGenomes = sum(study_table$estimatedNumberGenomes))
    return(study_summary)
}
