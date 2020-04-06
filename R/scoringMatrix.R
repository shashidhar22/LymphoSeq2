#' Bhattacharyya matrix
#' 
#' Calculates the Bhattacharyya coefficient of all pairwise comparison from a 
#' list of data frames.
#' 
#' @param productive.seqs A tibble of productive sequences generated 
#' by the LymphoSeq function productiveSeq.  "frequencyCount" and "aminoAcid" 
#' are a required columns.
#' @return A data frame of Bhattacharyya coefficients or Similarity scores calculated from all 
#' pairwise comparisons from a list of sample data frames.  Both metrics 
#' measure the amount of overlap between two samples.  The 
#' value ranges from 0 to 1 where 1 indicates the sequence frequencies are 
#' identical in the two samples and 0 indicates no shared frequencies.
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' study_table <- readImmunoSeq(path = file.path)
#' 
#' productive_aa <- productiveSeq(study_table, aggregate = "aminoAcid")
#' 
#' bhattacharyya_matrix <- scoringMatrix(productive_table = productive_aa, mode = "Bhattacharyya")
#' similarity_matrix <- scoringMatrix(productive_table = productive_aa, mode = "Similarity")
#' @seealso \code{\link{pairwisePlot}} for plotting results as a heat map.
#' @export
#' @import tidyverse
scoringMatrix <- function(productive_table, mode="Bhattacharyya") {
    sample_list <- productive_table %>% 
                            select(sample, aminoAcid, frequencyCount) %>% 
                            group_by(sample) %>% 
                            group_split()
    if (mode == "Bhattacharyya") {
        scoring_matrix <- list(sample_list, sample_list) %>% 
                                cross() %>% 
                                map(bhattacharyyaCoefficient) %>%
                                reduce(rbind) %>% 
                                pivot_wider(id_cols=sample1, names_from=sample2, values_from=bhattacharyya_coefficient)
    } else if ( mode == "Similarity") {
        scoring_matrix <- list(sample_list, sample_list) %>% 
                                cross() %>% 
                                map(similarityScore) %>%
                                reduce(rbind) %>% 
                                pivot_wider(id_cols=sample1, names_from=sample2, values_from=similarityScore)
    }
    row_names <- scoring_matrix$sample1
    scoring_matrix <- scoring_matrix %>% select(-sample1) %>% as.matrix()
    rownames(scoring_matrix) <- row_names
    return(scoring_matrix) 
}

#' Bhattacharyya coefficient
#' 
#' Calculates the Bhattacharyya coefficient of two samples.
#' 
#' @param sample_list A list of two tibble corresponding derived from the productiveSeq
#' function in LymphoSeq2. `frequencyCount`, `aminoAcid`, and `sample` columns are necessary
#' for the calcualtion of the Bhattacharyya coefficient.
#' 
#' @return A tibble with one row and three columns sample1, sample2, bhattacharyyaCoefficient
#'
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' study_table <- readImmunoSeq(path = file.path)
#' 
#' productive.aa <- productiveSeq(study_table, aggregate = "aminoAcid")
#' 
#' sample_list <- productive.aa %>% 
#'                filter(sample %in% c("TRB_Unsorted_32", "TRB_Unsorted_83")) %>%
#'                group_by(sample) %>% group_split()
#' 
#' bhattacharyyaCoefficient(sample_list)
#' @seealso \code{\link{scoringMatrix}}
#' @export
bhattacharyyaCoefficient <- function(sample_list) {
    sample1 <- sample_list[[1]]
    sample2 <- sample_list[[2]]
    sample_merged <- full_join(sample1, sample2, by="aminoAcid", suffix = c("_p", "_q")) %>% 
                     mutate(frequencyCount_p = replace_na(frequencyCount_p, 0), 
                            frequencyCount_q = replace_na(frequencyCount_q, 0))              
    s <- sample_merged$frequencyCount_p * sample_merged$frequencyCount_q
    bc <- sum(sqrt(s))
    bhattacharyya_coefficient <- tibble(sample1=sample1$sample[1], 
                                        sample2=sample2$sample[1],
                                        bhattacharyya_coefficient=bc)
    return(bhattacharyya_coefficient)
}
#' Similarity score
#' 
#' Calculates the similarity score of two samples.
#' 
#' @param sample1 A data frame consisting of frequencies of antigen receptor 
#' sequences.  "aminoAcid" and "count" are a required columns.
#' @param sample2 A data frame consisting of frequencies of antigen receptor 
#' sequences.  "aminoAcid" and "count" are a required columns.
#' @return Returns the similarity score, a measure of the amount of 
#' overlap between two samples.  The value ranges from 0 to 1 where 1 indicates 
#' the sequence frequencies are identical in the two samples and 0 
#' indicates no shared frequencies.
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' study_table <- readImmunoSeq(path = file.path)
#' 
#' productive.aa <- productiveSeq(study_table, aggregate = "aminoAcid")
#' 
#' sample_list <- productive.aa %>% 
#'                filter(sample %in% c("TRB_Unsorted_32", "TRB_Unsorted_83")) %>%
#'                group_by(sample) %>% group_split()
#' 
#' similarityScore(sample_list)
#' @seealso \code{\link{scoringMatrix}}
#' @export
similarityScore <- function(sample_list) {
    sample1 <- sample_list[[1]]
    sample2 <- sample_list[[2]]
    s1 <- sample1 %>% filter(aminoAcid %in% sample2$aminoAcid) %>% summarise(total = sum(`count`)) %>% as.integer()
    s2 <- sample2 %>% filter(aminoAcid %in% sample1$aminoAcid) %>% summarise(total = sum(`count`)) %>% as.integer()
    score <- (sample1 + sample2)/ (sum(sample1$`count`) + sum(sample2$`count`))
    similarity_score <- tibble(sample1=sample1$sample[1], sample2=sample2$sample[1], similarityScore=score)
    return(score)
}