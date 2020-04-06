#' Clonal relatedness
#' 
#' Calculates the clonal relatedness for each sample in a study.
#' 
#' @param study_table A tibble of unproductive or productive nucleotide 
#' sequences or productive nucleotide sequences.  Nucleotide and count are 
#' required columns.
#' @param editDistance An integer giving the minimum edit distance that the 
#' sequence must be less than or equal to. See details below.
#' @details  Clonal relatedness is the proportion of nucleotide sequences that
#' are related by a defined edit distance threshold.  The value ranges from 0 to
#' 1 where 0 indicates no sequences are related and 1 indicates all sequences 
#' are related.
#' 
#' Edit distance is a way of quantifying how dissimilar two sequences 
#' are to one another by counting the minimum number of operations required to 
#' transform one sequence into the other. For example, an edit distance of 0 
#' means the sequences are identical and an edit distance of 1 indicates that 
#' the sequences different by a single amino acid or nucleotide.
#' @return Returns a tibble with the calculated clonal relatedness for each sample.
#' @examples
#' file.path <- system.file("extdata", "IGH_sequencing", package = "LymphoSeq")
#' 
#' study_table <- readImmunoSeq(path = file.path)
#' 
#' clonal_relatedness <- clonalRelatedness(study_table = study_table, editDistance = 10)
#' 
#' # Merge results with clonality table
#' clonality <- clonality(file.list = file.list)
#' merged <- merge(clonality, clonal.relatedness)
#' 
#' @export
#' @importFrom stringdist stringdist
clonalRelatedness <- function(study_table, editDistance = 10){
    clonal_relatedness  <- study_table %>%
                           group_by(sample) %>%
                           group_split() %>%
                           map(~getRelatedness(.x, editDistance)) %>% 
                           reduce(rbind)
}

#' Calculate relatedness
#' 
#' Calculates the clonal relatedness of a sample.
#' 
#' @param sample_table A tibble of unproductive or productive nucleotide 
#' sequences or productive nucleotide sequences.  Nucleotide and count are 
#' required columns.
#' @param editDistance An integer giving the minimum edit distance that the 
#' sequence must be less than or equal to. See details below.
#' 
#' @export
#' @importFrom stringdist stringdist 
#' @import tidyverse  
getRelatedness <- function(sample_table, editDistance=10) {
    sample <- sample_table$sample[1]
    top_seq <- sample_table %>% 
               arrange(desc(count)) %>% 
               top_n(1, wt=count) %>% 
               select(nucleotide) %>% 
               as.character()
    seq_distance <- stringdist::stringdist(top_seq, sample_table$nucleotide)
    related_seq <- seq_distance[seq_distance <= editDistance]
    relatedness <- length(related_seq)/nrow(sample_table)
    related_tbl <- tibble(sample=sample, clonalRelatedness=relatedness)
    return(related_tbl)
}
