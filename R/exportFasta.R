#' Export sequences in fasta format
#' 
#' Export nucleotide or amino acid sequences in fasta format.
#' 
#' @param samlpe_table A tibble consisting of antigen receptor sequences 
#' imported by the LymphoSeq function readImmunoSeq.
#' @param type A character vector indicating whether "aminoAcid" or "nucleotide" sequences
#' should be exported.  If "aminoAcid" is specified, then run productiveSeqs first.
#' @param names A character vector of one or more column names to name the sequences.
#' If "rank" is specified, then the rank order of the sequences by frequency is used.
#' @return Exports fasta files to the working directory.
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' study_table <- readImmunoSeq(path = file.path)
#' 
#' exportFasta(study_table = study_table, type = "nucleotide", names = c("rank", "aminoAcid", "count"))
#' 
#' productive_aa <- productiveSeq(study_table = study_table, aggregate = "aminoAcid")
#' 
#' exportFasta(list = productive_aa, type = "aminoAcid", names = "frequencyCount")
#' @export
#' @import tidyverse
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings AAStringSet
#' @importFrom Biostrings writeXStringSet
exportFasta <- function(study_table, type = "nucleotide", names = c("rank", "aminoAcid", "count")) {
    if (type == "nucleotide") {
        study_table <- study_table %>% 
                       arrange(sample, desc(frequencyCount)) %>% 
                       rowid_to_column %>% 
                       rename(rowid = rank) %>%
                       mutate(sequences = nucleotide) %>%
                       unite(fasta_name, names)
    } else if (type == "aminoAcid") {
        study_table <- productiveSeq(study_table)
        study_table <- productive_aa %>% 
                       arrange(sample, desc(frequencyCount)) %>% 
                       rowid_to_column %>% 
                       rename(rowid = rank) %>%
                       mutate(sequences = aminoAcid) %>%
                       unite(fasta_name, names) 
    }
    study_table %>% group_by(sample) %>% group_split() %>% map(~writeFasta(.x, type))
    message(paste("Fasta files saved to", getwd()))
}

writeFasta <- function(sample_table, type) {
    sample <- sample_table$sample[1]
    if (type == "nucleotide"){
        fasta <- Biostrings::DNAStringSet(sample_table$sequences)
    } else if (type == "aminoAcid") {
        fasta <- Biostrings::AAStringSet(sample_table$sequences)
    }
    names(fasta) <- sample_table$fasta_name
    Biostrings::writeXStringSet(fasta, paste(sample, "fasta", sep="."))
}
