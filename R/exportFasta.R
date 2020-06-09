#' Export sequences in fasta format
#' 
#' Export junction or amino acid sequences in fasta format.
#' 
#' @param sample_table A tibble consisting of antigen receptor sequences 
#' imported by the LymphoSeq function readImmunoSeq.
#' @param type A character vector indicating whether "junction_aa" or "junction" sequences
#' should be exported.  If "junction_aa" is specified, then run productiveSeqs first.
#' @param names A character vector of one or more column names to name the sequences.
#' If "rank" is specified, then the rank order of the sequences by frequency is used.
#' @return Exports fasta files to the working directory.
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' study_table <- readImmunoSeq(path = file.path)
#' 
#' exportFasta(study_table = study_table, type = "junction", names = c("rank", "junction_aa", "duplicate_count"))
#' 
#' productive_aa <- productiveSeq(study_table = study_table, aggregate = "junction_aa")
#' 
#' exportFasta(list = productive_aa, type = "junction_aa", names = "duplicate_frequency")
#' @export
#' @import tidyverse
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings AAStringSet
#' @importFrom Biostrings writeXStringSet
exportFasta <- function(study_table, type = "junction", 
                        names = c("rank", "junction_aa", "duplicate_count")) {
    if (type == "junction") {
        study_table <- study_table %>% 
                       arrange(repertoire_id, desc(duplicate_frequency)) %>% 
                       rowid_to_column %>% 
                       rename(rowid = rank) %>%
                       mutate(sequences = junction) %>%
                       unite(fasta_name, names)
    } else if (type == "junction_aa") {
        study_table <- productiveSeq(study_table)
        study_table <- productive_aa %>% 
                       arrange(repertoire_id, desc(duplicate_frequency)) %>% 
                       rowid_to_column %>% 
                       rename(rowid = rank) %>%
                       mutate(sequences = junction_aa) %>%
                       unite(fasta_name, names) 
    }
    study_table %>% 
    group_by(repertoire_id) %>% 
    group_split() %>% 
    map(~writeFasta(.x, type))
    message(paste("Fasta files saved to", getwd()))
}

writeFasta <- function(sample_table, type) {
    repertoire_id <- sample_table$repertoire_id[1]
    if (type == "junction"){
        fasta <- Biostrings::DNAStringSet(sample_table$sequences)
    } else if (type == "junction_aa") {
        fasta <- Biostrings::AAStringSet(sample_table$sequences)
    }
    names(fasta) <- sample_table$fasta_name
    Biostrings::writeXStringSet(fasta, paste(repertoire_id, "fasta", sep="."))
}
