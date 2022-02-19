#' Find kmers and its counts
#'
#' Calculate counts of kmers in the query nucleotide sequence
#'
#' @param study_table A tibble outputed from the readImmunoSeq function.
#' @param k The length of kmers.
#' @return A tibble with two columns: the kmers and its respective count.
#' @examples
#' file_path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq2")
#'
#' stable <- readImmunoSeq(path = file_path)
#'
#' countKmer(study_table = stable, k = 5)
#'
#' @export
#' @import tidyverse Biostrings

countKmer <- function(study_table, k) {
    seq <- dplyr::pull(study_table, "junction") %>%
            stats::na.omit()
    kmer_counts <- purrr::map(seq,
        function(x) Biostrings::oligonucleotideFrequency(Biostrings::DNAString(x), k))
    kmer_counts <- purrr::map(kmer_counts, function(x) x[x > 0])
    kmer_counts <- purrr::map(kmer_counts,
        function(x) data.frame(Kmer = names(x), Count = unname(x))) %>%
                    dplyr::bind_rows() %>%
                    dplyr::group_by(Kmer) %>%
                    dplyr::summarise(Count = sum(Count))
    return(kmer_counts)
}