#' Find kmers and its counts
#'
#' Calculate counts of kmers in the query nucleotide sequence
#'
#' @param study_table A tibble outputed from the readImmunoSeq function.
#' @param k The length of kmers.
#' @param separate A boolean value. TRUE to separate the counts of each k-mer
#' by repertoire_ids. FALSE to show cumulative counts instead.
#' @return A tibble with the k-mer and its counts. The counts can be cumulative
#' counts of entire study_table or counts by repertoire_id.
#' @examples
#' file_path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq2")
#'
#' stable <- readImmunoSeq(path = file_path)
#'
#' kmer_table <- countKmer(study_table = stable, k = 5, separate = TRUE)
#'
#' @export
#' @import tidyverse Biostrings

countKmer <- function(study_table, k, separate = TRUE) {
    if (separate) {
        repertoire_id_names <- study_table %>%
                                dplyr::pull(repertoire_id) %>%
                                unique()

        kmer_counts <- study_table %>%
                        dplyr::group_by(repertoire_id) %>%
                        dplyr::group_split() %>%
                        purrr::map(~calculateCounts(.x, k))
        kmer_counts <- purrr::map2(kmer_counts, c(repertoire_id_names), ~dplyr::rename(.x, !!.y := Count)) %>%
                        purrr::reduce(inner_join, by = "Kmer")
        return(kmer_counts)
    } else {
        kmer_counts <- calculateCounts(study_table, k)
        return(kmer_counts)
    }
}

calculateCounts <- function(study_table, k) {
    seq <- dplyr::pull(study_table, "junction") %>%
            stats::na.omit()
    kmer_counts <- purrr::map(seq,
        function(x) Biostrings::oligonucleotideFrequency(Biostrings::DNAString(x), k))
    kmer_counts <- purrr::map(kmer_counts,
        function(x) data.frame(Kmer = names(x), Count = unname(x))) %>%
                    dplyr::bind_rows() %>%
                    dplyr::group_by(Kmer) %>%
                    dplyr::summarise(Count = sum(Count))
    return(kmer_counts)
}

#' Plot kmer distributions
#'
#' Plot kmer distributions by repertoire id
#'
#' @param kmer_table A tibble output from the countKmer function.
#' @param top The number of top kmers to show
#' @return A stacked bar chart showing kmer distributions
#' @examples
#' file_path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq2")
#'
#' stable <- readImmunoSeq(path = file_path)
#'
#' kmer_table <- countKmer(study_table = stable, k = 5, separate = TRUE)
#' 
#' kmer_distributions <- kmerPlot(kmer_table, top = 10)
#'
#' @export
#' @import tidyverse Biostrings

kmerPlot <- function(kmer_table, top = 10) {
    kmer_table <- kmer_table %>%
                    tidyr::pivot_longer(!Kmer, names_to = "repertoire_id", values_to = "count")
    rep_num <- length(unique(kmer_table$repertoire_id))
    kmer_table <- kmer_table %>%
                    dplyr::group_by(Kmer) %>%
                    dplyr::mutate(total = sum(count)) %>%
                    dplyr::arrange(desc(total)) %>%
                    head(top * rep_num)

    bar <- ggplot(kmer_table, aes(fill = repertoire_id, y = count, x = Kmer)) +
                geom_bar(position = "stack", stat = "identity") + 
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    return(bar)
}