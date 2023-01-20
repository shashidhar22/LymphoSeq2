#' Clonal relatedness
#'
#' Calculates the clonal relatedness for each repertoire in a study.
#'
#' @param study_table A tibble of raw or productive nucleotide junction sequences.
#' `junction` and `duplicate_count` are required columns.
#' @param edit_distance An integer giving the minimum edit distance that the
#' sequence must be less than or equal to. See details below.
#' @details  Clonal relatedness is the proportion of junction sequences that
#' are related by a defined edit distance threshold.  The value ranges from 0 to
#' 1 where 0 indicates no sequences are related and 1 indicates all sequences
#' are related.
#'
#' Edit distance is a way of quantifying how dissimilar two sequences
#' are to one another by counting the minimum number of operations required to
#' transform one sequence into the other. For example, an edit distance of 0
#' means the sequences are identical and an edit distance of 1 indicates that
#' the sequences differ by a single amino acid or junction.
#' @return Returns a tibble with the calculated clonal relatedness for each
#' repertoire_id.
#' @examples
#' file_path <- system.file("extdata", "IGH_sequencing", package = "LymphoSeq2")
#'
#' study_table <- readImmunoSeq(path = file_path)
#'
#' clonal_relatedness <- clonalRelatedness(study_table, edit_distance = 10)
#'
#' # Merge results with clonality table
#' clonality <- clonality(study_table)
#' merged <- dplyr::full_join(clonality, clonal_relatedness, by = "repertoire_id")
#'
#' @export
#' @import magrittr
clonalRelatedness <- function(study_table, edit_distance = 10) {
  clonal_relatedness <- study_table %>%
    dplyr::group_by(repertoire_id) %>%
    dplyr::group_split() %>%
    purrr::map(~ getRelatedness(.x, edit_distance)) %>%
    dplyr::bind_rows()
  return(clonal_relatedness)
}

#' Calculate relatedness
#'
#' Calculates the clonal relatedness of a repertoire.
#'
#' @inheritParams clonalRelatedness
#' @export
#' @import magrittr
getRelatedness <- function(sample_table, edit_distance = 10) {
  top_seq <- sample_table %>%
    topSeqs(top = 1) %>%
    dplyr::pull(junction) %>%
    unique()
  seq_distance <- sample_table %>%
    dplyr::mutate(ed = stringdist::stringdist(top_seq, junction))
  related_seq <- seq_distance %>%
    dplyr::filter(ed <= edit_distance)
  relatedness <- related_seq %>%
    dplyr::group_by(repertoire_id) %>%
    dplyr::summarize(relatedness = length(unique(junction)) / base::nrow(sample_table))
  return(relatedness)
}
