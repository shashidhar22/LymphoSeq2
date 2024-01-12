#' Bhattacharyya, Similarity, Sorensen, or PSI matrix
#'
#' Calculates the Bhattacharyya coefficient, Similarity score, Sorensen Index, or
#' Percent Similarity Index of all pairwise comparison from a list of data frames.
#'
#' @param productive_table A tibble of productive sequences generated
#' by the LymphoSeq function [productiveSeq()].  "duplicate_frequency" and
#' "junction_aa" are a required columns.
#' @param mode The mode to use for calculating pairwise similarity. Can take the
#' values "Bhattacharyya", "Similarity", "Sorensen", or "PSI". Default is
#' "Bhattacharyya".
#' @return A data frame of Bhattacharyya coefficients, Similarity scores,
#' Sorensen Index, or Percent Similarity Index calculated from all pairwise
#' comparisons from a list of repertoire_id data frames. Both metrics measure
#' the amount of overlap between two samples. The value ranges from 0 to 1 where
#' 1 indicates the sequence frequencies are identical in the two samples and 0
#' indicates no shared frequencies.
#' @examples
#' file_path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq2")
#' study_table <- LymphoSeq2::readImmunoSeq(path = file_path, threads = 1)
#' study_table <- LymphoSeq2::topSeqs(study_table, top = 100)
#' amino_table <- LymphoSeq2::productiveSeq(study_table, aggregate = "junction_aa")
#' bhattacharyya_matrix <- LymphoSeq2::scoringMatrix(
#'   productive_table = amino_table,
#'   mode = "Bhattacharyya"
#' )
#' similarity_matrix <- LymphoSeq2::scoringMatrix(
#'   productive_table = amino_table,
#'   mode = "Similarity"
#' )
#' sorensen_matrix <- LymphoSeq2::scoringMatrix(
#'   productive_table = amino_table,
#'   mode = "Sorensen"
#' )
#' psi_matrix <- LymphoSeq2::scoringMatrix(
#'   productive_table = amino_table,
#'   mode = "PSI"
#' )
#' @seealso [LymphoSeq2::pairwisePlot()] for plotting results as a heat map.
#' @export
#' @import magrittr
scoringMatrix <- function(productive_table, mode = "Bhattacharyya") {
  sample_list <- productive_table %>%
    dplyr::select(
      repertoire_id, junction_aa, duplicate_frequency,
      duplicate_count
    ) %>%
    dplyr::group_by(repertoire_id) %>%
    dplyr::group_split()
  if (mode == "Bhattacharyya") {
    scoring_matrix <- list(sample_list, sample_list) %>%
      purrr::cross() %>%
      purrr::map(bhattacharyyaCoefficient) %>%
      dplyr::bind_rows() %>%
      tidyr::pivot_wider(
        id_cols = sample1,
        names_from = sample2,
        values_from = bhattacharyya_coefficient
      )
  } else if (mode == "Similarity") {
    scoring_matrix <- list(sample_list, sample_list) %>%
      purrr::cross() %>%
      purrr::map(similarityScore) %>%
      dplyr::bind_rows() %>%
      tidyr::pivot_wider(
        id_cols = sample1,
        names_from = sample2,
        values_from = similarityScore
      )
  } else if (mode == "Sorensen") {
    scoring_matrix <- list(sample_list, sample_list) %>%
      purrr::cross() %>%
      purrr::map(sorensenIndex) %>%
      dplyr::bind_rows() %>%
      tidyr::pivot_wider(
        id_cols = sample1,
        names_from = sample2,
        values_from = sorensenIndex
      )
  } else if (mode == "PSI") {
    scoring_matrix <- list(sample_list, sample_list) %>%
      purrr::cross() %>%
      purrr::map(percentSI) %>%
      dplyr::bind_rows() %>%
      tidyr::pivot_wider(
        id_cols = sample1,
        names_from = sample2,
        values_from = percentSI
      )
  }
  row_names <- scoring_matrix$sample1
  scoring_matrix <- scoring_matrix %>%
    dplyr::select(-sample1) %>%
    as.matrix()
  rownames(scoring_matrix) <- row_names
  return(scoring_matrix)
}

#' Bhattacharyya coefficient
#'
#' Calculates the Bhattacharyya coefficient of two samples.
#'
#' @param sample_list A list of two tibble corresponding derived from the
#' `productiveSeq()`function in LymphoSeq2. "duplicate_frequency",
#' "junction_aa", and "repertoire_id" columns are necessary for the calculation
#' of the Bhattacharyya coefficient.
#'
#' @return A tibble with one row and three columns sample1, sample2,
#' bhattacharyya coefficient
#'
#' @seealso [LymphoSeq2::scoringMatrix()]
#' @import magrittr
#' @export
bhattacharyyaCoefficient <- function(sample_list) {
  sample1 <- sample_list[[1]]
  sample2 <- sample_list[[2]]
  sample_merged <- dplyr::full_join(sample1, sample2,
    by = "junction_aa",
    suffix = c("_p", "_q")
  ) %>%
    dplyr::mutate(
      duplicate_frequency_p = tidyr::replace_na(duplicate_frequency_p, 0),
      duplicate_frequency_q = tidyr::replace_na(duplicate_frequency_q, 0)
    )
  s <- sample_merged$duplicate_frequency_p * sample_merged$duplicate_frequency_q
  bc <- base::sum(base::sqrt(s))
  bhattacharyya_coefficient <- tibble::tibble(
    sample1 = sample1$repertoire_id[1],
    sample2 = sample2$repertoire_id[1],
    bhattacharyya_coefficient = bc
  )
  return(bhattacharyya_coefficient)
}
#' Similarity score
#'
#' Calculates the similarity score of two samples.
#'
#' @param sample_list A list of two tibble corresponding derived from the
#' [productiveSeq()] function in LymphoSeq2. "duplicate_frequency",
#' "junction_aa", and "repertoire_id" columns are necessary for the calculation
#' of the Bhattacharyya coefficient.
#' @return Returns the similarity score, a measure of the amount of
#' overlap between two samples.  The value ranges from 0 to 1 where 1 indicates
#' the sequence frequencies are identical in the two samples and 0
#' indicates no shared frequencies.
#' @seealso [LymphoSeq2::scoringMatrix()]
#' @import magrittr
#' @export
similarityScore <- function(sample_list) {
  sample1 <- sample_list[[1]]
  sample2 <- sample_list[[2]]
  s1 <- sample1 %>%
    dplyr::filter(junction_aa %in% sample2$junction_aa) %>%
    dplyr::summarise(total = sum(duplicate_count)) %>%
    base::as.integer()
  s2 <- sample2 %>%
    dplyr::filter(junction_aa %in% sample1$junction_aa) %>%
    dplyr::summarise(total = sum(duplicate_count)) %>%
    base::as.integer()
  score <- (s1 + s2) / (base::sum(sample1$duplicate_count) + base::sum(sample2$duplicate_count))
  similarity_score <- tibble::tibble(
    sample1 = sample1$repertoire_id[1],
    sample2 = sample2$repertoire_id[1],
    similarityScore = score
  )
  return(similarity_score)
}
#' Sorensen index
#'
#' Calculates the Sorensen index between two groups of repertoires. Similar to
#' a Jaccard index, Sorensen index gives a greater weight to shared sequences
#' over unique sequences.
#'
#' @param sample_list A list of two tibble corresponding derived from the
#' [productiveSeq()] function in LymphoSeq2. "duplicate_frequency",
#' "junction_aa", and "repertoire_id" columns are necessary for the calculation
#' of the Bhattacharyya coefficient.
#' @return Returns the similarity score, a measure of the amount of
#' overlap between two samples.  The value ranges from 0 to 1 where 1 indicates
#' the sequence frequencies are identical in the two samples and 0
#' indicates no shared frequencies.
#' @seealso [LymphoSeq2::scoringMatrix()]
#' @import magrittr
#' @export
sorensenIndex <- function(sample_list) {
  sample1 <- sample_list[[1]]
  sample2 <- sample_list[[2]]
  intersection <- dplyr::inner_join(sample1, sample2, by = "junction_aa", suffix = c("_1", "_2"))
  unique_sample1 <- dplyr::anti_join(sample1, sample2, by = "junction_aa")
  unique_sample2 <- dplyr::anti_join(sample2, sample1, by = "junction_aa")
  a <- intersection %>%
    dplyr::pull(junction_aa) %>%
    unique() %>%
    length()
  b <- unique_sample1 %>%
    dplyr::pull(junction_aa) %>%
    unique() %>%
    length()
  c <- unique_sample2 %>%
    dplyr::pull(junction_aa) %>%
    unique() %>%
    length()
  sorensen_index <- (2 * a) / ((2 * a) + b + c)
  sorensen_score <- tibble::tibble(
    sample1 = sample1$repertoire_id[1],
    sample2 = sample2$repertoire_id[1],
    sorensenIndex = sorensen_index
  )
  return(sorensen_score)
}
#' Percent similarity index
#'
#' Calculates the Percent similarity index between two groups of repertoires.
#' Percent similarity index, not only compares the number of similar and
#' dissimilar species present between two sites, but also incorporate abundance.
#'
#' @param sample_list A list of two tibble corresponding derived from the
#' [productiveSeq()] function in LymphoSeq2. "duplicate_frequency",
#' "junction_aa", and "repertoire_id" columns are necessary for the calculation
#' of the Bhattacharyya coefficient.
#' @return Returns the similarity score, a measure of the amount of
#' overlap between two samples.  The value ranges from 0 to 1 where 1 indicates
#' the sequence frequencies are identical in the two samples and 0
#' indicates no shared frequencies.
#' @seealso [LymphoSeq2::scoringMatrix()]
#' @import magrittr
#' @export
percentSI <- function(sample_list) {
  sample1 <- sample_list[[1]]
  sample2 <- sample_list[[2]]
  combined <- dplyr::full_join(sample1, sample2, by = "junction_aa", suffix = c("_1", "_2")) %>%
    dplyr::select(junction_aa, duplicate_frequency_1, duplicate_frequency_2) %>%
    dplyr::mutate(
      duplicate_frequency_1 = tidyr::replace_na(duplicate_frequency_1, 0),
      duplicate_frequency_2 = tidyr::replace_na(duplicate_frequency_2, 0)
    ) %>%
    dplyr::group_by(junction_aa, duplicate_frequency_1, duplicate_frequency_2) %>%
    dplyr::summarize(min_count = min(duplicate_frequency_1, duplicate_frequency_2)) %>%
    dplyr::ungroup()
  min_sum <- combined %>%
    dplyr::pull(min_count) %>%
    sum()
  sum_sample1 <- combined %>%
    dplyr::pull(duplicate_frequency_1) %>%
    sum()
  sum_sample2 <- combined %>%
    dplyr::pull(duplicate_frequency_2) %>%
    sum()
  percent_si <- 200 * sum(min_sum) / (sum_sample1 + sum_sample2)
  psi_score <- tibble::tibble(
    sample1 = sample1$repertoire_id[1],
    sample2 = sample2$repertoire_id[1],
    percentSI = percent_si
  )
  return(psi_score)
}
