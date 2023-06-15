#' Select productive sequences
#'
#'@description
#' [productiveSeq()] Select productive nucleotide/amino acid CDR3 sequences
#' from a tibble containing raw AIRR formatted data. Aggregation of the raw data
#' is either done on the productive CDR3 amino acid sequence (junction_aa) or
#' the productive CDR3 nucleotide sequence (junction).  If "junction_aa"
#' is selected, then resulting tibble will display the most frequently observed.
#' V, D, J gene that were associated with the formation of the productive CDR3
#' amino acid sequence. If "junction" is selected then all columns in the
#' original list will be present in the outputted list.  The difference in
#' output is due to the fact that the same amino acid CDR3 sequence may be
#' encoded by multiple unique junction sequences with differing V, D, and J
#' genes.
#'
#' @param study_table A tibble consisting antigen receptor sequencing
#' data imported by the LymphoSeq2 function [readImmunoSeq()]. "junction_aa",
#' "duplicate_count", and "duplicate_frequency" are required columns.
#' @param aggregate Indicates whether the values of "duplicate_count" and
#' "duplicate_frequency" should be aggregated by amino acid or junction sequence.
#' Acceptable values are "junction_aa" or "junction".
#' @param prevalence A Boolean value 
#'  * `TRUE` : Add a new column the study table giving the prevalence of each CDR3 amino acid
#' sequence in 55 healthy donor peripheral blood samples.  
#'  * `FALSE` (the default): Do not add prevelance information
#' @return Returns a list of data frames of productive amino acid sequences with
#' recomputed values for "duplicate_count", "duplicate_frequency".
#' A productive sequences is defined as a sequences
#' that is in frame and does not have an early stop codon.
#' @examples
#' file_path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq2")
#' study_table <- LymphoSeq2::readImmunoSeq(path = file_path)
#' study_table <- LymphoSeq2::topSeqs(study_table, top = 100)
#' amino_table <- LymphoSeq2::productiveSeq(
#'   study_table = study_table,
#'   aggregate = "junction_aa",
#'   prevalence = TRUE
#' )
#' nucleotide_table <- LymphoSeq2::productiveSeq(
#'   study_table = study_table,
#'   aggregate = "junction",
#'   prevalence = FALSE
#' )
#' @export
#' @import magrittr
productiveSeq <- function(study_table, aggregate = "junction_aa", prevalence = FALSE) {
  if (aggregate == "junction" & prevalence) {
    stop("In order to add prevalence to your list of data frames, aggregate must be equal 'junction_aa'.",
      call. = FALSE
    )
  }
  nsample <- study_table %>%
    dplyr::pull(repertoire_id) %>%
    base::unique() %>%
    base::length()
  progress_bar <- progress::progress_bar$new(
    format = "Subsetting productive sequences [:bar] :percent eta: :eta",
    total = nsample, clear = FALSE, width = 60
  )
  progress_bar$tick(0)
  agg_table <- study_table %>%
    dplyr::group_by(repertoire_id) %>%
    dplyr::group_split() %>%
    purrr::map(~ aggreateSeq(.x, aggregate, prevalence, progress_bar)) %>%
    dplyr::bind_rows()
  return(agg_table)
}
#' Group productive sequences by repertoire
#' 
#' @keywords internal
#' @inheritParams productiveSeq
#' @param progress_bar Progress bar
#' @noRd
aggreateSeq <- function(study_table, aggregate, prevalence, progress_bar) {
  progress_bar$tick()
  if (aggregate == "junction") {
    study_table <- study_table %>%
      dplyr::filter(reading_frame == "in-frame") %>%
      dplyr::mutate(
        vdj_comb_call = stringr::str_glue("{v_call};{j_call};{d_call}"),
        vdj_comb_family = stringr::str_glue("{v_family};{j_family};{d_family}")
      ) %>%
      dtplyr::lazy_dt()
    study_table <- study_table %>%
      dplyr::group_by(junction, vdj_comb_call) %>%
      dplyr::mutate(vdj_comb_count = sum(duplicate_count)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(junction) %>%
      dplyr::arrange(desc(duplicate_count), desc(vdj_comb_count)) %>%
      dplyr::summarize(
        repertoire_id = dplyr::first(repertoire_id),
        junction_aa = dplyr::first(junction_aa),
        duplicate_count = base::sum(duplicate_count),
        v_call = dplyr::first(v_call),
        reading_frame = dplyr::first(reading_frame),
        j_call = dplyr::first(j_call),
        d_call = dplyr::first(d_call),
        v_family = dplyr::first(v_family),
        d_family = dplyr::first(d_family),
        j_family = dplyr::first(j_family)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(duplicate_frequency = duplicate_count / base::sum(duplicate_count)) %>%
      dplyr::select(
        repertoire_id, junction, junction_aa, v_call, d_call, j_call, v_family,
        d_family, j_family, reading_frame, duplicate_count, duplicate_frequency
      ) %>%
      dplyr::as_tibble()
  } else if (aggregate == "junction_aa") {
    study_table <- study_table %>%
      dplyr::filter(reading_frame == "in-frame") %>%
      dplyr::mutate(
        vdj_comb_call = stringr::str_glue("{v_call};{j_call};{d_call}"),
        vdj_comb_family = stringr::str_glue("{v_family};{j_family};{d_family}")
      ) %>%
      dtplyr::lazy_dt()
    study_table <- study_table %>%
      dplyr::group_by(junction_aa, vdj_comb_call) %>%
      dplyr::mutate(vdj_comb_count = sum(duplicate_count)) %>%
      dplyr::ungroup() %>%
      dplyr::group_by(junction_aa) %>%
      dplyr::arrange(desc(duplicate_count), desc(vdj_comb_count)) %>%
      dplyr::summarize(
        repertoire_id = dplyr::first(repertoire_id),
        duplicate_count = base::sum(duplicate_count),
        v_call = dplyr::first(v_call),
        reading_frame = dplyr::first(reading_frame),
        j_call = dplyr::first(j_call),
        d_call = dplyr::first(d_call),
        v_family = dplyr::first(v_family),
        d_family = dplyr::first(d_family),
        j_family = dplyr::first(j_family)
      ) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(duplicate_frequency = duplicate_count / base::sum(duplicate_count)) %>%
      dplyr::select(
        repertoire_id, junction_aa, v_call, d_call, j_call, v_family, d_family,
        j_family, reading_frame, duplicate_count, duplicate_frequency
      ) %>%
      dplyr::as_tibble()
  }
  if (prevalence) {
    prev_table <- LymphoSeq2::prevalenceTRB
    study_table <- dplyr::left_join(study_table, prev_table, by = c("junction_aa" = "aminoAcid")) %>%
      dplyr::mutate(prevalence = tidyr::replace_na(prevalence, 0))
  }
  return(study_table)
}
