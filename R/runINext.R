#' Run iNEXT on repertoire_ids
#'
#' Given a repertoire_id table, for each generate rarefaction curves to estimate
#' repertoire diversity. The method used to generate the rarefaction curve
#' is derived from Chao et al., (2014) using the iNEXT library
#'
#' @param sample_table A tibble consisting antigen receptor sequencing
#' data imported by the LymphoSeq2 function readImmunoSeq. "junction_aa",
#' "duplicate_count", and "duplicate_frequency" are required columns.
#' @examples
#' file_path <- system.file("extdata", "TCRB_sequencing",
#'  package = "LymphoSeq2")
#' study_table <- LymphoSeq2::readImmunoSeq(path = file_path, threads = 1) |>
#'   LymphoSeq2::topSeqs(top = 100)
#' amino_table <- LymphoSeq2::productiveSeq(study_table,
#'   aggregate = "junction_aa",
#'   prevalence = TRUE
#' )
#' amino_table <- amino_table |>
#'   dplyr::filter(repertoire_id == "TRB_Unsorted_1320")
#' rarefaction_table <- LymphoSeq2::runINext(amino_table)
#' @export
runINext <- function(sample_table) {
  repertoire_id <- sample_table$repertoire_id[1]
  rarefaction_tables <- sample_table |>
    dplyr::group_by(junction_aa, repertoire_id) |>
    dplyr::summarise(duplicate_count = sum(duplicate_count)) |>
    tidyr::pivot_wider(
      names_from = repertoire_id,
      id_cols = junction_aa,
      values_from = duplicate_count,
      values_fill = list(duplicate_count = 0)
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-junction_aa) |>
    as.matrix()

  rarefaction_tables <- iNEXT::iNEXT(rarefaction_tables,
    q = 0,
    datatype = "abundance",
    endpoint = 100000,
    se = TRUE,
    conf = 0.95,
    nboot = 10
  )
  rarefaction_tables <- rarefaction_tables$iNextEst$size_based |>
    dplyr::as_tibble() |>
    dplyr::rename(repertoire_id = Assemblage)
  return(rarefaction_tables)
}
