#' Calculate Lorenz curve
#' 
#' Calculate a Lorenz curve derived from the frequency of the amino acid sequences.
#' 
#' @param sample_table A tibble for a single repertoire_id generated using the LymphoSeq 
#' function readImmunoSeq or productiveSeq.  "duplicate_frequency" is a required column.
#' @return Returns a Lorenz curve tibble.
getLorenz <- function(sample_table) {
  repertoire_id <- sample_table$repertoire_id[1]
  lc <- ineq::Lc(sample_table$duplicate_frequency)
  lctbl <- tibble::tibble(L = lc$L, p = lc$p) %>% 
    dplyr::mutate(repertoire_id = repertoire_id)
  return(lctbl)
}
