#' Common sequences scatter plot
#'
#' Creates a scatter plot of just the sequences in common between two samples.
#'
#' @param sample1 A name of a repertoire_id in a list of data frames generated
#'  by the LymphoSeq2 function [productiveSeq()].
#' @param sample2 A name of a repertoire_id in a list of data frames generated
#'  by the LymphoSeq function [productiveSeq()].
#' @param amino_table A tibble of productive amino acid sequences produced by
#'  the function [productiveSeq()] containing the samples to be compared.
#' @param show A character vector specifying whether only the common sequences
#'  should be shown or all sequences.  Available options are "common" or "all".
#' @return Returns a frequency scatter plot of two samples showing only the
#' shared sequences.
#' @details The plot is made using the package ggplot2 and can be reformatted
#' using ggplot2 functions. See examples below.
#' @seealso An excellent resource for examples on how to reformat a ggplot can
#'  be found in the R Graphics Cookbook online
#' (\url{http://www.cookbook-r.com/Graphs/}).
#' @seealso [LymphoSeq2::productiveSeq()], [LymphoSeq2::commonSeqs()],
#' [LymphoSeq2::commonSeqsVenn()], [LymphoSeq2::commonSeqsBar()]
#' @examples
#' file_path <- system.file("extdata", "TCRB_sequencing",
#'  package = "LymphoSeq2")
#' study_table <- readImmunoSeq(path = file_path, threads = 1)
#' study_table <- LymphoSeq2::topSeqs(study_table, top = 100)
#' amino_table <- productiveSeq(study_table = study_table,
#'  aggregate = "junction_aa")
#' commonSeqsPlot("TRB_Unsorted_32", "TRB_Unsorted_83",
#'   amino_table = amino_table
#' )
#' # Change the X and Y axis to log-10 scale
#' commonSeqsPlot("TRB_Unsorted_32", "TRB_Unsorted_83",
#'   amino_table = amino_table
#' ) +
#'   ggplot2::scale_x_log10() +
#'   ggplot2::scale_y_log10() +
#'   ggplot2::annotation_logticks(sides = "bl")
#' @export
commonSeqsPlot <- function(sample1, sample2, amino_table, show = "common") {
  # Check if tibble is contains unproductive sequences
  if (show == "common") {
    common <- LymphoSeq2::commonSeqs(
      repertoire_ids = c(sample1, sample2),
      study_table = amino_table
    )
    plot <- ggplot2::ggplot(
      data = common,
      ggplot2::aes_string(
        x = as.name(names(common)[2]),
        y = as.name(names(common)[3]), label = "junction_aa"
      )
    ) +
      ggplot2::geom_point() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = paste(sample1, "frequency"),
        y = paste(sample2, "frequency")
      )
  }
  if (show == "all") {
    stable1 <- amino_table |>
      dplyr::filter(repertoire_id == sample1)
    stable2 <- amino_table |>
      dplyr::filter(repertoire_id == sample2)
    all <- dplyr::full_join(stable1, stable2, by = "junction_aa") |>
      dplyr::mutate_all(~ replace(., is.na(.), 0))
    all <- all |>
      dplyr::recode(
        duplicate_frequency.x = sample1,
        duplicate_frequency.y = sample2
      )
    plot <- ggplot2::ggplot(
      data = all,
      ggplot2::aes_string(x = sample1, y = sample2, label = "junction_aa")
    ) +
      ggplot2::geom_point() +
      ggplot2::theme_minimal() +
      ggplot2::labs(
        x = paste(sample1, "frequency"),
        y = paste(sample2, "frequency")
      )
  }
  return(plot)
}
