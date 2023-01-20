#' Cumulative frequency bar plot of top sequences
#'
#' Create a cumulative frequency bar plot of a specified number of top
#' sequences.
#'
#' @param study_table A study tibble imported using the LymphoSeq2 function
#' `readImmunoSeq()` or `productiveSeq()`
#' @param top The number of top sequences to be colored in the bar plot.  All
#' other, less frequent sequences are colored violet.
#' @return Returns a cumulative frequency bar plot of the top sequences.
#' @details The plot is made using the package ggplot2 and can be reformatted
#' using ggplot2 functions.  See examples below.
#' @seealso An excellent resource for examples on how to reformat a ggplot can
#' be found in the R Graphics Cookbook online (\url{http://www.cookbook-r.com/Graphs/}).
#' @examples
#' file_path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq2")
#' study_table <- readImmunoSeq(path = file_path)
#' amino_table <- productiveSeq(study_table = study_table, aggregate = "junction_aa")
#' topSeqsPlot(study_table = amino_table, top = 10)
#' # Display the number of sequences at the top of bar plot and add a title
#' n <- as.character(nrow(study_table))
#' topSeqsPlot(study_table = amino_table, top = 10) +
#'   ggplot2::annotate("text", x = 1:length(n), y = 105, label = n, color = "black") +
#'   ggplot2::expand_limits(y = c(0, 110)) + ggplot2::ggtitle("Top sequences") +
#'   ggplot2::scale_x_discrete(limits = names(n))
#' @export
#' @import magrittr
topSeqsPlot <- function(study_table, top = 10) {
  dominant <- study_table %>%
    dplyr::group_by(repertoire_id) %>%
    dplyr::arrange(desc(duplicate_frequency)) %>%
    dplyr::slice_head(n = top) %>%
    dplyr::select(repertoire_id, junction_aa, duplicate_frequency) %>%
    dplyr::mutate(
      Sequence = 1:n(),
      Sequence = as.factor(Sequence)
    ) %>%
    dplyr::rename(Frequency = duplicate_frequency) %>%
    dplyr::ungroup()

  subdominant <- dominant %>%
    dplyr::group_by(repertoire_id) %>%
    dplyr::summarize(
      Frequency = 1 - sum(Frequency),
      junction_aa = "All other sequences"
    ) %>%
    dplyr::select(repertoire_id, junction_aa, Frequency) %>%
    dplyr::mutate(Sequence = as.factor(11)) %>%
    dplyr::ungroup()
  topfreq <- bind_rows(dominant, subdominant) %>%
    dplyr::arrange(repertoire_id, Sequence, desc(Frequency)) %>%
    dplyr::mutate(Frequency = Frequency * 100)
  getPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))
  sample_order <- subdominant %>%
    dplyr::arrange(Frequency) %>%
    dplyr::select(repertoire_id) %>%
    dplyr::pull()
  ggplot2::ggplot(topfreq, ggplot2::aes_string(
    x = "repertoire_id",
    y = "Frequency",
    fill = "Sequence",
    label = "Frequency"
  ),
  text = "junction_aa"
  ) +
    ggplot2::geom_bar(stat = "identity") +
    ggplot2::scale_x_discrete(limits = sample_order) +
    ggplot2::scale_fill_manual(values = getPalette(top + 1)) +
    ggplot2::theme_classic() +
    ggplot2::scale_y_continuous(expand = c(0, 0)) +
    ggplot2::theme(legend.position = "none") +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
      axis.text.y = ggplot2::element_text(size = 10)
    ) +
    labs(x = "", y = "Frequency (%)")
}
