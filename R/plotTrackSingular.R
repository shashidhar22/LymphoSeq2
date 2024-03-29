#' Highlight clones
#'
#' Create alluvial plots highlighting each sequence in amino acid list
#'
#' @param clone_table A tibble of productive amino acid sequences to highlight
#' generated by LymphoSeq2 function [cloneTrack()]
#' @examples
#' file_path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq2")
#' study_table <- LymphoSeq2::readImmunoSeq(path = file_path, threads = 1)
#' study_table <- LymphoSeq2::topSeqs(study_table, top = 100)
#' amino_table <- LymphoSeq2::productiveSeq(study_table, aggregate = "junction_aa")
#' top_seq <- LymphoSeq2::topSeqs(amino_table, top = 3)
#' clone_table <- LymphoSeq2::cloneTrack(
#'   study_table = top_seq,
#'   sample_list = c("TRB_CD8_949", "TRB_CD8_CMV_369")
#' )
#' LymphoSeq2::plotTrackSingular(clone_table)
#' @return A list of alluvial plots highlighting a single sequence
#' @details The plot is made using the package ggplot2 and can be reformatted
#' using ggplot2 functions.
#' @export
plotTrackSingular <- function(clone_table) {
  alist <- clone_table |>
    dplyr::pull(junction_aa) |>
    unique()
  plots <- alist |> purrr::map(~ highlightPlot(.x, clone_table))
  return(plots)
}
#' Highlight specific sequences
#'
#' @describeIn plotTrackSingular Highlight a specific amino acid sequence
#' @param aseq CDR3 amino acid sequence to highlight
#' @param clone_table  tibble of productive amino acid sequences to highlight
#' generated by LymphoSeq function cloneTrack
highlightPlot <- function(aseq, clone_table) {
  pal_table <- clone_table |>
    dplyr::select(junction_aa) |>
    dplyr::distinct() |>
    dplyr::mutate(color = dplyr::if_else(junction_aa == aseq, "#7570b3", "grey"))
  cpal <- pal_table |> dplyr::pull(color)
  names(cpal) <- pal_table |> dplyr::pull(junction_aa)
  tplot <- plotTrack(clone_table) +
    ggplot2::scale_fill_manual(values = cpal, breaks = aseq, name = "CDR3 Sequence") +
    ggplot2::theme(legend.position = "bottom")
  return(tplot)
}
