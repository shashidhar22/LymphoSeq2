#' Chord diagram of VJ or DJ gene associations
#'
#' Creates a chord diagram showing VJ or DJ gene associations from one or more
#' samples.
#'
#' @param study_table A tibble consisting of frequencies of antigen receptor
#' sequences.  "v_family", "j_family", and if applicable, "d_family"
#' are required columns.  Using output from the LymphoSeq2 function topSeqs is
#' recommended.
#' @param association A character vector of gene families to associate. Options
#' include "VJ" or "DJ".
#' @param colors A character vector of 2 colors corresponding to the V/D and J
#' gene colors respectively.
#' @details The size of the ribbons connecting VJ or DJ genes correspond to the
#' number of samples or number of sequences that make up that recombination
#' event. The thicker the ribbon, the higher the frequency of the recombination.
#' @return Returns a chord diagram showing VJ or DJ gene associations from one or
#' more samples.
#' @seealso \code{\link{topSeqs}}
#' @examples
#' file_path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq2")
#'
#' stable <- readImmunoSeq(file_path)
#'
#' ntable <- productiveSeq(stable, aggregate = "junction")
#'
#' top_seqs <- topSeqs(ntable, top = 1)
#'
#' chordDiagramVDJ(top_seqs, association = "VJ", colors = c("red", "blue"))
#'
#'
#' @export
#' @importFrom circlize colorRamp2 chordDiagram
#' @import magrittr
chordDiagramVDJ <- function(study_table, association = "VJ", colors = c("red", "blue")) {
    if (association == "VJ") {
        if (!base::all(c("v_family", "j_family") %in% base::colnames(study_table))) {
            stop("The source data frame does not contain the required columns 'v_family' and 'j_family'.")
        }
        vj <- study_table %>%
            dplyr::select(v_family, j_family) %>%
            dplyr::mutate(v_family = tidyr::replace_na(v_family, "unresolved"),
                j_family = tidyr::replace_na(j_family, "unresolved"))
        vj <- vj %>%
            dplyr::group_by(v_family, j_family) %>%
            dplyr::summarize(duplicate_count = dplyr::n()) %>%
            dplyr::mutate(duplicate_count = tidyr::replace_na(duplicate_count, 0)) %>%
            dplyr::rename(from = v_family,
                to = j_family,
                value = duplicate_count) %>%
            dplyr::ungroup() %>%
            dplyr::mutate(value = as.integer(value))
        vcolors <- base::rep(colors[1], base::length(base::unique(vj$from)))
        jcolors <- base::rep(colors[2], base::length(base::unique(vj$to)))
        ribbon.color <- circlize::colorRamp2(range(vj$value), c("grey", "black"))
        circlize::chordDiagram(vj,
            annotationTrack = c("grid", "name"),
            grid.col = c(vcolors, jcolors),
            col = ribbon.color)
    }
    if (association == "DJ") {
        if (!base::all(c("d_family", "j_family") %in% base::colnames(study_table))) {
            stop("The source data frame does not contain the required columns 'd_family' and 'j_family'.")
        }
        dj <- study_table %>%
            dplyr::select(d_family, j_family) %>%
            dplyr::mutate(d_family = tidyr::replace_na(d_family, "Unresolved"),
                j_family = tidyr::replace_na(j_family, "Unresolved"))
        dj <- dj %>%
            dplyr::group_by(d_family, j_family) %>%
            dplyr::summarize(duplicate_count = dplyr::n()) %>%
            dplyr::mutate(duplicate_count = tidyr::replace_na(duplicate_count, 0)) %>%
            dplyr::mutate(duplicate_count = tidyr::replace_na(duplicate_count, 0)) %>%
            dplyr::rename(from = d_family,
                to = j_family,
                value = duplicate_count)
        dcolors <- base::rep(colors[1], base::length(base::unique(dj$from)))
        jcolors <- base::rep(colors[2], base::length(base::unique(dj$to)))
        ribbon.color <- circlize::colorRamp2(range(dj$value), c("grey", "black"))
        circlize::chordDiagram(dj,
            annotationTrack = c("grid", "name"),
            grid.col = c(dcolors, jcolors),
            col = ribbon.color)
    }
}
