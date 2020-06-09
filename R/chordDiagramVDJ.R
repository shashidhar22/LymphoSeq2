#' Chord diagram of VJ or DJ gene associations
#' 
#' Creates a chord diagram showing VJ or DJ gene associations from one or more 
#' samples.
#' 
#' @param repertoire_id A tibble consisting of frequencies of antigen receptor 
#' sequences.  "v_family", "j_family", and if applicable, "d_family" 
#' are a required columns.  Using output from the LymphoSeq function topSeqs is
#' recommended.
#' @param association A character vector of gene familes to associate.  Options 
#' include "VJ" or "DJ".
#' @param colors A character vector of 2 colors corresponding to the V/D and J 
#' gene colors respectively.
#' @details The size of the ribbons connecting VJ or DJ genes correspond to the 
#' number of samples or number of sequences that make up that recombination 
#' event.  The thicker the ribbon, the higher the frequency of the recombination.
#' @return Returns a chord diagram showing VJ or DJ gene associations from one or 
#' more samples.
#' @seealso \code{\link{topSeqs}}
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' study_table <- readImmunoSeq(path = file.path)
#' 
#' productive_nt <- productiveSeq(study_table = study_table, aggregate = "junction")
#' 
#' top_seqs <- topSeqs(productive_table = productive_nt, top = 1)
#' 
#' chordDiagramVDJ(repertoire_id = top_seqs, association = "VJ", colors = c("red", "blue"))
#' 
#'
#' @export
#' @importFrom circlize colorRamp2 chordDiagram
#' @import tidyverse
chordDiagramVDJ <- function(study_table, association = "VJ", colors = c("red", "blue")) {
    if (association == "VJ") {
        if (!all(c("v_family", "j_family") %in% colnames(study_table))) {
            stop("The source data frame does not contain the required columns 'v_family' and 'j_family'.")
        }
        vj <- study_table %>% 
              select(v_family, j_family) %>% 
              mutate(v_family = replace_na(v_family, "Unresolved"), j_family = replace_na(j_family, "Unresolved"))
        vj <- vj %>% 
              group_by(v_family, j_family) %>% 
              summarize(duplicate_count = n()) %>% 
              pivot_wider(id_cols=v_family, names_from = j_family, values_from = duplicate_count) 
        row_names <- vj$v_family 
        vj <- vj %>% select(-v_family) %>% as.matrix()
        rownames(vj) <- row_names
        ribbon.color <- circlize::colorRamp2(range(vj), c("grey", "black"))
        circlize::chordDiagram(vj, 
                               annotationTrack = c("grid", "name"), 
                               grid.col = c(rep(colors[1], dim(vj)[1]), rep(colors[2], dim(vj)[2])), 
                               col = ribbon.color)
    }
    if (association == "DJ") {
        if (!all(c("d_family", "j_family") %in% colnames(study_table))) {
            stop("The source data frame does not contain the required columns 'd_family' and 'j_family'.")
        }
        dj <- study_table %>% 
              select(d_family, j_family) %>% 
              mutate(d_family = replace_na(d_family, "Unresolved"), j_family = replace_na(j_family, "Unresolved"))
        dj <- dj %>% 
              group_by(d_family, j_family) %>% 
              summarize(duplicate_count = n()) %>% 
              pivot_wider(id_cols=d_family, names_from = j_family, values_from = duplicate_count)
        row_names <- dj$d_family 
        dj <- dj %>% select(-d_family) %>% as.matrix()
        rownames(dj) <- row_names
        ribbon.color <- circlize::colorRamp2(range(dj), c("grey", "black"))
        circlize::chordDiagram(dj, 
                               annotationTrack = c("grid", "name"), 
                               grid.col = c(rep(colors[1], dim(dj)[1]), rep(colors[2], dim(dj)[2])), 
                               col = ribbon.color)
    }
}