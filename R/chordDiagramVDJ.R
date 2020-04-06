#' Chord diagram of VJ or DJ gene associations
#' 
#' Creates a chord diagram showing VJ or DJ gene associations from one or more 
#' samples.
#' 
#' @param sample A tibble consisting of frequencies of antigen receptor 
#' sequences.  "vFamilyName", "jFamilyName", and if applicable, "dFamilyName" 
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
#' productive_nt <- productiveSeq(study_table = study_table, aggregate = "nucleotide")
#' 
#' top_seqs <- topSeqs(productive_table = productive_nt, top = 1)
#' 
#' chordDiagramVDJ(sample = top_seqs, association = "VJ", colors = c("red", "blue"))
#' 
#'
#' @export
#' @importFrom circlize colorRamp2 chordDiagram
#' @import tidyverse
chordDiagramVDJ <- function(study_table, association = "VJ", colors = c("red", "blue")) {
    if (association == "VJ") {
        if (!all(c("vFamilyName", "jFamilyName") %in% colnames(study_table))) {
            stop("The source data frame does not contain the required columns 'vFamilyName' and 'jFamilyName'.")
        }
        vj <- study_table %>% 
              select(vFamilyName, jFamilyName) %>% 
              mutate(vFamilyName = replace_na(vFamilyName, "Unresolved"), jFamilyName = replace_na(jFamilyName, "Unresolved"))
        vj <- vj %>% 
              group_by(vFamilyName, jFamilyName) %>% 
              summarize(count = n()) %>% 
              pivot_wider(id_cols=vFamilyName, names_from = jFamilyName, values_from = count) 
        row_names <- vj$vFamilyName 
        vj <- vj %>% select(-vFamilyName) %>% as.matrix()
        rownames(vj) <- row_names
        ribbon.color <- circlize::colorRamp2(range(vj), c("grey", "black"))
        circlize::chordDiagram(vj, 
                               annotationTrack = c("grid", "name"), 
                               grid.col = c(rep(colors[1], dim(vj)[1]), rep(colors[2], dim(vj)[2])), 
                               col = ribbon.color)
    }
    if (association == "DJ") {
        if (!all(c("dFamilyName", "jFamilyName") %in% colnames(study_table))) {
            stop("The source data frame does not contain the required columns 'dFamilyName' and 'jFamilyName'.")
        }
        dj <- study_table %>% 
              select(dFamilyName, jFamilyName) %>% 
              mutate(dFamilyName = replace_na(dFamilyName, "Unresolved"), jFamilyName = replace_na(jFamilyName, "Unresolved"))
        dj <- dj %>% 
              group_by(dFamilyName, jFamilyName) %>% 
              summarize(count = n()) %>% 
              pivot_wider(id_cols=dFamilyName, names_from = jFamilyName, values_from = count)
        row_names <- dj$dFamilyName 
        dj <- dj %>% select(-dFamilyName) %>% as.matrix()
        rownames(dj) <- row_names
        ribbon.color <- circlize::colorRamp2(range(dj), c("grey", "black"))
        circlize::chordDiagram(dj, 
                               annotationTrack = c("grid", "name"), 
                               grid.col = c(rep(colors[1], dim(dj)[1]), rep(colors[2], dim(dj)[2])), 
                               col = ribbon.color)
    }
}