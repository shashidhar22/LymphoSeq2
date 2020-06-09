#' Pairwise comparison plot
#' 
#' Creates a heat map from a similarity or Bhattacharyya matrix.
#' 
#' @param matrix A similarity or Bhattacharyya matrix produced by the LymphoSeq2 
#' scoringMatrix function.
#' @return A pairwise comparison heat map.
#' @details The plot is made using the package ggplot2 and can be reformatted
#' using ggplot2 functions.  See examples below.
#' @seealso An excellent resource for examples on how to reformat a ggplot can 
#' be found in the R Graphics Cookbook online (\url{http://www.cookbook-r.com/Graphs/}).
#' The functions to create the similarity or Bhattacharyya matrix can be found 
#' here: \code{\link{similarityMatrix}} and \code{\link{bhattacharyyaMatrix}}
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' study_table <- readImmunoSeq(path = file.path)
#' 
#' productive_aa <- productiveSeq(study_table = study_table, aggregate = "junction_aa")
#' 
#' similarity_matrix <- scoringMatrix(productive_table = productive_aa, mode="Similarity")
#' 
#' pairwisePlot(matrix = similarity_matrix)
#' 
#' bhattacharyya_matrix <- scoringMatrix(productive_table = productive_aa, mode="Bhattacharyya")
#' 
#' pairwisePlot(matrix = bhattacharyya_matrix)
#' 
#' # Change plot color, title legend, and add title
#' pairwisePlot(matrix = similarity.matrix) + 
#'    ggplot2::scale_fill_gradient(low = "#deebf7", high = "#3182bd") + 
#'    ggplot2::labs(fill = "Similarity score") + ggplot2::ggtitle("Figure Title")
#' @export
#' @import ggplot2
#' @import tidyverse
#' @importFrom reshape melt.data.frame
#' @importFrom stats na.omit
pairwisePlot <- function(matrix) {
    i <- 1
    l <- length(matrix)
    p <- l - 1
    for (i in 1:p) {
        j <- i + 1
        matrix[i, j:l] <- NA
    }
    matrix$repertoire_ids <- rownames(matrix)
    melt <- reshape::melt.data.frame(matrix, id.vars = "repertoire_ids")
    melt <- na.omit(melt)
    melt$variable <- factor(melt$variable, levels = rownames(matrix))
    melt$repertoire_ids <- factor(melt$repertoire_ids, levels = rev(rownames(matrix)))
    names(melt) = c("repertoire_id.x", "repertoire_id.y", "Score")
    ggplot(data = melt, aes_string(x = "repertoire_id.x", y = "repertoire_id.y", fill = "Score")) + 
    ggplot2::geom_tile() + 
    ggplot2::scale_fill_gradient(low = "#fee8c8", high = "#e34a33") + 
    ggplot2::theme_classic() + 
    ggplot2::labs(x = "", y = "", fill = "") + 
    ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
}