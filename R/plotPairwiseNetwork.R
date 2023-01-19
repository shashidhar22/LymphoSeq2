#' Pairwise comparison network graph
#' 
#' Creates a network graph from a Bhattacharyya, Similarity, Sorensen, or 
#' PSI matrix.
#;
#' @param matrix A Bhattacharyya, Similarity, Sorensen, or PSI matrix produced
#' by the LymphoSeq2 `scoringMatrix()` function.
#' @return A network graph visualizing pairwise comparisons. The thicker the
#' line connecting two nodes, the greater the similarity.
#' @examples
#' 
#' file_path <- system.file("extdata", "TCRB_sequencing",
#'   package = "LymphoSeq2")
#' study_table <- LymphoSeq2::readImmunoSeq(path = file_path)
#' amino_table <- LymphoSeq2::productiveSeq(study_table = study_table,
#'   aggregate = "junction_aa")
#' matrix <- LymphoSeq2::scoringMatrix(amino_table, mode = "Similarity")
#' network_graph <- LymphoSeq2::pairwiseNetwork(matrix)
#' @export
pairwiseNetwork <- function(matrix) {
    g <- igraph::graph_from_adjacency_matrix(matrix, mode = "undirected",
                                             weighted = TRUE, diag = FALSE)

    plot(g, layout = igraph::layout.circle,
            edge.width = igraph::E(g)$weight * 5,
            vertex.size = 10,
            vertex.label.cex = 1,
            vertex.label.dist = 2)
}