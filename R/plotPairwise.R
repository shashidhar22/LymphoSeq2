#' Pairwise comparison plot
#'
#' Creates a heat map from a Bhattacharyya, Similarity, Sorensen, or PSI matrix.
#'
#' @param matrix A Bhattacharyya, Similarity, Sorensen, or PSI matrix produced
#' by the LymphoSeq2 `scoringMatrix()` function.
#' @return A pairwise comparison heat map.
#' @details The plot is made using the package ggplot2 and can be reformatted
#' using ggplot2 functions.  See examples below.
#' @seealso An excellent resource for examples on how to reformat a ggplot can
#' be found in the R Graphics Cookbook online (\url{http://www.cookbook-r.com/Graphs/}).
#' The functions to create the similarity matrices can be found
#' here: [LymphoSeq2::scoringMatrix()]
#' @examples
#' file_path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq2")
#' study_table <- LymphoSeq2::readImmunoSeq(path = file_path)
#' amino_table <- LymphoSeq2::productiveSeq(study_table, aggregate = "junction_aa")
#' # Plot similarity using Similarity score
#' similarity_matrix <- LymphoSeq2::scoringMatrix(
#'   productive_table = amino_table,
#'   mode = "Similarity"
#' )
#' LymphoSeq2::pairwisePlot(matrix = similarity_matrix)
#' # Plot similarity using Bhattacharyya score
#' bhattacharyya_matrix <- LymphoSeq2::scoringMatrix(
#'   productive_table = amino_table,
#'   mode = "Bhattacharyya"
#' )
#' LymphoSeq2::pairwisePlot(matrix = bhattacharyya_matrix)
#' # Change plot color, title legend, and add title
#' LymphoSeq2::pairwisePlot(matrix = similarity_matrix) +
#'   ggplot2::scale_fill_gradient(low = "#deebf7", high = "#3182bd") +
#'   ggplot2::labs(fill = "Similarity score") +
#'   ggplot2::ggtitle("Pairwise similarity score")
#' @export
#' @import magrittr
pairwisePlot <- function(matrix) {
  samples <- rownames(matrix)
  matrix[base::lower.tri(matrix)] <- NA
  matrix <- matrix %>%
    tibble::as_tibble() %>%
    dplyr::mutate(repertoire_id = samples) %>%
    dplyr::select(repertoire_id, dplyr::everything()) %>%
    tidyr::pivot_longer(-repertoire_id,
      names_to = "repertoire_id_y",
      values_to = "score", values_drop_na = TRUE
    ) %>%
    dplyr::arrange(repertoire_id, repertoire_id_y)
  ggplot2::ggplot(data = matrix, ggplot2::aes_string(
    x = "repertoire_id", y = "repertoire_id_y",
    fill = "score"
  )) +
    ggplot2::geom_tile() +
    ggplot2::geom_text(
      ggplot2::aes(repertoire_id, repertoire_id_y,
        label = sprintf("%0.2f", round(score, digits = 2))
      ),
      color = "black", size = 4
    ) +
    ggplot2::scale_fill_gradient(low = "#fee8c8", high = "#e34a33") +
    ggplot2::theme_classic() +
    ggplot2::labs(x = "", y = "", fill = "") +
    ggplot2::scale_y_discrete(
      limits =
        rev(levels(as.factor(matrix$repertoire_id_y)))
    ) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(
      angle = 90, hjust = 1,
      vjust = 0.5
    ))
}
