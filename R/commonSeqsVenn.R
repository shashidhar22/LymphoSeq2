#' Common sequences Venn diagram
#'
#' Creates a Venn diagram comparing the number of common sequences in two or
#' three repertoire_ids using native ggplot2 implementation.
#'
#' @param repertoire_ids A character vector of two or three names of
#'  repertoire_ids in [productiveSeq()] table to compare.
#' @param amino_table A tibble of amino acid sequences generated
#'  by the function [productiveSeq()].
#' @return Returns a ggplot2 Venn diagram of the number of common sequences between
#'  two or three repertoire_ids.
#' @seealso [LymphoSeq2::productiveSeq()], [LymphoSeq2::commonSeqs()],
#' [LymphoSeq2::commonSeqsPlot()], [LymphoSeq2::commonSeqsBar()]
#' @examples
#' file_path <- system.file("extdata", "TCRB_sequencing",
#'  package = "LymphoSeq2")
#' study_table <- LymphoSeq2::readImmunoSeq(path = file_path, threads = 1)
#' study_table <- LymphoSeq2::topSeqs(study_table, top = 100)
#' amino_table <- LymphoSeq2::productiveSeq(
#'   study_table = study_table,
#'   aggregate = "junction_aa"
#' )
#' # Plot a triple Venn diagram
#' LymphoSeq2::commonSeqsVenn(
#'   repertoire_ids = c(
#'     "TRB_Unsorted_0",
#'     "TRB_Unsorted_32", "TRB_Unsorted_83"
#'   ),
#'   amino_table = amino_table
#' )
#' # Plot a double Venn diagram
#' LymphoSeq2::commonSeqsVenn(repertoire_ids = c(
#'   "TRB_Unsorted_0",
#'   "TRB_Unsorted_32"
#' ), amino_table = amino_table)
#' @export
commonSeqsVenn <- function(repertoire_ids, amino_table) {
  if (base::length(repertoire_ids) > 3 | base::length(repertoire_ids) < 2) {
    stop("Please enter 2 or 3 repertoire_ids.")
  }

  # Extract sequences for each repertoire
  seqs_list <- lapply(repertoire_ids, function(id) {
    amino_table |>
      dplyr::filter(repertoire_id == id) |>
      dplyr::pull(junction_aa)
  })
  names(seqs_list) <- repertoire_ids

  if (base::length(repertoire_ids) == 2) {
    return(draw_venn2_native(seqs_list, repertoire_ids))
  } else if (base::length(repertoire_ids) == 3) {
    return(draw_venn3_native(seqs_list, repertoire_ids))
  }
}

#' Draw 2-way Venn diagram with ggplot2
#'
#' @param seqs_list Named list of sequence vectors
#' @param labels Character vector of set names
#' @return ggplot object
#' @keywords internal
draw_venn2_native <- function(seqs_list, labels) {
  set1 <- seqs_list[[1]]
  set2 <- seqs_list[[2]]

  # Calculate counts
  only1 <- length(setdiff(set1, set2))
  only2 <- length(setdiff(set2, set1))
  both <- length(intersect(set1, set2))

  # Circle parameters
  r <- 1  # radius
  d <- 0.6 * r  # distance between centers

  # Create circles
  theta <- seq(0, 2 * pi, length.out = 100)

  circle1 <- data.frame(
    x = -d/2 + r * cos(theta),
    y = r * sin(theta),
    set = labels[1]
  )

  circle2 <- data.frame(
    x = d/2 + r * cos(theta),
    y = r * sin(theta),
    set = labels[2]
  )

  # Label positions
  label_data <- data.frame(
    x = c(-d/2 - r/2, d/2 + r/2, 0),
    y = c(0, 0, 0),
    label = c(only1, only2, both),
    type = c("only1", "only2", "both")
  )

  # Set labels
  set_labels <- data.frame(
    x = c(-d/2 - r, d/2 + r),
    y = c(r + 0.3, r + 0.3),
    label = labels
  )

  # Create plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = circle1, ggplot2::aes(x = x, y = y),
                         fill = "#3288bd", alpha = 0.5, color = "#3288bd", linewidth = 2) +
    ggplot2::geom_polygon(data = circle2, ggplot2::aes(x = x, y = y),
                         fill = "#d53e4f", alpha = 0.5, color = "#d53e4f", linewidth = 2) +
    ggplot2::geom_text(data = label_data, ggplot2::aes(x = x, y = y, label = label),
                      size = 5, fontface = "bold") +
    ggplot2::geom_text(data = set_labels, ggplot2::aes(x = x, y = y, label = label),
                      size = 4, fontface = "bold") +
    ggplot2::coord_fixed() +
    ggplot2::theme_void()

  return(p)
}

#' Draw 3-way Venn diagram with ggplot2
#'
#' @param seqs_list Named list of sequence vectors
#' @param labels Character vector of set names
#' @return ggplot object
#' @keywords internal
draw_venn3_native <- function(seqs_list, labels) {
  set1 <- seqs_list[[1]]
  set2 <- seqs_list[[2]]
  set3 <- seqs_list[[3]]

  # Calculate all intersections
  only1 <- length(setdiff(setdiff(set1, set2), set3))
  only2 <- length(setdiff(setdiff(set2, set1), set3))
  only3 <- length(setdiff(setdiff(set3, set1), set2))
  n12_not3 <- length(setdiff(intersect(set1, set2), set3))
  n13_not2 <- length(setdiff(intersect(set1, set3), set2))
  n23_not1 <- length(setdiff(intersect(set2, set3), set1))
  n123 <- length(Reduce(intersect, list(set1, set2, set3)))

  # Circle parameters
  r <- 1
  angle_offset <- pi / 6

  # Create three circles
  theta <- seq(0, 2 * pi, length.out = 100)

  circle1 <- data.frame(
    x = r * cos(angle_offset) + r * cos(theta),
    y = r * sin(angle_offset) + r * sin(theta),
    set = labels[1]
  )

  circle2 <- data.frame(
    x = r * cos(pi - angle_offset) + r * cos(theta),
    y = r * sin(pi - angle_offset) + r * sin(theta),
    set = labels[2]
  )

  circle3 <- data.frame(
    x = r * cos(3 * pi / 2) + r * cos(theta),
    y = r * sin(3 * pi / 2) + r * sin(theta),
    set = labels[3]
  )

  # Label positions (approximate)
  label_data <- data.frame(
    x = c(1.2, -1.2, 0, 0.5, -0.5, 0, 0),
    y = c(0.8, 0.8, -1.3, 0.2, 0.2, -0.3, 0),
    label = c(only1, only2, only3, n12_not3, n13_not2, n23_not1, n123)
  )

  # Set labels
  set_labels <- data.frame(
    x = c(1.8, -1.8, 0),
    y = c(1.2, 1.2, -1.8),
    label = labels
  )

  # Create plot
  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(data = circle1, ggplot2::aes(x = x, y = y),
                         fill = "#3288bd", alpha = 0.4, color = "#3288bd", linewidth = 2) +
    ggplot2::geom_polygon(data = circle2, ggplot2::aes(x = x, y = y),
                         fill = "#abdda4", alpha = 0.4, color = "#abdda4", linewidth = 2) +
    ggplot2::geom_polygon(data = circle3, ggplot2::aes(x = x, y = y),
                         fill = "#d53e4f", alpha = 0.4, color = "#d53e4f", linewidth = 2) +
    ggplot2::geom_text(data = label_data, ggplot2::aes(x = x, y = y, label = label),
                      size = 5, fontface = "bold") +
    ggplot2::geom_text(data = set_labels, ggplot2::aes(x = x, y = y, label = label),
                      size = 4, fontface = "bold") +
    ggplot2::coord_fixed() +
    ggplot2::theme_void()

  return(p)
}
