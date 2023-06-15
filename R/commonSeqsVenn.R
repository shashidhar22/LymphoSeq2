#' Common sequences Venn diagram
#'
#' Creates a Venn diagram comparing the number of common sequences in two or
#' three repertoire_ids.
#'
#' @param repertoire_ids A character vector of two or three names of repertoire_ids in
#' [productiveSeq()] table to compare.
#' @param amino_table A tibble of amino acid sequences generated
#' by the function [productiveSeq()].
#' @return Returns a a Venn diagram of the number of common sequences between
#' two or three repertoire_ids.
#' @seealso [LymphoSeq2::productiveSeq()], [LymphoSeq2::commonSeqs()],
#' [LymphoSeq2::commonSeqsPlot()], [LymphoSeq2::commonSeqsBar()]
#' @examples
#' file_path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq2")
#' study_table <- LymphoSeq2::readImmunoSeq(path = file_path)
#' study_table <- LymphoSeq2::topSeqs(study_table, top = 100)
#' amino_table <- LymphoSeq2::productiveSeq(study_table = study_table,
#'   aggregate = "junction_aa")
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
#' @import magrittr
commonSeqsVenn <- function(repertoire_ids, amino_table) {
  if (base::length(repertoire_ids) > 3 | base::length(repertoire_ids) < 2) {
    stop("Please enter 2 or 3 repertoire_ids.")
  }
  if (base::length(repertoire_ids) == 2) {
    a <- amino_table %>%
      dplyr::filter(repertoire_id == repertoire_ids[[1]])
    b <- amino_table %>%
      dplyr::filter(repertoire_id == repertoire_ids[[2]])
    grid::grid.newpage()
    venn <- VennDiagram::draw.pairwise.venn(
      area1 = length(a$junction_aa),
      area2 = length(b$junction_aa),
      cross.area = length(intersect(a$junction_aa, b$junction_aa)),
      category = c(repertoire_ids[1], repertoire_ids[2]),
      cat.fontfamily = rep("sans", 2), fontfamily = rep("sans", 3),
      fill = c("#3288bd", "#d53e4f"), cat.pos = c(0, 0),
      cat.dist = rep(0.025, 2), cex = 1, cat.cex = 0.7, lwd = rep(2, 2)
    )
    grid::grid.draw(venn)
  }
  if (base::length(repertoire_ids) == 3) {
    a <- amino_table %>%
      dplyr::filter(repertoire_id == repertoire_ids[[1]])
    b <- amino_table %>%
      dplyr::filter(repertoire_id == repertoire_ids[[2]])
    c <- amino_table %>%
      dplyr::filter(repertoire_id == repertoire_ids[[3]])
    grid::grid.newpage()
    venn <- VennDiagram::draw.triple.venn(
      area1 = length(a$junction_aa),
      area2 = length(b$junction_aa), area3 = length(c$junction_aa),
      n12 = length(intersect(a$junction_aa, b$junction_aa)),
      n23 = length(intersect(b$junction_aa, c$junction_aa)),
      n13 = length(intersect(a$junction_aa, c$junction_aa)),
      n123 = length(Reduce(intersect, list(
        a$junction_aa, b$junction_aa,
        c$junction_aa
      ))),
      category = c(
        repertoire_ids[1], repertoire_ids[2],
        repertoire_ids[3]
      ),
      cat.fontfamily = rep("sans", 3), fontfamily = rep("sans", 7),
      fill = c("#3288bd", "#abdda4", "#d53e4f"),
      cat.pos = c(0, 0, 180), cat.dist = rep(0.025, 3), cex = 1,
      cat.cex = 0.7, lwd = rep(2, 3)
    )
    grid::grid.draw(venn)
  }
}
