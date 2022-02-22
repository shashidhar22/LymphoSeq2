#' Common sequences Venn diagram
#' 
#' Creates a Venn diagram comparing the number of common sequences in two or 
#' three repertoire_ids.
#' 
#' @param repertoire_ids A character vector of two or three names of repertoire_ids in 
#' productiveSeq table to compare.
#' @param productive_aa A tibble of amino acid sequences generated
#' by the LymphoSeq function productiveSeq.
#' @return Returns a a Venn diagram of the number of common sequences between
#' two or three repertoire_ids.
#' @seealso \code{\link{commonSeqs}}
#' @examples
#' file_path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq2")
#' 
#' stable <- readImmunoSeq(path = file_path)
#' 
#' atable <- productiveSeq(study_table = stable, aggregate = "junction_aa")
#' 
#' # Plot a triple Venn diagram
#' commonSeqsVenn(repertoire_ids = c("TRB_Unsorted_0", 
#'    "TRB_Unsorted_32", "TRB_Unsorted_83"), 
#'    productive_aa = atable)
#' 
#' # Plot a double Venn diagram
#' commonSeqsVenn(repertoire_ids = c("TRB_Unsorted_0", 
#'    "TRB_Unsorted_32"), productive_aa = atable)
#' 
#' @export
commonSeqsVenn <- function(repertoire_ids, productive_aa) {
    if (base::length(repertoire_ids) > 3 | base::length(repertoire_ids) < 2) {
        stop("Please enter 2 or 3 repertoire_ids.")
    }
    if (base::length(repertoire_ids) == 2) {
        a <- productive_aa %>% 
             dplyr::filter(repertoire_id == repertoire_ids[[1]])
        b <- productive_aa %>% 
             dplyr::filter(repertoire_id == repertoire_ids[[2]])
        grid::grid.newpage()
        venn <- VennDiagram::draw.pairwise.venn(area1 = length(a$junction_aa), 
                                                area2 = length(b$junction_aa), 
                                                cross.area = length(intersect(a$junction_aa, 
                                                                              b$junction_aa)), 
                                                category = c(repertoire_ids[1], 
                                                             repertoire_ids[2]), 
                                                cat.fontfamily = rep("sans", 2), 
                                                fontfamily = rep("sans", 3), 
                                                fill = c("#3288bd", "#d53e4f"), 
                                                cat.pos = c(0, 0),
                                                cat.dist = rep(0.025, 2),
                                                cex = 1, 
                                                cat.cex = 0.7,
                                                lwd = rep(2, 2))
        grid::grid.draw(venn)
    }
    if (base::length(repertoire_ids) == 3) {
        a <- productive_aa %>% 
             dplyr::filter(repertoire_id == repertoire_ids[[1]])
        b <- productive_aa %>% 
             dplyr::filter(repertoire_id == repertoire_ids[[2]])
        c <- productive_aa %>% 
             dplyr::filter(repertoire_id == repertoire_ids[[3]])
        grid::grid.newpage()
        venn <- VennDiagram::draw.triple.venn(area1 = length(a$junction_aa), 
                                              area2 = length(b$junction_aa), 
                                              area3 = length(c$junction_aa), 
                                              n12 = length(intersect(a$junction_aa, 
                                                                     b$junction_aa)), 
                                              n23 = length(intersect(b$junction_aa, 
                                                                     c$junction_aa)), 
                                              n13 = length(intersect(a$junction_aa, 
                                                                     c$junction_aa)), 
                                              n123 = length(Reduce(intersect, 
                                                                   list(a$junction_aa, 
                                                                        b$junction_aa, 
                                                                        c$junction_aa))), 
                                              category = c(repertoire_ids[1], 
                                                           repertoire_ids[2], 
                                                           repertoire_ids[3]), 
                                              cat.fontfamily = rep("sans", 3), 
                                              fontfamily = rep("sans", 7), 
                                              fill = c("#3288bd", "#abdda4", "#d53e4f"), 
                                              cat.pos = c(0, 0, 180), 
                                              cat.dist = rep(0.025, 3),
                                              cex = 1, 
                                              cat.cex = 0.7,
                                              lwd = rep(2, 3))
        grid::grid.draw(venn)
    }
} 
