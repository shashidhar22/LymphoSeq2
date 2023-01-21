#' Visualize multiple sequence alignment of CDR3 sequences
#'
#' Generate MSA alignment figures from the results of `alignSeq()`
#' @param msa An msa object obtained from `alignSeq()` function in LymphoSeq2.
#' @return Multiple sequence alignment plot.
#' @seealso The function utilizes ggmsa package for visualizations. Further
#' details on ggmsa can be found at the link mentioned below.
#' \url{https://cran.r-project.org/web/packages/ggmsa/vignettes/ggmsa.html}
#' @examples
#' file_path <- system.file("extdata", "IGH_sequencing", package = "LymphoSeq2")
#' study_table <- LymphoSeq2::readImmunoSeq(path = file_path)
#' nucleotide_table <- LymphoSeq2::productiveSeq(study_table, aggregate = "junction")
#' msa <- LymphoSeq2::alignSeq(nucleotide_table,
#'   repertoire_id = "IGH_MVQ92552A_BL",
#'   type = "junction_aa", method = "ClustalW"
#' )
#' LymphoSeq2::plotAlignment(msa)
#' @import ggmsa
#' @export
plotAlignment <- function(msa) {
  if (class(msa)[1] == "MsaDNAMultipleAlignment") {
    msa <- Biostrings::DNAMultipleAlignment(msa)
    ggmsa::ggmsa(msa, font = NULL, color = "Chemistry_NT", seq_name = FALSE) + 
      ggmsa::geom_msaBar()  
  } else if (class(msa)[1] == "MsaAAMultipleAlignment") {
    msa <- Biostrings::AAMultipleAlignment(msa)
    ggmsa::ggmsa(msa, font = NULL, color = "Chemistry_AA", seq_name = FALSE) + 
      ggmsa::geom_msaBar()
  }
}
