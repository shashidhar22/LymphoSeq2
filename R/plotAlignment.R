#' Visualize multiple sequence alignment of CDR3 sequences
#'
#' Generate MSA alignment figures from the results of [alignSeq()]
#'
#' @param msa An msa object obtained from [alignSeq()] function in LymphoSeq2.
#' @return Multiple sequence alignment plot.
#' @seealso The function utilizes ggmsa package for visualizations. Further
#' details on ggmsa can be found at the link mentioned below.
#' \url{https://cran.r-project.org/web/packages/ggmsa/vignettes/ggmsa.html}
#' @examples
#' library(ggmsa)
#' file_path <- system.file("extdata", "IGH_sequencing", package = "LymphoSeq2")
#' study_table <- LymphoSeq2::readImmunoSeq(path = file_path, threads = 1)
#' study_table <- LymphoSeq2::topSeqs(study_table, top = 100)
#' nucleotide_table <- LymphoSeq2::productiveSeq(study_table,
#'  aggregate = "junction")
#' msa <- LymphoSeq2::alignSeq(nucleotide_table,
#'   repertoire_id = "IGH_MVQ92552A_BL",
#'   type = "junction_aa", method = "ClustalW"
#' )
#' LymphoSeq2::plotAlignment(msa)
#' @export
plotAlignment <- function(msa) {
  if (class(msa)[1] == "MsaDNAMultipleAlignment") {
    msa <- Biostrings::DNAMultipleAlignment(msa)
    names(msa@unmasked) <- paste(names(msa@unmasked),
      seq(1:length(names(msa@unmasked))),
      sep = "_"
    )
    ggmsa::ggmsa(msa, font = NULL, color = "Chemistry_NT") +
      ggmsa::geom_msaBar()
  } else {
    msa <- Biostrings::AAMultipleAlignment(msa)
    names(msa@unmasked) <- paste(names(msa@unmasked),
      seq(1:length(names(msa@unmasked))),
      sep = "_"
    )
    ggmsa::ggmsa(msa, font = NULL, color = "Chemistry_AA") +
    ggmsa::geom_msaBar()
  }
}
