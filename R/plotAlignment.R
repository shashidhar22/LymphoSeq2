#' Visualize multiple sequence alignment of CDR3 sequences
#' 
#' Generate MSA alignment figures from the results of alignSeq
#' @param msa An msa object obtained from alignSeq function in LymphoSeq2.
#' @return Multiple sequence alignment plot.
#' @seealso The function utilizes ggmsa package for visualizations. Further details 
#' on ggmsa can be found at the link mentioned below.
#' \url{https://cran.r-project.org/web/packages/ggmsa/vignettes/ggmsa.html}
#' @examples
#' file_path <- system.file("extdata", "IGH_sequencing", package = "LymphoSeq")
#' 
#' file_list <- readImmunoSeq(path = file_path)
#' 
#' productive_nt <- productiveSeq(file_list = file_list, aggregate = "junction")
#' 
#' msa <- alignSeq(list = productive_nt, repertoire_id = "IGH_MVQ92552A_BL", type = "junction", 
#'          method = "ClustalW", output = "console")
#' plotAlignment(msa)
#' @export
plotAlignment <- function(msa) {
  ggmsa::ggmsa(msa, font = NULL, color = "Chemistry_NT")
}
