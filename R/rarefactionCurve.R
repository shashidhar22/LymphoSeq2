#' Plot rarefication and extrapolation curves for samples
#' 
#' Given a study table, for each sample plot rarefaction curves to estimate 
#' repertoire diversity. The method used to generate the rarefaction curve
#' is derived from Chao et al., (2014) using the iNEXT library
#'
#' @param study_table A tibble consisting antigen receptor sequencing
#' data imported by the LymphoSeq2 function readImmunoSeq. "junction_aa",
#' "duplicate_count", and "duplicate_frequency" are required columns.
#' @seealso \code{\link{runINext}}
#' @examples
#' file_path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq2")
#' stable <- readImmunoSeq(path = file_path)
#' plotRarefactionCurve(stable)
#'
#' @export
plotRarefactionCurve <- function(study_table) {
    rarefaction_tables <- study_table %>% 
                          dplyr::group_by(repertoire_id) %>% 
                          dplyr::group_split() %>% 
                          purrr::map(runINext) %>% 
                          dplyr::bind_rows()
    rarefaction_tables <- rarefaction_tables %>% 
                          dplyr::mutate(method = recode(method, observed = "Interpolated", 
                                        interpolated = "Interpolated", extrapolated = "Extrapolated")) 
    rarefaction_curves <- ggplot2::ggplot(rarefaction_tables, aes(x=m, y=qD, fill=repertoire_id)) + 
                          ggplot2::geom_line(aes(linetype=method, color=repertoire_id), size=1.5) +  
                          ggplot2::geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL), alpha=0.5) +
                          ggplot2::scale_linetype_manual(values=c("dashed", "solid"),  
                                                         labels=c("Extrapolated", "Intrapolated")) + 
                          ggplot2::theme_classic() +
                          ggplot2::xlab("Total number of sequences") +
                          ggplot2::ylab("TCR diversity") + 
                          ggplot2::labs(fill = "Sample", color="Sample", linetype="Method") 
    return(rarefaction_curves)
}