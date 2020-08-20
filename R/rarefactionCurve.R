#' Plot rarefication and extrapolation curves for samples
#' 
#' Given a study table, for each sample plot rarefaction curves to estimate 
#' repertoire diversity. The method used to generate the rarefaction curve
#' is derived from Chao et al., (2014) using the iNEXT library
#'
#' @param study_table A tibble consisting antigen receptor sequencing 
#' data imported by the LymphoSeq function readImmunoSeq. "aminoAcid", "count", and 
#' "frequencyCount" are required columns.
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' study_table <- readImmunoSeq(path = file.path)
#' 
#' productive_aa <- productiveSeq(study_table = study_table, 
#'   aggregate = "aminoAcid", prevalence = TRUE)
#'
#' rarefaction_curve <- rarefactionCurve(study_table)
#'
#' @export
#' @import tidyverse
#' @import purrr

plotRarefactionCurve <- function(study_table) {
    progress <- dplyr::progress_estimated(length(unique(study_table$sample)))
    rarefaction_tables <- study_table %>% 
                          dplyr::group_by(repertoire_id) %>% 
                          dplyr::group_split() %>% 
                          dplyr::map(~runINext(.x, progress)) %>% 
                          dplyr::bind_rows()
    rarefaction_tables <- rarefaction_tables %>% 
                          dplyr::mutate(method = recode(method, observed = "Interpolated", 
                                        interpolated = "Interpolated", extrapolated = "Extrapolated")) 
    rarefaction_curves <- ggplot2::ggplot(rarefaction_tables, aes(x=m, y=qD, fill=sample)) + 
                          ggplot2::geom_line(aes(linetype=method, color=sample), size=1.5) +  
                          ggplot2::geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL), alpha=0.5) +
                          ggplot2::scale_linetype_manual(values=c("dashed", "solid"),  
                                                         labels=c("Extrapolated", "Intrapolated")) + 
                          ggplot2::theme_classic() +
                          ggplot2::xlab("Total number of sequences") +
                          ggplot2::ylab("TCR diversity") + 
                          ggplot2::labs(fill = "Sample", color="Sample", linetype="Method") 
    return(rarefaction_curves)
}
