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

rarefactionCurve <- function(study_table) {
    progress <- progress_estimated(length(unique(study_table$sample)))
    rarefaction_tables <- study_table %>% group_by(sample) %>% group_split() %>% map(~runINext(.x, progress)) %>% bind_rows()
    rarefaction_tables <- rarefaction_tables %>% 
                          mutate(method = recode(method, observed = "Interpolated", interpolated = "Interpolated", extrapolated = "Extrapolated")) 
    rarefaction_curves <- ggplot(rarefaction_tables, aes(x=m, y=qD, fill=sample)) + 
                          geom_line(aes(linetype=method, color=sample), size=1.5) +  
                          geom_ribbon(aes(ymin = qD.LCL, ymax = qD.UCL), alpha=0.5) +
                          scale_linetype_manual(values=c("dashed", "solid")) + 
                          theme_classic() +
                          xlab("Total number of sequences") +
                          ylab("TCR diversity") + 
                          labs(fill = "Sample", color="Sample", linetype="Method") 
    return(rarefaction_curves)
}
