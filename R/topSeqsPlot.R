#' Cumulative frequency bar plot of top sequences
#' 
#' Create a cumulative frequency bar plot of a specified number of top 
#' sequences.
#' 
#' @param study_table A study tibble imported using the LymphoSeq function readImmunoSeq 
#' or productiveSeq.
#' @param top The number of top sequences to be colored in the bar plot.  All 
#' other, less frequent sequences are colored violet.
#' @return Returns a cumulative frequency bar plot of the top sequences.
#' @details The plot is made using the package ggplot2 and can be reformatted
#' using ggplot2 functions.  See examples below.
#' @seealso An excellent resource for examples on how to reformat a ggplot can 
#' be found in the R Graphics Cookbook online (\url{http://www.cookbook-r.com/Graphs/}).
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' study_table <- readImmunoSeq(path = file.path)
#' 
#' productive_aa <- productiveSeq(study_table = study_table, aggregate = "junction_aa")
#' 
#' topSeqsPlot(study_table = productive_aa, top = 10)
#' 
#' # Display the number of sequences at the top of bar plot and add a title
#' n <- as.character(nrow(study_table))
#' 
#' topSeqsPlot(study_table = productive_aa, top = 10) + 
#'    ggplot2::annotate("text", x = 1:length(file.list), y = 105, label = n, color = "black") +
#'    ggplot2::expand_limits(y = c(0, 110)) + ggplot2::ggtitle("Figure Title") + 
#'    ggplot2::scale_x_discrete(limits = names(file.list))
#' @export
#' @import ggplot2
#' @import tidyverse
#' @importFrom RColorBrewer brewer.pal
#' @importFrom plyr llply
topSeqsPlot <- function(study_table, top = 10) {
    dominant <- study_table %>% 
                dplyr::group_by(repertoire_id) %>% 
                dplyr::arrange(desc(duplicate_frequency)) %>% 
                dplyr::top_n(top, wt=duplicate_frequency) %>%
                dplyr::select(repertoire_id, junction_aa, duplicate_frequency) %>%
                dplyr::ungroup()
                
    
    subdominant <- dominant %>% 
                   dplyr::group_by(repertoire_id) %>% 
                   dplyr::summarize(duplicate_frequency = 100 - sum(duplicate_frequency), 
                                    junction_aa = "All other sequences") %>%
                   dplyr::select(repertoire_id, junction_aa, duplicate_frequency)
    topfreq <- dplyr::bind_rows(dominant,subdominant) %>% 
               dplyr::rename(repertoire_id=repertoire_id, 
                             junction_aa= junction_aa, 
                             Frequency=duplicate_frequency) %>%
               dplyr::group_by(repertoire_id) %>% 
               dplyr::mutate(Sequence = 1:n()) %>%
               dplyr::mutate(Sequence = as.factor(Sequence)) %>%
               dplyr::ungroup() %>%
               dplyr::arrange(repertoire_id, Sequence, desc(Frequency))
    getPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))
    sample_order <- subdominant %>%
                    dplyr::arrange(duplicate_frequency) %>% 
                    dplyr::select(repertoire_id) %>% 
                    dplyr::pull()
    ggplot2::ggplot(topfreq, aes_string(x = "repertoire_id", 
                               y = "Frequency", 
                               fill = "Sequence", 
                               label = "Frequency", 
                               text = "junction_aa")) + 
    ggplot2::geom_bar(stat = "identity") + 
    ggplot2::scale_x_discrete(limits = sample_order) + 
    ggplot2::scale_fill_manual(values = getPalette(top + 1)) + 
    ggplot2::theme_classic() + 
    ggplot2::scale_y_continuous(expand = c(0, 0)) + 
    ggplot2::theme(legend.position = "none") + 
    ggplot2::theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), 
                   axis.text.y = element_text(size = 10)) + labs(x = "", y = "Frequency (%)")
} 
