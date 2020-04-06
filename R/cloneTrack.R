#' Clone tracking plot
#' 
#' Creates line plot tracking amino acid frequencies across multiple samples
#' 
#' @param sequence_matrix A sequence matrix generated from the LymphoSeq 
#' function seqMatrix.
#' @param map An optional character vector of one or more sample names contained 
#' in the productive_aa tibble.  If the same sequence appears in multiple mapped 
#' samples, then it will be assigned the label of the first listed sample only.
#' @param productive.aa A list of data frames of productive amino acid sequences 
#' containing the samples to be mapped.  This parameter is only required if 
#' sequence mapping is performed.
#' @param label An optional character vector of one or more labels used to 
#' annotate the mapped sequences.  The order of the labels must match the order 
#' of the samples listed in map.
#' @param track An optional character vector of one or more amino acid sequences 
#' to track.
#' @param unassigned A Boolean value indicating whether or not to draw the lines 
#' of sequences not being mapped or tracked.  If TRUE then the unassigned 
#' sequences are drawn.  If FALSE, then the unassigned sequences are not drawn.
#' @return Returns a line plot showing the amino acid frequencies across 
#' multiple samples in the sequence matrix where each line represents one 
#' unique sequence.
#' @details The plot is made using the package ggplot2 and can be reformatted
#' using ggplot2 functions.  See examples below.
#' @seealso An excellent resource for examples on how to reformat a ggplot can 
#' be found in the R Graphics Cookbook online (\url{http://www.cookbook-r.com/Graphs/}).
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' study_table <- readImmunoSeq(path = file.path)
#' 
#' productive_aa <- productiveSeq(study_table = study_table, aggregate = "aminoAcid")
#' 
#' top_freq <- topFreq(productive_table = productive_aa, percent = 0.1)
#' 
#' sequence_matrix <- seqMatrix(productive_table = productive_aa, sequences = top.freq$aminoAcid)
#' 
#' # Track clones without mapping or tracking specific sequences
#' cloneTrack(sequence_matrix = sequence_matrix)
#' 
#' # Track top 20 clones mapping to the CD4 and CD8 samples
#' cloneTrack(sequence_matrix = sequence_matrix, productive_table =  productive_aa, 
#'    map = c("TRB_CD4_949", "TRB_CD8_949"), label = c("CD4", "CD8"), 
#'    track = top_freq$aminoAcid[1:20], unassigned = TRUE) 
#' 
#' # Track the top 10 clones from top.freq
#' cloneTrack(sequence_matrix = sequence_matrix, productive_table = productive_aa, 
#'    track = top_freq$aminoAcid[1:10], unassigned = FALSE) 
#' 
#' # Track clones mapping to the CD4 and CD8 samples while ignoring all others
#' cloneTrack(sequence_matrix = sequence_matrix, productive_table = productive_aa, 
#'    map = c("TRB_CD4_949", "TRB_CD8_949"), label = c("CD4", "CD8"), 
#'    unassigned = FALSE) 
#' 
#' # Track clones mapping to the CD4 and CD8 samples and track 2 specific sequences
#' cloneTrack(sequence_matrix = sequence_matrix, productive_table = productive_aa, 
#'    map = c("TRB_CD4_949", "TRB_CD8_949"), label = c("CD4", "CD8"), 
#'    track = c("CASSPPTGERDTQYF", "CASSQDRTGQYGYTF"), unassigned = FALSE)
#' 
#' # Reorder the x axis, change the axis labels, convert to log scale, and add title
#' x.limits <- c("TRB_Unsorted_0", "TRB_Unsorted_32", 
#'    "TRB_Unsorted_83", "TRB_Unsorted_949", "TRB_Unsorted_1320")
#' 
#' sequence.matrix <- sequence.matrix[ ,c("aminoAcid", x.limits)]
#'    
#' clone.track <- cloneTrack(sequence.matrix = sequence.matrix, 
#'    productive.aa = productive.aa, track = top.freq$aminoAcid[1:10], unassigned = FALSE) 
#' 
#' x.labels <- c("Day 0", "Day 32", "Day 83", "Day 949", "Day 1320")
#' 
#' clone.track + 
#'    ggplot2::scale_x_discrete(expand = c(0,0), labels = x.labels) + 
#'    ggplot2::scale_y_log10() + ggplot2::annotation_logticks(sides = "l") + 
#'    ggplot2::ggtitle("Figure Title")
#' @export
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom reshape melt.data.frame
#' @import tidyverse
cloneTrack <- function(sequence_matrix, sample_map = "none", productive_aa, 
                       label = "none", sequence_track = "none", unassigned = TRUE) {
    amino_list <- rownames(sequence_matrix)
    if (label == "none" & sample_map != "none"){
        label <- sample_map
        names(label) <- sample_map
    } else {
        names(label) <- sample_map
    }
    search_space <- sequence_matrix %>% 
                    as_tibble() %>% 
                    select(-numberSamples) %>%
                    pivot_longer(names_to = "Sample", values_to="Frequency") %>%
                    mutate(Map = if_else(sample %in% sample_map, label[sample], "Unassigned")) %>%
                    mutate(Map = if_else(aminoAcid %in% sequence_track, aminoAcid, Map)) 

    if (unassigned == FALSE) {
        search_space <- search_space %>% filter(sample_guide != "Unassigned")
    }
    getPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
    plot <- ggplot2::ggplot(sequence.melt, 
                            aes_string(x = "Sample", 
                                       y = "Frequency", 
                                       group = "aminoAcid", 
                                       color = "Map")) + 
            geom_line() + 
            theme_minimal() + 
            scale_x_discrete(expand = c(0, 0)) + 
            scale_y_continuous(expand = c(0, 0)) + 
            labs(x = "", y = "Frequency (%)", color = "") + 
            scale_colour_manual(values = c(getPalette(length(label) + 
                                                      length(track) + 1))) + 
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
    return(plot)
}
  