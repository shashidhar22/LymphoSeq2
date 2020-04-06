#' Lorenz curve
#' 
#' Plots a Lorenz curve derived from the frequency of the amino acid sequences.
#' 
#' @param samples A character vector of sample names in list.
#' @param study_table A tibble generated using the LymphoSeq function readImmunoSeq 
#' or productiveSeq.  "frequencyCount" is a required column.
#' @return Returns a Lorenz curve.
#' @details The Gini coefficient is an alternative metric used to calculate 
#' repertoire diversity and is derived from the Lorenz curve.  The Lorenz curve 
#' is drawn such that x-axis represents the cumulative percentage of unique 
#' sequences and the y-axis represents the cumulative percentage of reads.  A 
#' line passing through the origin with a slope of 1 reflects equal frequencies 
#' of all sequences.  The Gini coefficient is the ratio of the area between the 
#' line of equality and the observed Lorenz curve over the total area under the 
#' line of equality.
#' 
#' The plot is made using the package ggplot2 and can be reformatted
#' using ggplot2 functions.  See examples below.
#' @seealso An excellent resource for examples on how to reformat a ggplot can 
#' be found in the R Graphics Cookbook online (\url{http://www.cookbook-r.com/Graphs/}).
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' study_table <- readImmunoSeq(path = study_table)
#' 
#' samples <- study_table %>% pull(sample) %>% unique()
#'
#' lorenzCurve(samples = samples), study_table = study_table)
#' 
#' productive_aa <- productiveSeq(study_table = study_table, aggregate = "aminoAcid")
#' 
#' samples <- productive_aa %>% pull(sample) %>% unique()
#' 
#' lorenzCurve(samples = samples, study_table = productive_aa)
#'
#' # Change the legend labels, line colors, and add a title
#' samples <- c("TRB_Unsorted_0", "TRB_Unsorted_32", 
#'    "TRB_Unsorted_83", "TRB_Unsorted_949", "TRB_Unsorted_1320")
#' 
#' lorenz_curve <- lorenzCurve(samples = samples, study_table = productive_aa)
#' 
#' labels <- c("Day 0", "Day 32", "Day 83", "Day 949", "Day 1320")
#' 
#' colors <- c("navyblue", "red", "darkgreen", "orange", "purple")
#' 
#' lorenz_curve + ggplot2::scale_color_manual(name = "Samples", breaks = samples, 
#'    labels = labels, values = colors) + ggplot2::ggtitle("Figure Title")
#' @export
#' @import ggplot2
#' @import tidyverse
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ineq Lc
lorenzCurve <- function(samples, study_table) {
    lorenz <- study_table %>% 
              group_by(sample) %>% 
              group_split() %>%
              map(getLorenz) %>%
              row_bind()
    getPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
    plot <- ggplot2::ggplot(lorenz, aes_string(x = "p", y = "L", color = "sample")) + 
        geom_line(size = 1) + 
        theme_minimal() + 
        scale_color_manual(values = getPalette(length(samples) + 1)) + 
        scale_y_continuous(expand = c(0, 0)) + 
        scale_x_continuous(expand = c(0,0)) + 
        geom_abline(intercept = 0, slope = 1, color = "grey", linetype = 2) + 
        coord_fixed() + 
        labs(x = "Cumulative percentage of unique sequences", 
             y = "Cumulative percentage of reads", color = "")
    return(plot)
}
#' Calculate Lorenz curve
#' 
#' Calculate a Lorenz curve derived from the frequency of the amino acid sequences.
#' 
#' @param sample_table A tibble for a single sample generated using the LymphoSeq 
#' function readImmunoSeq or productiveSeq.  "frequencyCount" is a required column.
#' @return Returns a Lorenz curve tibble.
getLorenz <- function(sampl_table) {
    sample <- sample_table$sample[1]
    lc <- ineq::Lc(sample_table$frequencyCount)
    lctbl <- tibble(L = lc$L, p = lc$p) %>% add_column(sample = sample)
    return(lctbl)
}