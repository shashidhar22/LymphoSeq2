#' Lorenz curve
#' 
#' Plots a Lorenz curve derived from the frequency of the amino acid sequences.
#' 
#' @param repertoire_ids A character vector of repertoire_id names in list.
#' @param study_table A tibble generated using the LymphoSeq2 function 
#' [readImmunoSeq()] or [productiveSeq()]. "duplicate_frequency" is a required 
#' column.
#' @return Returns a Lorenz curve.
#' @details The Gini coefficient is an alternative metric used to calculate 
#' repertoire diversity and is derived from the Lorenz curve.  The Lorenz curve 
#' is drawn such that x-axis represents the cumulative percentage of unique 
#' sequences and the y-axis represents the cumulative percentage of reads.  A 
#' line passing through the origin with a slope of 1 reflects equal frequencies 
#' of all sequences.  The Gini coefficient is the ratio of the area between the 
#' line of equality and the observed Lorenz curve over the total area under the 
#' line of equality.
#' The plot is made using the package ggplot2 and can be reformatted
#' using ggplot2 functions.  See examples below.
#' @seealso An excellent resource for examples on how to reformat a ggplot can 
#' be found in the R Graphics Cookbook online
#' (\url{http://www.cookbook-r.com/Graphs/}).
#' @examples
#' library(magrittr)
#' file_path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq2")
#' # Plot Lorenz curve with raw data
#' study_table <- LymphoSeq2::readImmunoSeq(path = file_path, threads = 1) %>% 
#' LymphoSeq2::topSeqs(top = 100)
#' repertoire_ids <- study_table %>% 
#'   dplyr::pull(repertoire_id) %>% 
#'   unique()
#' LymphoSeq2::lorenzCurve(repertoire_ids = repertoire_ids, study_table = study_table)
#' # Plot Lorenz curve with productive amino acid sequences
#' amino_table <- LymphoSeq2::productiveSeq(study_table = study_table, 
#'   aggregate = "junction_aa")
#' repertoire_ids <- amino_table %>%
#'   dplyr::pull(repertoire_id) %>%
#'   unique()
#' LymphoSeq2::lorenzCurve(repertoire_ids = repertoire_ids, study_table = amino_table)
#' # Change the legend labels, line colors, and add a title
#' repertoire_ids <- c("TRB_Unsorted_0", "TRB_Unsorted_32", 
#'    "TRB_Unsorted_83", "TRB_Unsorted_949", "TRB_Unsorted_1320")
#' lorenz_curve <- LymphoSeq2::lorenzCurve(repertoire_ids = repertoire_ids, 
#'   study_table = amino_table)
#' labels <- c("Day 0", "Day 32", "Day 83", "Day 949", "Day 1320")
#' colors <- c("navyblue", "red", "darkgreen", "orange", "purple",
#'              "yellow", "pink", "lightgreen", "cyan", "maroon")
#' lorenz_curve + 
#' ggplot2::scale_color_manual(name = "repertoire_ids", 
#'                             breaks = repertoire_ids, 
#'                             labels = labels, values = colors) + 
#' ggplot2::ggtitle("Lorenz curve")
#' @export
#' @import magrittr
lorenzCurve <- function(repertoire_ids, study_table) {
    lorenz <- study_table %>% 
              dplyr::group_by(repertoire_id) %>% 
              dplyr::group_split() %>%
              purrr::map(LymphoSeq2::getLorenz) %>%
              dplyr::bind_rows()
    getPalette <- grDevices::colorRampPalette(
      RColorBrewer::brewer.pal(9, "Set1"))
    plot <- ggplot2::ggplot(lorenz, ggplot2::aes_string(x = "p", y = "L", 
                                               color = "repertoire_id")) + 
            ggplot2::geom_line(size = 1) + 
            ggplot2::theme_minimal() + 
            ggplot2::scale_color_manual(values = getPalette(
              length(repertoire_ids) + 1)) + 
            ggplot2::scale_y_continuous(expand = c(0, 0)) + 
            ggplot2::scale_x_continuous(expand = c(0,0)) + 
            ggplot2::geom_abline(intercept = 0, slope = 1, color = "grey", 
                                 linetype = 2) + 
            ggplot2::coord_fixed() + 
            ggplot2::labs(x = "Cumulative percentage of unique sequences", 
                          y = "Cumulative percentage of reads", color = "")
    return(plot)
}
