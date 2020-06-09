#' Align mutliple sequences
#' 
#' Perform multiple sequence alignment using one of three methods and output results to the 
#' console or as a pdf file.  One may perform the alignment of all amino acid or junction
#' sequences in a single repertoire_id.  Alternatively, one may search for a given sequence 
#' within a list of samples using an edit distance threshold.
#' 
#' @param study_table A tibble consisting of antigen receptor sequences 
#' imported by the LymphoSeq function readImmunoSeq.
#' @param repertoire_id A character vector indicating the name of the repertoire_id in the productive
#' sequence list. 
#' @param sequence A character vector of one ore more amino acid or junction 
#' CDR3 sequences to search.
#' @param editDistance An integer giving the minimum edit distance that the 
#' sequence must be less than or equal to.  See details below.
#' @param type A character vector indicating whether "junction_aa" or "junction" sequences
#' should be aligned.  If "junction_aa" is specified, then run productiveSeqs first.
#' @param method A character vector indicating the multiple sequence alignment method to 
#' be used.  Refer to the Bioconductor msa package for more details.  Options incude 
#' "ClustalW", "ClustalOmega", and "Muscle".
#' @param output A character vector indicating where the multiple sequence alignemnt should be
#' printed.  Options include "console" or "pdf".  If "pdf" is selected, the file is saved to
#' the working directory.  For "pdf" to work, Texshade must be installed.  Refer to the 
#' Bioconductor package msa installation instructions for more details.
#' @details Edit distance is a way of quantifying how dissimilar two sequences 
#' are to one another by counting the minimum number of operations required to 
#' transform one sequence into the other.  For example, an edit distance of 0 
#' means the sequences are identical and an edit distance of 1 indicates that 
#' the sequences different by a single amino acid or junction.
#' 
#' @return Performs a multiple sequence alignemnt and prints to the console or saves a pdf to 
#' the working directory.
#' @seealso If having trouble saving pdf files, refer to Biconductor package msa for
#' installation instructions
#' \url{http://bioconductor.org/packages/release/bioc/vignettes/msa/inst/doc/msa.pdf}
#' @examples
#' file.path <- system.file("extdata", "IGH_sequencing", package = "LymphoSeq")
#' 
#' file.list <- readImmunoSeq(path = file.path)
#' 
#' productive.nt <- productiveSeq(file.list = file.list, aggregate = "junction")
#' 
#' alignSeq(list = productive.nt, repertoire_id = "IGH_MVQ92552A_BL", type = "junction", 
#'          method = "ClustalW", output = "console")
#' @export
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings AAStringSet
#' @import msa
#' @import tidyverse
alignSeq = function(study_table, sample_list = NULL, sequence_list = NULL, editDistance = 15, output = "console", type = "junction", method = "ClustalOmega") {
    if(!is.null(sequence_list) & is.null(sample_list)){
        search_table <- searchSeq(study_table = study_table, sequence = sequence_list, editDistance = editDistance, type = type, match = "partial")
        if(is.null(search_table)){
            stop("There are no sequences to be aligned", call. = FALSE)
        }
    }
    if(!is.null(sequence_list) & !is.null(sample_list)){
        study_table %>% study_table %>% filter(repertoire_id %in% sample_list)
        search_table <- searchSeq(list = list[repertoire_id], sequence = sequence_list, editDistance = editDistance, type = type, match = "partial")
        if(is.null(search_table)){
            stop("There are no sequences to be aligned", call. = FALSE)
        }
    }
    if(is.null(sequence_list)){
        search_table <- study_table %>% filter(repertoire_id %in% sample_list)
    }
    if(nrow(search_table) > 150){
        if(is.null(sequence_list)){
            search_table <- search_table %>% top_n(150) 
        } else {
            search_table <- search_table %>% group_by(searchSequence) %>% top_n(150) %>% ungroup()
        }
        message("Only 150 sequences sampled equally from each search group will be selected")
    }
    if(type == "junction"){
        search_table <- search_table %>% filter(nchar(junction) > 15)
        string_list <- Biostrings::DNAStringSet(search_table$junction)
    }
    if(type == "junction_aa"){
        search_table <- search_table %>% filter(nchar(junction_aa) > 3)
        string_list <- Biostrings::AAStringSet(search_table$junction_aa)
    }
    if(nrow(search_table) < 3){
        stop("There are less than 3 sequences to be aligned", call. = FALSE)
    }
    if(!is.null(sequence_list)){
        names(string_list) <- paste(search_table$repertoire_id)
    } else {
        names(string_list) <- paste(search_table$v_family, search_table$d_family, 
                                    search_table$j_family, search_table$duplicate_count, sep="_")
    }
    names(string_list) <- gsub(names(string_list), pattern = "IGH|IGL|IGK|TCRB|TCRA", replacement = "")
    names(string_list) <- gsub(names(string_list), pattern = "NA|unresolved", replacement = "UNR")
    alignment <- msa::msa(string_list, method = method)
    if(output == "console"){
        print(alignment, show = "complete")
    }
    if(output == "pdf"){
        file.name <- "Sequence_aligment.pdf"
        msa::msaPrettyPrint(alignment, output = "pdf", file = file.name, 
                            showNumbering = "right", showNames = "left", showConsensus = "top",  
                            shadingMode = "similar", shadingColors = "blues",
                            paperWidth = ncol(alignment)*0.1 + 5, paperHeight = nrow(alignment)*0.1 + 5,
                            showLogo = "bottom", showLogoScale = "right", logoColors = "rasmol",
                            psFonts = TRUE, showLegend = TRUE, askForOverwrite = FALSE,
                            furtherCode=c("\\defconsensus{.}{lower}{upper}","\\showruler{1}{top}"))
       message(paste(file.name, "saved to", getwd()))
    }
}
