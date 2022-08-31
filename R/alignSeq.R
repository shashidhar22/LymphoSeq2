#' Align multiple sequences
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
#' @param sequence_list A character vector of one ore more amino acid or junction 
#' CDR3 sequences to search.
#' @param edit_distance An integer giving the minimum edit distance that the 
#' sequence must be less than or equal to.  See details below.
#' @param type A character vector indicating whether "junction_aa" or "junction" sequences
#' should be aligned.  If "junction_aa" is specified, then run productiveSeqs first.
#' @param method A character vector indicating the multiple sequence alignment method to 
#' be used.  Refer to the Bioconductor msa package for more details.  Options include 
#' "ClustalW", "ClustalOmega", and "Muscle".
#' @param top The number of top sequences to perform alignment.
#' @details Edit distance is a way of quantifying how dissimilar two sequences 
#' are to one another by counting the minimum number of operations required to 
#' transform one sequence into the other.  For example, an edit distance of 0 
#' means the sequences are identical and an edit distance of 1 indicates that 
#' the sequences different by a single amino acid or junction.
#' @return Performs a multiple sequence alignment object.
#' @seealso If having trouble saving pdf files, refer to Bioconductor package msa for
#' installation instructions
#' \url{http://bioconductor.org/packages/release/bioc/vignettes/msa/inst/doc/msa.pdf}
#' @examples
#' file_path <- system.file("extdata", "IGH_sequencing", package = "LymphoSeq2")
#' 
#' stable <- readImmunoSeq(path = file_path)
#' 
#' ntable <- productiveSeq(stable, aggregate = "junction")
#' 
#' alignSeq(ntable, repertoire_id = "IGH_MVQ92552A_BL", type = "junction", 
#'          method = "ClustalW")
#' @export
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings AAStringSet
#' @import msa
#' @import magrittr
alignSeq <- function(study_table, repertoire_ids = NULL, 
                    sequence_list = NULL, edit_distance = 15, 
                    type = "junction", method = "ClustalOmega", 
                    top = 150) {
    #If the sequence list is not null and repertoire ids are NULL,
    #set repertoire_ids to all samples and align all sequences against 
    #sequence list
    if(!is.null(sequence_list) & is.null(repertoire_ids)){
        search_table <- LymphoSeq2::searchSeq(study_table = study_table, 
            sequence = sequence_list, 
            edit_distance = edit_distance, 
            seq_type = type, 
            match = "partial")
        repertoire_ids <- study_table %>%
            dplyr::pull(repertoire_id) %>%
            base::unique()
        if(is.null(search_table)){
            stop("There are no sequences to be aligned", call. = FALSE)
        }
    } else if(!is.null(sequence_list) & !is.null(repertoire_ids)){
        study_table <- study_table %>% 
            dplyr::filter(repertoire_id %in% repertoire_ids)
        search_table <- LymphoSeq2::searchSeq(study_table = study_table,
            sequence = sequence_list, 
            edit_distance =  edit_distance, 
            seq_type = type, 
            match = "partial")
        if(is.null(search_table)){
            stop("There are no sequences to be aligned", call. = FALSE)
        }
    } else if(is.null(sequence_list) & !is.null(repertoire_ids)){
        search_table <- study_table %>% 
                        dplyr::filter(repertoire_id %in% repertoire_ids)
    } else if(is.null(sequence_list) & is.null(repertoire_ids)) {
        search_table <- study_table
        repertoire_ids <- study_table %>% 
            dplyr::pull(repertoire_id) %>%
            base::unique()
    }
    #Select only the top sequences to perform alignment
    if(nrow(search_table) > top){
        if(is.null(sequence_list)){
            search_table <- search_table %>% 
                            LymphoSeq2::topSeqs(top = top)
        } else {
            search_table <- search_table %>% 
                            dplyr::filter(!!base::as.symbol(type) %in% searchSequence) %>% 
                            LymphoSeq2::topSeqs(top = top)
        }
        message("Only 150 sequences sampled equally from each search group will be selected")
    }
    #Select the string set according to the type provided by the user
    #and create a DNA/AAStringSet with this
    if(type == "junction"){
        search_table <- search_table %>% 
                        dplyr::filter(base::nchar(junction) > 15)
        string_list <- search_table %>% 
            dplyr::pull(junction) %>%
            Biostrings::DNAStringSet()
    }
    if(type == "junction_aa"){
        search_table <- search_table %>% 
                        dplyr::filter(base::nchar(junction_aa) > 3)
        string_list <- search_table %>% 
            dplyr::pull(junction_aa) %>%
            Biostrings::AAStringSet()
    }
    #Stop if the number of sequences are less than three
    if(nrow(search_table) < 3){
        stop("There are less than 3 sequences to be aligned", call. = FALSE)
    }
    #Prepare a string list to name annotate the sequences with repertoire_ids
    #if provided or use gene family names with sequence counts 
    if(!is.null(repertoire_ids)){
        names(string_list) <- search_table %>%
            dplyr::pull(repertoire_id) 
    } else {
        names(string_list) <- search_table %>% 
            dplyr::mutate(seq_name = base::paste(v_family, d_family, j_family, 
                    duplicate_count, sep="_"),
                seq_name = stringr::str_replace(seq_name, 
                    "IGH|IGL|IGK|TCRB|TCRA", ""),
                seq_name = stringr::str_replace(seq_name, "unresolved", "UNR"),
                seq_name = stringr::str_replace_na(seq_name, 
                    replacement = "UNR")) %>%
            dplyr::pull(seq_name)
    }
    #Perform multiple sequence alignment using the method described by the user
    base::set.seed(12357)
    alignment <- msa::msa(string_list, method = method)
    return(alignment)
}


