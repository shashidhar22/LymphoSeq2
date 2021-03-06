#' Read ImmunoSeq files
#'
#' \code{readImmunoSeq} Imports tab-separated value (.tsv) files exported by the Adaptive 
#' Biotechnologies ImmunoSEQ analyzer, BGI IR-SEQ, MiXCR and stores 
#' them as MiAIRR compliant tibble. 
#' 
#' @details May import tab-delimited files containing antigen receptor 
#' sequencing from with the following header set.  
#' 
#' | MiAIRR field ID | Type1           | Type2                       | Type3               | Type4                          | Type5                               | Type6             |
#' | --------------- | --------------- | --------------------------- | ------------------- | ------------------------------ | ----------------------------------- | ----------------- |
#' | junction        | nucleotide      | nucleotide                  | nucleotide          | nucleotide.CDR3.in.lowercase . | nucleotide( CDR3   in   lowercase ) | nSeqCDR3          |
#' | junction_aa     | aminoAcid       | aminoAcid                   | aminoAcid           | aminoAcid.CDR3.in.lowercase .  | aminoAcid( CDR3   in   lowercase )  | aaSeqCDR3         |
#' | duplicate_count | count ( reads ) | count ( templates / reads ) | count ( templates ) | cloneCount                     | cloneCount                          | cloneCount        |
#' | v_call          | vGeneName       | vGeneName                   | vGeneName           | vGene                          | vGene                               | allVHitsWithScore |
#' | d_call          | dGeneName       | dGeneName                   | dGeneName           | dGene                          | dGene                               | allDHitsWithScore |
#' | j_call          | jGeneName       | jGeneName                   | jGeneName           | jGene                          | jGene                               | allJHitsWithScore |
#' 
#' 
#' @param path Path to AIRR-seq files. Can be a file path, directory path or a vector of file paths
#' @param recursive Recursively search directory for AIRR-seq files. Boolean value
#' @return Returns a tibble with MiAIRR headers and repertoire_id
#'
#' @examples
#' file_path <- base::system.file("extdata", "TCRB_study", package = "LymphoSeq2")
#' stable <- LymphoSeq2::readImmunoSeq(path = file_path, 
#'                                          recursive = FALSE)
#' file_list <- base::list.files(file_path, 
#'                               full.name = TRUE,
#'                               all.files = FALSE,
#'                               pattern = ".tsv",
#'                               include.dirs = FALSE)
#' stable <- LymphoSeq2::readImmunoSeq(path = file_list[1],
#'                                          recursive = FALSE)
#' stable <- LymphoSeq2::readImmunoSeq(path = file_list[1:2],
#'                                          recursive = FALSE)
#' @export
#' @rdname readImmunoSeq
readImmunoSeq <- function(path, recursive = FALSE) {
    if (length(path) > 1) {
        file_paths <- path
    } else if (file_test("-d", path)) {
        file_paths <- list.files(path, 
                                 full.names = TRUE, 
                                 all.files = FALSE, 
                                 recursive = recursive, 
                                 pattern = ".tsv|.txt|.tsv.gz", 
                                 include.dirs = FALSE)
    } else {
        file_paths <- c(path)
    }
    file_num <- length(file_paths)
    file_info <- file.info(file_paths)
    file_paths <- rownames(file_info)[file_info$size > 0]
    if(file_num != length(file_paths)){
        warning("One or more of the files you are trying to import has no sequences and will be ignored.", 
                call. = FALSE)
    }
    progress_bar <- progress::progress_bar$new(format = "Reading AIRR-Seq files [:bar] :percent eta: :eta",
                                           total = length(file_paths), clear = FALSE, width = 60)
    progress_bar$tick(0)
    file_list <- file_paths %>% 
                 purrr::map(~ readFiles(.x, progress_bar)) %>% 
                 dplyr::bind_rows()
    progress_bar$terminate()
    return(file_list)
}

#' Identify file type:
#'
#' AIRRSeq results have various formats, this functions identifies six formats 
#' 
#' @param clone_file .tsv file containing results from AIRRSeq pipeline
#'
#' @return List containing file type and standardized header names 
#'
getFileType <- function(clone_file) {
    type1 <- c("nucleotide", "aminoAcid", "count (reads)", "frequencyCount (%)",
               "vGeneName", "dGeneName", "jGeneName", "vFamilyName", 
               "dFamilyName", "jFamilyName", "sequenceStatus", 
               "estimatedNumberGenomes")
        
    type2 <- c("nucleotide", "aminoAcid", "count (templates/reads)", 
               "frequencyCount (%)", "vGeneName", "dGeneName", "jGeneName", 
               "vFamilyName", "dFamilyName", "jFamilyName", "sequenceStatus", 
               "estimatedNumberGenomes")
    
    type3 <- c("nucleotide", "aminoAcid", "count (templates)", 
               "frequencyCount (%)", "vGeneName", "dGeneName", "jGeneName", 
               "vMaxResolved", "dMaxResolved", "jMaxResolved", "sequenceStatus", 
               "estimatedNumberGenomes")
    
    type4 <- c("nucleotide", "aminoAcid", "count (templates)", 
               "frequencyCount (%)", "vGeneName", "dGeneName", "jGeneName", 
               "vFamilyName", "dFamilyName", "jFamilyName", "sequenceStatus", 
                "estimatedNumberGenomes")

    type5 <- c("nucleotide", "aminoAcid", "count", "frequencyCount", 
               "vGeneName", "dGeneName", "jGeneName", "vFamilyName", 
               "dFamilyName", "jFamilyName", "sequenceStatus", 
               "estimatedNumberGenomes")
    
    type6 <- c("nucleotide.CDR3.in.lowercase.", "aminoAcid.CDR3.in.lowercase.", 
               "cloneCount", "clonefrequency....", "vGene", "dGene", "jGene", 
               "fuction")

    type7 <- c("nucleotide(CDR3 in lowercase)", "aminoAcid(CDR3 in lowercase)", 
               "cloneCount", "clonefrequency (%)", "vGene", "dGene", "jGene", 
               "fuction")
    
    type8 <- c("cloneCount", "cloneFraction", "allVHitsWithScore", 
               "allDHitsWithScore", "allJHitsWithScore", "nSeqCDR3", 
               "aaSeqCDR3")

    type9 <- c("duplicate_count", "v_call", "d_call", "j_call", "junction", 
               "junction_aa")

    type10 <- c("nucleotide", "aminoAcid", "count", 
               "frequencyCount (%)", "vGeneName", "dGeneName", "jGeneName", 
               "vFamilyName", "dFamilyName", "jFamilyName", "sequenceStatus", 
                "estimatedNumberGenomes")

    type11 <- c("nucleotide", "aminoAcid", "count", 
               "frequencyCount", "vGeneName", "dGeneName", "jGeneName", 
               "vFamilyName", "dFamilyName", "jFamilyName", "sequenceStatus", 
                "estimatedNumberGenomes")


    columns <- invisible(colnames(readr::read_tsv(clone_file, 
                                                  n_max = 1, 
                                                  col_types = readr::cols())))
    if (all(type1 %in% columns)) {
        file_type <- "type1"
        header_list <- readr::cols_only(nucleotide = "c", 
                                        aminoAcid = "c", 
                                        `count (reads)` = "i", 
                                        vGeneName = "c", 
                                        dGeneName = "c", 
                                        jGeneName = "c")
    } else if (all(type2 %in% columns)) {
        file_type <- "type2"
        header_list <- readr::cols_only(nucleotide = "c", 
                                        aminoAcid = "c", 
                                        `count (templates/reads)` = "i", 
                                        vGeneName = "c", 
                                        dGeneName = "c", 
                                        jGeneName = "c")
    } else if (all(type3 %in% columns)) {
        file_type <- "type3"
        header_list <- readr::cols_only(nucleotide = "c", 
                                       aminoAcid = "c", 
                                       `count (templates)` = "i", 
                                       vMaxResolved = "c", 
                                       dMaxResolved = "c", 
                                       jMaxResolved = "c")
    } else if (all(type4 %in% columns)) {
        file_type <- "type4"
        header_list <- readr::cols_only(nucleotide = "c", 
                                        aminoAcid = "c", 
                                        `count (templates)` = "i", 
                                        vFamilyName = "c", 
                                        dGeneName = "c", 
                                        jGeneAllele = "c")
    } else if (all(type5 %in% columns)) {
        file_type <- "type5"
        header_list <- readr::cols_only(nucleotide = "c", 
                                        aminoAcid = "c", 
                                        `count` = "i", 
                                        vGeneName = "c", 
                                        dGeneName = "c", 
                                        jGeneName = "c")
    } else if (all(type6 %in% columns)) {
        file_type <- "type6"
        header_list <- readr::cols_only(`nucleotide.CDR3.in.lowercase.` = "c", 
                                        `aminoAcid.CDR3.in.lowercase.` = "c", 
                                        cloneCount = "i", 
                                        `clonefrequency....` = "d",
                                        vGene = "c", 
                                        dGene = "c", 
                                        jGene = "c")
    } else if (all(type7 %in% columns)) {
        file_type <- "type7"
        header_list <- readr::cols_only(`nucleotide(CDR3 in lowercase)` = "c", 
                                        `aminoAcid(CDR3 in lowercase)` = "c", 
                                        cloneCount = "i", 
                                        `clonefrequency (%)` = "d", 
                                        vGene = "c", 
                                        dGene = "c", 
                                        jGene = "c")
    } else if (all(type8 %in% columns)) {
        file_type <- "type8"
        header_list <- readr::cols_only(cloneCount = 'd', 
                                        cloneFraction = 'd', 
                                        allVHitsWithScore = 'c', 
                                        allDHitsWithScore = 'c', 
                                        allJHitsWithScore = 'c', 
                                        nSeqCDR3 = 'c', 
                                        aaSeqCDR3 = 'c')
    } else if (all(type9 %in% columns)) {
        file_type <- "type9"
        header_list <- readr::cols_only(duplicate_count = 'd', 
                                        v_call = 'c', 
                                        d_call = 'c', 
                                        j_call = 'c', 
                                        junction = 'c', 
                                        junction_aa = 'c')
    } else if (all(type10 %in% columns)) {
        file_type <- "type4"
        header_list <- readr::cols_only(nucleotide = "c", 
                                        aminoAcid = "c", 
                                        count = "i", 
                                        vFamilyName = "c", 
                                        dGeneName = "c", 
                                        jGeneAllele = "c")
    } else if (all(type11 %in% columns)) {
        file_type <- "type4"
        header_list <- readr::cols_only(nucleotide = "c", 
                                        aminoAcid = "c", 
                                        count = "i", 
                                        vFamilyName = "c", 
                                        dGeneName = "c", 
                                        jGeneAllele = "c")
    }
    ret_val <- list(file_type, header_list)
    return(ret_val)
}

#' @section Get standard headers:
#'
#' Retrives MiAIRR standard compliant fields from the clone files
#'
#' @param file_type file type from getFileType method
#'
#' @return List of standard header names 
#' 
#' @rdname readImmunoSeq
getStandard <- function(file_type) {
    type1 <- c(junction = "nucleotide", junction_aa = "aminoAcid", 
               duplicate_count = "count (reads)", v_call = "vGeneName", 
               d_call = "dGeneName", j_call = "jGeneName")

    type2 <- c(junction = "nucleotide", junction_aa = "aminoAcid", 
               duplicate_count = "count (templates/reads)",
               v_call = "vGeneName", d_call = "dGeneName", j_call = "jGeneName")
    
    type3 <- c(junction = "nucleotide", junction_aa = "aminoAcid", 
               duplicate_count = "count (templates)", 
               v_call = "vMaxResolved", d_call = "dMaxResolved", j_call = "jMaxResolved")
    

    type4 <- c(junction = "nucleotide", junction_aa = "aminoAcid", 
               duplicate_count = "count (templates)", 
               v_call = "vFamilyName", d_call = "dGeneName", j_call = "jGeneAllele")

    type5 <- c(junction = "nucleotide", junction_aa = "aminoAcid", 
               duplicate_count = "count", v_call = "vGeneName", d_call = "dGeneName", 
               j_call = "jGeneName")

    type6 <- c(junction = "nucleotide.CDR3.in.lowercase.", junction_aa = "aminoAcid.CDR3.in.lowercase.", 
               duplicate_count = "cloneCount", v_call = "vGene", d_call = "dGene", j_call = "jGene")   
    
    type7 <- c(junction = "nucleotide(CDR3 in lowercase)", junction_aa = "aminoAcid(CDR3 in lowercase)", 
               duplicate_count = "cloneCount", v_call = "vGene", d_call = "dGene", j_call = "jGene")   
    
    type8 <- c(junction = "nSeqCDR3", junction_aa = "aaSeqCDR3", duplicate_count = "cloneCount",
               v_call = "allVHitsWithScore", j_call = "allJHitsWithScore", d_call = "allDHitsWithScore")

    type9 <- c(junction = "junction", junction_aa = "junction_aa", duplicate_count = "duplicate_count",
               v_call = "v_call", j_call = "j_call", d_call = "d_call")

    type10 <- c(junction = "nucleotide", junction_aa = "aminoAcid", 
               duplicate_count = "count", 
               v_call = "vFamilyName", d_call = "dGeneName", j_call = "jGeneAllele")

    type11 <- c(junction = "nucleotide", junction_aa = "aminoAcid", 
               duplicate_count = "count", 
               v_call = "vFamilyName", d_call = "dGeneName", j_call = "jGeneAllele")

    type_hash <- list("type1"=type1, "type2"=type2, "type3"=type3, "type4"=type4,
                      "type5"=type5, "type6"=type6, "type7"=type7, "type8"=type8,
                      "type9"=type9, "type10"=type10, "type11"=type11)
    return(type_hash[[file_type]])
}

#' @section Read clone file from path:
#' 
#' Given the path to a single AIRRSeq clone file, generate MiAIRR compliant tibble.
#'
#' @param clone_file .tsv file containing results from AIRRSeq pipeline
#'
#' @return Tibble in MiAIRR format 
#'
#' @rdname readImmunoSeq
readFiles <- function(clone_file, progress_bar) {
    progress_bar$tick()
    file_info <-getFileType(clone_file)
    file_type <- file_info[[1]]
    header_list <- file_info[[2]]
    col_std <- getStandard(file_type)
    col_old = header_list
    file_names <- tools::file_path_sans_ext(basename(clone_file))
    options(readr.show_progress = FALSE)
    clone_frame <- readr::read_tsv(clone_file, 
                                   col_types = col_old,
                                   na = c("", "NA", "Nan", "NaN", " "), 
                                   trim_ws = FALSE)
    clone_frame <- clone_frame %>% 
                   dplyr::rename(!!!col_std)
    clone_frame <- clone_frame %>% 
                   dplyr::mutate(reading_frame = dplyr::if_else(stringr::str_detect(junction_aa, "\\*") | is.na(junction_aa), 
                                                                "out-of-frame", 
                                                                "in-frame"),
                                 junction = dplyr::if_else(stringr::str_detect(junction, "[acgt]+"), 
                                                           toupper(stringr::str_extract(junction, "[acgt]+")), 
                                                           junction),
                                 junction_aa = dplyr::if_else(stringr::str_detect(junction_aa, "[acdefghiklmnpqrstvwy]+"),
                                                              toupper(stringr::str_extract(junction_aa, 
                                                                                           "[acdefghiklmnpqrstvwy]+")), 
                                                              junction_aa))
    clone_frame <- clone_frame %>%
                   dplyr::mutate(v_call = dplyr::if_else(stringr::str_detect(v_call, 
                                                                             "(TRB|TCRB|IGH|IGL|IGK|TCRA)V\\d+-\\d+\\*\\d+"),
                                                         stringr::str_extract(v_call, 
                                                                              "(TRB|TCRB|IGH|IGL|IGK|TCRA)V\\d+-\\d+\\*\\d+"), 
                                                         v_call),
                                 j_call = dplyr::if_else(stringr::str_detect(j_call, "(TRB|TCRB|IGH|IGL|IGK|TCRA)J\\d+-\\d+\\*\\d+"),
                                                         stringr::str_extract(j_call, "(TRB|TCRB|IGH|IGL|IGK|TCRA)J\\d+-\\d+\\*\\d+"), 
                                                         j_call),
                                 d_call = dplyr::if_else(stringr::str_detect(d_call, "(TRB|TCRB|IGH|IGL|IGK|TCRA)D\\d+-?\\d+\\*\\d+"),
                                                         stringr::str_extract(d_call, "(TRB|TCRB|IGH|IGL|IGK|TCRA)D\\d+-?\\d+\\*\\d+"), 
                                                         d_call))
    clone_frame <- clone_frame %>%
                   dplyr::mutate(v_family = dplyr::if_else(stringr::str_detect(v_call, "(TRB|TCRB|IGH|IGL|IGK|TCRA)V"),
                                                           stringr::str_extract(v_call, "(TRB|TCRB|IGH|IGL|IGK|TCRA)V\\d+"), 
                                                          "unresolved"),
                                 j_family = dplyr::if_else(stringr::str_detect(j_call, "(TRB|TCRB|IGH|IGL|IGK|TCRA)J"),
                                                           stringr::str_extract(j_call, "(TRB|TCRB|IGH|IGL|IGK|TCRA)J\\d+"), 
                                                           "unresolved"),
                                 d_family = dplyr::if_else(stringr::str_detect(d_call, "(TRB|TCRB|IGH|IGL|IGK|TCRA)D"),
                                                           stringr::str_extract(d_call, "(TRB|TCRB|IGH|IGL|IGK|TCRA)D\\d+"), 
                                                           "unresolved"))
    rep_col <- c(repertoire_id = file_names)
    clone_frame <- clone_frame %>% 
                   dplyr::group_by(junction, junction_aa, v_call, j_call, d_call, v_family, j_family, d_family, reading_frame) %>% 
                   dplyr::summarize(duplicate_count = sum(duplicate_count)) %>% 
                   dplyr::ungroup() %>%
                   tibble::add_column(!!!rep_col[!names(rep_col) %in% names(.)]) %>%
                   dplyr::mutate(duplicate_frequency = duplicate_count/sum(duplicate_count)) %>%
                   dplyr::select(repertoire_id, everything())
                 
    return(clone_frame)
}



