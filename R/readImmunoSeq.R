#' Read ImmunoSeq files
#'
#' Imports tab-separated value (.tsv) files exported by the Adaptive 
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
#' @md
#' @name readImmunoSeq
NULL
#' @describeIn readImmunoSeq Read a set of files 
#' @param path Path to the directory containing tab-delimited files.  Only 
#' files with the extension .tsv are imported.  The names of the data frames are 
#' the same as names of the files.
#' 
#' @return Returns a tibble with MiAIRR headers and repertoire_id
#'
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' study_table <- readImmunoSeq(path = file.path, 
#'                              recursive = FALSE)
#' @export
#' @import tidyverse dplyr
#' @rdname readImmunoSeq
readImmunoSeq <- function(path, recursive = FALSE) {
    file_paths <- list.files(path, full.names = TRUE, all.files = FALSE, 
                             recursive = recursive, pattern = ".tsv|.txt|.tsv.gz", 
                             include.dirs = FALSE)
    file_num <- length(file_paths)
    file_info <- file.info(file_paths)
    file_paths <- rownames(file_info)[file_info$size > 0]
    if(file_num != length(file_paths)){
        warning("One or more of the files you are trying to import has no sequences and will be ignored.", 
                call. = FALSE)
    }
    progress <- dplyr::progress_estimated(length(file_paths))
    file_list <- file_paths %>% purrr::map(~ readFiles(.x, progress)) %>% 
                 bind_rows()
    progress$stop()
    return(file_list)
}

#' @section Identify file type:
#'
#' AIRRSeq results have various formats, this functions identifies six formats 
#' 
#' @param clone_file .tsv file containing results from AIRRSeq pipeline
#'
#' @return List containing file type and standardized header names 
#'
#' @import tidyverse
#' @export
#' @rdname readImmunoSeq
getFileType <- function(clone_file) {
    type1 <- c("nucleotide", "aminoAcid", "count (reads)", "frequencyCount (%)", "vGeneName", 
        "dGeneName", "jGeneName", "vFamilyName", "dFamilyName", "jFamilyName", "sequenceStatus", 
        "estimatedNumberGenomes")
        
    type2 <- c("nucleotide", "aminoAcid", "count (templates/reads)", "frequencyCount (%)", "vGeneName", 
        "dGeneName", "jGeneName", "vFamilyName", "dFamilyName", "jFamilyName", "sequenceStatus", 
        "estimatedNumberGenomes")
    
    type3 <- c("nucleotide", "aminoAcid", "count (templates)", "frequencyCount (%)", "vGeneName", 
        "dGeneName", "jGeneName", "vFamilyName", "dFamilyName", "jFamilyName", "sequenceStatus", 
        "estimatedNumberGenomes")

    type4 <- c("nucleotide", "aminoAcid", "count", "frequencyCount", "vGeneName", 
        "dGeneName", "jGeneName", "vFamilyName", "dFamilyName", "jFamilyName", "sequenceStatus", 
        "estimatedNumberGenomes")
    
    type5 <- c("nucleotide.CDR3.in.lowercase.", "aminoAcid.CDR3.in.lowercase.", "cloneCount", 
        "clonefrequency....", "vGene", "dGene", "jGene", "fuction")

    type6 <- c("nucleotide(CDR3 in lowercase)", "aminoAcid(CDR3 in lowercase)", "cloneCount",
        "clonefrequency (%)", "vGene", "dGene", "jGene", "fuction")
    
    type7 <- c("cloneCount", "cloneFraction", "allVHitsWithScore", "allDHitsWithScore",
        "allJHitsWithScore", "nSeqCDR3", "aaSeqCDR3")

    type8 <- c("duplicate_count", "v_call", "d_call", "j_call", "junction", "junction_aa")

    columns <- invisible(colnames(readr::read_tsv(clone_file, n_max=1, col_types = cols())))
    if (all(type1 %in% columns)) {
        file_type <- "type1"
        header_list <- readr::cols_only(nucleotide = "c", aminoAcid = "c", `count (reads)` = "i", 
                vGeneName = "c", dGeneName = "c", jGeneName = "c")
    } else if (all(type2 %in% columns)) {
        file_type <- "type2"
        header_list <- readr::cols_only(nucleotide = "c", aminoAcid = "c", `count (templates/reads)` = "i", 
                vGeneName = "c", dGeneName = "c", jGeneName = "c")
    } else if (all(type3 %in% columns)) {
        file_type <- "type3"
        header_list <- readr::cols_only(nucleotide = "c", aminoAcid = "c", `count (templates)` = "i", 
                vGeneName = "c", dGeneName = "c", jGeneName = "c")
    } else if (all(type4 %in% columns)) {
        file_type <- "type4"
        header_list <- readr::cols_only(nucleotide = "c", aminoAcid = "c", `count` = "i", 
                vGeneName = "c", dGeneName = "c", jGeneName = "c")
    } else if (all(type5 %in% columns)) {
        file_type <- "type5"
        header_list <- readr::cols_only(`nucleotide.CDR3.in.lowercase.` = "c", 
                `aminoAcid.CDR3.in.lowercase.` = "c", cloneCount = "i", `clonefrequency....` = "d",
                vGene = "c", dGene = "c", jGene = "c")
    } else if (all(type6 %in% columns)) {
        file_type <- "type6"
        header_list <- readr::cols_only(`nucleotide(CDR3 in lowercase)` = "c", 
                `aminoAcid(CDR3 in lowercase)` = "c", cloneCount = "i", 
                `clonefrequency (%)` = "d", vGene = "c", dGene = "c", jGene = "c")
    } else if (all(type7 %in% columns)) {
        file_type <- "type7"
        header_list <- readr::cols_only(cloneCount = 'd', cloneFraction = 'd', 
                allVHitsWithScore = 'c', allDHitsWithScore = 'c', 
                allJHitsWithScore = 'c', nSeqCDR3 = 'c', aaSeqCDR3 = 'c')
    } else if (all(type8 %in% columns)) {
        file_type <- "type8"
        header_list <- readr::cols_only(duplicate_count = 'd', 
                v_call = 'c', d_call = 'c', j_call = 'c', junction = 'c', junction_aa = 'c')
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
#' @import tidyverse
#' @export
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
               v_call = "vGeneName", d_call = "dGeneName", j_call = "jGeneName")

    type4 <- c(junction = "nucleotide", junction_aa = "aminoAcid", 
               duplicate_count = "count", v_call = "vGeneName", d_call = "dGeneName", 
               j_call = "jGeneName")

    type5 <- c(junction = "nucleotide.CDR3.in.lowercase.", junction_aa = "aminoAcid.CDR3.in.lowercase.", 
               duplicate_count = "cloneCount", v_call = "vGene", d_call = "dGene", j_call = "jGene")   
    
    type6 <- c(junction = "nucleotide(CDR3 in lowercase)", junction_aa = "aminoAcid(CDR3 in lowercase)", 
               duplicate_count = "cloneCount", v_call = "vGene", d_call = "dGene", j_call = "jGene")   
    
    type7 <- c(junction = "nSeqCDR3", junction_aa = "aaSeqCDR3", duplicate_count = "cloneCount",
               v_call = "allVHitsWithScore", j_call = "allJHitsWithScore", d_call = "allDHitsWithScore")

    type8 <- c(junction = "junction", junction_aa = "junction_aa", duplicate_count = "duplicate_count",
               v_call = "v_call", j_call = "j_call", d_call = "d_call")

    type_hash <- list("type1"=type1, "type2"=type2, "type3"=type3, "type4"=type4,
                      "type5"=type5, "type6"=type6, "type7"=type7, "type8"=type8)
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
#' @import tidyverse
#' @export
#' @rdname readImmunoSeq
readFiles <- function(clone_file, progress) {
    progress$tick()$print()
    file_info <- getFileType(clone_file)
    file_type <- file_info[[1]]
    header_list <- file_info[[2]]
    col_std <- getStandard(file_type)
    col_old = header_list
    file_names <- tools::file_path_sans_ext(basename(clone_file))
    options(readr.show_progress = FALSE)
    clone_frame <- readr::read_tsv(clone_file, col_types = col_old,
                                    na = c("", "NA", "Nan", "NaN"), trim_ws = TRUE)
    clone_frame <- clone_frame %>% 
                   dplyr::rename(!!!col_std)
    clone_frame <- clone_frame %>% 
                   dplyr::mutate(reading_frame = dplyr::if_else(stringr::str_detect(junction_aa, "\\*"), 
                                                                "out-of-frame", "in-frame"),
                                 junction = dplyr::if_else(stringr::str_detect(junction, "[acgt]+"), 
                                             toupper(stringr::str_extract(junction, "[acgt]+")), junction),
                                 junction_aa = dplyr::if_else(stringr::str_detect(junction_aa, "[acdefghiklmnpqrstvwy]+"),
                                                toupper(stringr::str_extract(junction_aa, "[acdefghiklmnpqrstvwy]+")), junction_aa))
    clone_frame <- clone_frame %>%
                   dplyr::mutate(v_call = dplyr::if_else(stringr::str_detect(v_call, "(TRB|TCRB)V\\d+-\\d+\\*\\d+"),
                                           stringr::str_extract(v_call, "(TRB|TCRB)V\\d+-\\d+\\*\\d+"), v_call),
                                 j_call = dplyr::if_else(stringr::str_detect(j_call, "(TRB|TCRB)J\\d+-\\d+\\*\\d+"),
                                           stringr::str_extract(j_call, "(TRB|TCRB)J\\d+-\\d+\\*\\d+"), j_call),
                                 d_call = dplyr::if_else(stringr::str_detect(d_call, "(TRB|TCRB)D\\d+\\*\\d+"),
                                           stringr::str_extract(d_call, "(TRB|TCRB)D\\d+\\*\\d+"), d_call))
    clone_frame <- clone_frame %>%
                   dplyr::mutate(v_family = dplyr::if_else(stringr::str_detect(v_call, "(TRB|TCRB)V"),
                                             stringr::str_extract(v_call, "(TRB|TCRB)V\\d+"), "unrecognized"),
                                 j_family = dplyr::if_else(stringr::str_detect(j_call, "(TRB|TCRB)J"),
                                             stringr::str_extract(j_call, "(TRB|TCRB)J\\d+"), "unrecognized"),
                                 d_family = dplyr::if_else(stringr::str_detect(d_call, "(TRB|TCRB)D"),
                                             stringr::str_extract(d_call, "(TRB|TCRB)D\\d+"), "unrecognized"))
    clone_frame <- clone_frame %>% 
                   dplyr::group_by(junction, junction_aa, v_call, j_call, d_call, v_family, j_family, d_family, reading_frame) %>% 
                   dplyr::summarize(duplicate_count = sum(duplicate_count)) %>% 
                   dplyr::mutate(duplicate_frequency = duplicate_count / sum(duplicate_count)) %>% 
                   dplyr::ungroup() %>%
                   dplyr::mutate(repertoire_id = file_names) %>%
                   dplyr::select(repertoire_id, everything())
    return(clone_frame)
}



