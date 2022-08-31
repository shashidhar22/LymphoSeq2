
#' Read 10x Genomics files
#' 
#' Imports tab-separated value (.tsv) files exported by 10x Genomics
#' and stores them as MiAIRR compliant tibbles.
#' 
#' @param path Path to the directory containing .tsv files. Only
#' files with the .tsv extension are imported.
#' @return Returns a tibble with MiAIRR headers and repertoire_id
#' @export
#' @import magrittr
read10x <- function(path, recursive = FALSE) {
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
    if(file_num != length(file_paths)) {
        warning("One or more of the files you are trying to import has no sequences and will be ignored.", 
                call. = FALSE)
    }
    file_list <- file_paths %>%
                 purrr::map(~ standardize10x(.x)) %>%
                 dplyr::bind_rows()
    return(file_list)
}
#' Read 10x Genomics data containing alpha and beta chains
#' 
#' Read in 10x Genomics data and collapse alpha and beta
#' chains appropriately.
#' 
#' @param clone_file A single .txv file from 10x Genomics
#' @return A tibble with alpha and beta chains collapsed
#' @export
collpase_chains <- function(clone_file) {
    clone_data <- readr::read_tsv(clone_file, na = c("", "NA", "Nan", "NaN", "unresolved"))

    double_fields <- c("v_sequence_start", "v_sequence_end",
                        "d_sequence_start", "d_sequence_end",
                        "j_sequence_start", "j_sequence_end",
                        "c_sequence_start", "c_sequence_end",
                        "junction_length", "junction_aa_length",
                        "consensus_count", "duplicate_count")
    collapsed_data <- clone_data %>%
                        dplyr::group_by(cell_id) %>%
                        dplyr::group_split() %>%
                        purrr::map(~mostPrevalent(., clone_data)) %>%
                        dplyr::bind_rows() %>%
                        dtplyr::lazy_dt()
    collapsed_data <- collapsed_data %>%
                        dplyr::group_by(cell_id, clone_id) %>%
                        dplyr::summarise_all(~paste(., collapse = ".")) %>%
                        dplyr::as_tibble() %>%
                        dplyr::mutate_at(double_fields,
                                         function(x) as.double(x))

    return(collapsed_data)
}

#' Select the most frequent chain
#' 
#' [description]
#' 
#' @param barcode_data A tibble that holds data for one barcode
#' identifier
#' @param clone_data A tibble that holds data for one repertoire_id
#' (one file)
#' @param chain The chain to examine to select the most frequently
#' occurring one. Values given can only be "TRA" or "TRB" to indicate
#' alpha or beta chain respectively.
#' @return A tibble with one row of data that contains the most
#' frequently occuring chain.
selectChain <- function(barcode_data, clone_data, chain = "TRA") {
    if (chain == "TRA") {
        opp_chain <- "TRB"
    } else {
        opp_chain <- "TRA"
    }
    vdj_TRx <- barcode_data %>%
                dtplyr::lazy_dt() %>%
                dplyr::filter(grepl(chain, v_call)) %>%
                dplyr::select(c(junction, v_call, j_call, d_call)) %>%
                dplyr::as_tibble()

    if (nrow(vdj_TRx) > 1) {
        TRx_count <- purrr::map(1:nrow(vdj_TRx),
                                function(x) {
                                    clone_lazy_data <- clone_data %>%
                                                        dtplyr::lazy_dt()
                                    clone_lazy_data %>%
                                        dplyr::filter(junction == vdj_TRx$junction[x],
                                                    v_call == vdj_TRx$v_call[x],
                                                    d_call == vdj_TRx$d_call[x] | is.na(d_call),
                                                    j_call == vdj_TRx$j_call[x]) %>%
                                        dplyr::as_tibble() %>%
                                        nrow()})
        barcode_data <- barcode_data %>%
                            dtplyr::lazy_dt() %>%
                            dplyr::filter(grepl(opp_chain, v_call) |
                                            v_call == vdj_TRx$v_call[which.max(TRx_count)]) %>%
                            dplyr::as_tibble()
    }
    return(barcode_data)
}

#' Selecting the alpha and beta chains
#' 
#' Select the most frequently occurring alpha and beta chains
#' for each barcode
#' 
#' @param barcode_data A tibble that holds data for one barcode
#' identifier
#' @param clone_data A tibble that holds data for one repertoire_id
#' (one file)
#' @return A tibble containing the most frequent alpha chain and
#' the most frequent beta chain.
mostPrevalent <- function(barcode_data, clone_data) {
    ab_chains <- barcode_data %>%
                    dplyr::pull(v_call)
    if (length(ab_chains) > 2) {
        barcode_data <- selectChain(barcode_data, clone_data, chain = "TRA")
        barcode_data <- selectChain(barcode_data, clone_data, chain = "TRB")
    }
    return(barcode_data)
}

standardize10x <- function(clone_file) {
    airr_headers_path <- system.file("extdata", "AIRR_fields.csv", package = "LymphoSeq2")
    empty_airr_frame <- readr::read_csv(airr_headers_path, trim_ws = TRUE)
    file_info <- collpase_chains(clone_file)
    if (ncol(file_info) == 144) {
        return(file_info)
    }
    clone_frame <- dplyr::right_join(empty_airr_frame, file_info)
    clone_frame <- clone_frame %>%
                dplyr::mutate(sequence_id = row_number(),
                    junction_length = stringr::str_length(junction),
                    junction_aa_length = stringr::str_length(junction_aa),
                    rev_comp = FALSE,
                    stop_codon = dplyr::if_else(stringr::str_detect(sequence, "\\*") | stringr::str_detect(sequence_aa, "\\*") | 
                                            is.na(sequence) | is.na(sequence_aa), TRUE, FALSE),
                    productive = dplyr::if_else(stringr::str_detect(sequence, "\\*") | stringr::str_detect(sequence_aa, "\\*") | 
                                            is.na(sequence) | is.na(sequence_aa), FALSE, TRUE),
                    v_call = dplyr::if_else(stringr::str_detect(v_call, "/"), "unresolved", v_call),
                    j_call = dplyr::if_else(stringr::str_detect(j_call, "/"), "unresolved", j_call),
                    complete_vdj = dplyr::if_else(is.na(v_call) | is.na(d_call) | is.na(j_call) | 
                                            stringr::str_detect(v_call, "unresolved") | stringr::str_detect(d_call, "unresolved") |
                                            stringr::str_detect(j_call, "unresolved"), FALSE, TRUE),
                    duplicate_frequency = duplicate_count / sum(duplicate_count),
                    reading_frame = dplyr::if_else(stringr::str_detect(sequence_aa, "\\*"), "out-of-frame", "in-frame"),
                    v_family = dplyr::if_else(stringr::str_detect(v_call, "(TRB|TCRB)V"),
                                            stringr::str_extract(v_call, "(TRB|TCRB)V\\d+"), "unrecognized"),
                    j_family = dplyr::if_else(stringr::str_detect(j_call, "(TRB|TCRB)J"),
                                            stringr::str_extract(j_call, "(TRB|TCRB)J\\d+"), "unrecognized"),
                    d_family = dplyr::if_else(stringr::str_detect(d_call, "(TRB|TCRB)D"),
                                            stringr::str_extract(d_call, "(TRB|TCRB)D\\d+"), "unrecognized"))

    return(clone_frame)
}