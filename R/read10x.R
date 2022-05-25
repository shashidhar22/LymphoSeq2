#' @param clone_data
#'
#' @return 
#'
#'
read10x <- function(clone_data) {
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
                        dplyr::mutate_at(double_fields, function(x) as.double(x))

    return(collapsed_data)
}

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

#' @param barcode_data
#' @param clone_data
#'
#' @return 
#'
#'
mostPrevalent <- function(barcode_data, clone_data) {
    ab_chains <- barcode_data %>%
                    dplyr::pull(v_call)
    if (length(ab_chains) > 2) {
        barcode_data <- selectChain(barcode_data, clone_data, chain = "TRA")
        barcode_data <- selectChain(barcode_data, clone_data, chain = "TRB")
    }
    return(barcode_data)
}