#' Import 10x Genomics single-cell VDJ data
#'
#' Load TCR or BCR sequencing data from 10x Genomics Cell Ranger VDJ output.
#' Automatically detects file format (AIRR standard or contig_annotations.csv)
#' and converts to MiAIRR-compliant format for use with LymphoSeq2.
#'
#' @param path Path to directory containing 10x files, or a vector of file paths.
#' Accepts .tsv, .csv, .txt, or .tsv.gz formats.
#' @param recursive Logical. Should subdirectories be searched for files?
#'  * `TRUE`: Recursively search all subdirectories
#'  * `FALSE` (default): Only search the specified directory
#'
#' @return A tibble with MiAIRR-formatted columns including:
#' * `cell_id`: Cell barcode identifier
#' * `repertoire_id`: Sample name (from filename)
#' * `junction`, `junction_aa`: CDR3 nucleotide and amino acid sequences
#' * `v_call`, `d_call`, `j_call`: V(D)J gene assignments
#' * `duplicate_count`: UMI count for the sequence
#' * `productive`: Whether the sequence is in-frame and functional
#'
#' @details
#' Supported file formats:
#' AIRR format (*_airr.tsv) - Standard AIRR rearrangement TSV from Cell Ranger.
#' Contig annotations (*_contig_annotations.csv) - Legacy Cell Ranger CSV format.
#' Legacy 10x - Older TSV formats (automatically detected and converted).
#'
#' Usage with single-cell data:
#' After importing, use merge_chains() to pair alpha/beta chains from the same cell.
#' For example: sc_data <- read10x("path/to/10x_output/") followed by
#' paired <- merge_chains(sc_data, mode = "strict").
#'
#' File naming: The repertoire_id is extracted from the filename. For example,
#' Sample1_airr.tsv becomes repertoire_id = "Sample1".
#' @export
read10x <- function(path, recursive = FALSE) {
  if (length(path) > 1) {
    file_paths <- path
  } else if (utils::file_test("-d", path)) {
    file_paths <- list.files(path,
      full.names = TRUE,
      all.files = FALSE,
      recursive = recursive,
      pattern = "\\.tsv$|\\.txt$|\\.csv$|\\.tsv\\.gz$",
      include.dirs = FALSE
    )
  } else {
    file_paths <- c(path)
  }

  # Filter empty files
  file_info <- file.info(file_paths)
  non_empty_files <- rownames(file_info)[file_info$size > 0]

  if (length(non_empty_files) == 0) {
    stop("No valid files found to import.", call. = FALSE)
  }

  if (length(file_paths) != length(non_empty_files)) {
    warning(
      "One or more files have no sequences and will be ignored.",
      call. = FALSE
    )
  }

  # Process files and combine
  file_list <- non_empty_files |>
    purrr::map(standardize10x) |>
    dplyr::bind_rows()

  return(file_list)
}

#' Detect 10x file format
#'
#' Determines if a file is AIRR-formatted or contig_annotations format
#'
#' @param file_path Path to the file
#' @return String: "airr", "contig_annotations", or "legacy"
detect_10x_format <- function(file_path) {
  # Check CSV files first
  if (grepl("\\.csv$", file_path, ignore.case = TRUE)) {
    headers <- readr::read_csv(file_path, n_max = 0, show_col_types = FALSE) |>
      names()
    if ("barcode" %in% headers && "contig_id" %in% headers) {
      return("contig_annotations")
    }
  }

  # Check TSV files
  headers <- readr::read_tsv(file_path, n_max = 0, show_col_types = FALSE) |>
    names()

  if (all(c("cell_id", "sequence_id", "junction", "v_call") %in% headers)) {
    return("airr")
  } else if ("barcode" %in% headers && "contig_id" %in% headers) {
    return("contig_annotations")
  } else {
    return("legacy")
  }
}

#' Read AIRR-formatted 10x file
#'
#' @param file_path Path to AIRR-formatted .tsv file
#' @return Tibble with AIRR data
read_airr_format <- function(file_path) {
  airr_data <- readr::read_tsv(file_path,
    na = c("", "NA", "Nan", "NaN", "unresolved", "None"),
    show_col_types = FALSE
  )

  # Extract repertoire_id from filename
  repertoire_id <- basename(file_path) |>
    stringr::str_remove("_airr\\.tsv$") |>
    stringr::str_remove("\\.tsv(\\.gz)?$")

  airr_data |>
    dplyr::mutate(repertoire_id = repertoire_id)
}

#' Read contig_annotations format
#'
#' @param file_path Path to contig_annotations.csv file
#' @return Tibble mapped to AIRR format
read_contig_annotations <- function(file_path) {
  contig_data <- readr::read_csv(file_path, show_col_types = FALSE)

  # Extract repertoire_id from filename
  repertoire_id <- basename(file_path) |>
    stringr::str_remove("_contig_annotations\\.csv$") |>
    stringr::str_remove("\\.csv$")

  # Map contig_annotations to AIRR fields
  contig_data |>
    dplyr::mutate(
      cell_id = barcode,
      sequence_id = contig_id,
      clone_id = raw_clonotype_id,
      productive = productive %in% c("true", TRUE),
      rev_comp = FALSE,
      v_call = v_gene,
      d_call = d_gene,
      j_call = j_gene,
      c_call = c_gene,
      junction = cdr3_nt,
      junction_aa = cdr3,
      junction_length = stringr::str_length(cdr3_nt),
      junction_aa_length = stringr::str_length(cdr3),
      duplicate_count = as.double(umis),
      consensus_count = as.double(reads),
      is_cell = is_cell %in% c("true", TRUE),
      repertoire_id = repertoire_id
    ) |>
    dplyr::select(
      cell_id, clone_id, sequence_id, productive, rev_comp,
      v_call, d_call, j_call, c_call, junction, junction_aa,
      junction_length, junction_aa_length, duplicate_count,
      consensus_count, is_cell, repertoire_id
    )
}

#' Standardize 10X data inputs
#'
#' @param clone_file A single file from 10x Genomics (AIRR .tsv or contig_annotations.csv)
#' @return A tibble standardized to AIRR format.
standardize10x <- function(clone_file) {
  format_type <- detect_10x_format(clone_file)

  switch(format_type,
    "airr" = read_airr_format(clone_file),
    "contig_annotations" = read_contig_annotations(clone_file),
    "legacy" = collapse_chains(clone_file)
  )
}

#' Read 10x Genomics data containing alpha and beta chains
#'
#' Read in 10x Genomics data and collapse alpha and beta
#' chains appropriately. This is for legacy 10x formats.
#'
#' @param clone_file A single .tsv file from 10x Genomics
#' @return A tibble with alpha and beta chains collapsed
#' @export
collapse_chains <- function(clone_file) {
  clone_data <- readr::read_tsv(clone_file,
    na = c("", "NA", "Nan", "NaN", "unresolved"),
    show_col_types = FALSE
  )

  # Use data.table for faster grouped operations
  collapsed_data <- clone_data |>
    dtplyr::lazy_dt() |>
    dplyr::group_by(cell_id) |>
    dplyr::group_modify(~ mostPrevalent(.x, clone_data)) |>
    dplyr::ungroup() |>
    dplyr::group_by(cell_id, clone_id) |>
    dplyr::summarise(
      dplyr::across(dplyr::everything(), ~ paste(., collapse = ".")),
      .groups = "drop"
    ) |>
    dplyr::as_tibble()

  # Convert numeric fields back to double
  numeric_fields <- c(
    "v_sequence_start", "v_sequence_end",
    "d_sequence_start", "d_sequence_end",
    "j_sequence_start", "j_sequence_end",
    "c_sequence_start", "c_sequence_end",
    "junction_length", "junction_aa_length",
    "consensus_count", "duplicate_count"
  )

  collapsed_data |>
    dplyr::mutate(
      dplyr::across(
        dplyr::any_of(numeric_fields),
        ~ as.double(.)
      )
    )
}

#' Select the most frequent chain
#'
#' @param barcode_data A tibble that holds data for one barcode identifier
#' @param clone_data A tibble that holds data for one repertoire_id (one file)
#' @param chain The chain to examine to select the most frequently
#' occurring one. Values can be "TRA" or "TRB" to indicate
#' alpha or beta chain respectively.
#' @return A tibble with one row of data that contains the most
#' frequently occurring chain.
selectChain <- function(barcode_data, clone_data, chain = "TRA") {
  opp_chain <- if (chain == "TRA") "TRB" else "TRA"

  vdj_TRx <- barcode_data |>
    dplyr::filter(grepl(chain, v_call)) |>
    dplyr::select(junction, v_call, j_call, d_call)

  if (nrow(vdj_TRx) > 1) {
    # Count occurrences of each chain variant in the entire dataset
    TRx_count <- purrr::map_int(
      seq_len(nrow(vdj_TRx)),
      function(x) {
        clone_data |>
          dplyr::filter(
            junction == vdj_TRx$junction[x],
            v_call == vdj_TRx$v_call[x],
            (is.na(d_call) & is.na(vdj_TRx$d_call[x])) | d_call == vdj_TRx$d_call[x],
            j_call == vdj_TRx$j_call[x]
          ) |>
          nrow()
      }
    )

    # Keep opposite chain and most frequent current chain
    most_freq_v_call <- vdj_TRx$v_call[which.max(TRx_count)]
    barcode_data <- barcode_data |>
      dplyr::filter(grepl(opp_chain, v_call) | v_call == most_freq_v_call)
  }

  barcode_data
}

#' Selecting the alpha and beta chains
#'
#' Select the most frequently occurring alpha and beta chains
#' for each barcode
#'
#' @param barcode_data A tibble that holds data for one barcode identifier
#' @param clone_data A tibble that holds data for one repertoire_id (one file)
#' @return A tibble containing the most frequent alpha chain and
#' the most frequent beta chain.
mostPrevalent <- function(barcode_data, clone_data) {
  n_chains <- nrow(barcode_data)

  if (n_chains > 2) {
    barcode_data <- barcode_data |>
      selectChain(clone_data, chain = "TRA") |>
      selectChain(clone_data, chain = "TRB")
  }

  barcode_data
}
