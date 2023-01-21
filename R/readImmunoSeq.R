#' Read ImmunoSeq files
#'
#' `readImmunoSeq()` Imports tab-separated value (.tsv) files exported by the
#' Adaptive Biotechnologies ImmunoSEQ analyzer, BGI IR-SEQ, MiXCR and stores
#' them as MiAIRR compliant tibble.
#'
#' @param path Path to the directory containing tab-delimited files. Only files
#' with the extension .tsv are imported. The names of the data frames are
#' the same as names of the files.
#' @param recursive If TRUE, the function will recursively search the input
#' directory for all .tsv files
#' @param threads Number of threads.
#' @return Returns a tibble with MiAIRR headers and repertoire_id
#' @examples
#' file_path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq2")
#' study_table <- LymphoSeq2::readImmunoSeq(path = file_path, recursive = FALSE)
#' @export
#' @import magrittr
readImmunoSeq <- function(path, recursive = FALSE, threads = parallel::detectCores() / 2) {
  Sys.setenv("VROOM_SHOW_PROGRESS" = "false")
  if (length(path) > 1) {
    file_paths <- path
  } else if (file_test("-d", path)) {
    file_paths <- list.files(path,
      full.names = TRUE,
      all.files = FALSE,
      recursive = recursive,
      pattern = ".tsv|.txt|.tsv.gz",
      include.dirs = FALSE
    )
  } else {
    file_paths <- c(path)
  }
  file_num <- length(file_paths)
  file_info <- file.info(file_paths)
  file_paths <- rownames(file_info)[file_info$size > 0]
  if (file_num != length(file_paths)) {
    warning("One or more of the files you are trying to import has no sequences and will be ignored.",
      call. = FALSE
    )
  }

  progress_bar <- progress::progress_bar$new(
    format = "Reading AIRR-Seq files [:bar] :current/:total (:percent) eta: :eta elapsed: :elapsed",
    total = length(file_paths), clear = FALSE, width = 100
  )
  progress_bar$tick(0)
  file_list <- file_paths %>%
    purrr::map(~ getStandard(.x, progress_bar, threads = threads)) %>%
    dplyr::bind_rows()
  progress_bar$terminate()
  return(file_list)
}

#' `getFileType()` retrieve the file type of the input TSV file
#' @keywords internal
#' @param clone_file A .tsv file to identify the file type
#' @return Returns "immunoSEQLegacy", "immunoSEQ", "10X", "BGI"
#' @import magrittr
#' @noRd
getFileType <- function(col_names) {
  colname_file <- system.file("extdata", "Accepted_file_types.csv", package = "LymphoSeq2")
  colname_table <- vroom::vroom(colname_file, show_col_types = FALSE)
  immunoseq <- colname_table %>%
    dplyr::pull(immunoseq_v3) %>%
    purrr::discard(is.na)
  tenx <- colname_table %>%
    dplyr::pull(tenx) %>%
    purrr::discard(is.na)
  bgi <- colname_table %>%
    dplyr::pull(bgi) %>%
    purrr::discard(is.na)

  if (all(col_names %in% immunoseq)) {
    file_type <- "immunoSEQ"
  } else if (all(col_names %in% tenx)) {
    file_type <- "10X"
  } else if (all(col_names %in% bgi)) {
    file_type <- "BGI"
  } else {
    file_type <- "immunoSEQLegacy"
  }
  return(file_type)
}

#' `getStandard()` Converts AIRR-Seq data into MiAIRR compatible format
#' @keywords internal
#' @param clone_file A .tsv file to read in and standardize its fields to be MiAIRR compliant
#' @param airr_fields A character vector of MiAIRR headers
#' @return Tibble of given data with MiAIRR fields
#'
#' @import magrittr
#' @noRd
getStandard <- function(clone_file, progress, threads) {
  progress$tick()
  airr_headers_path <- system.file("extdata", "AIRR_fields.csv", package = "LymphoSeq2")
  airr_fields <- vroom::vroom(airr_headers_path,
    trim_ws = TRUE,
    show_col_types = FALSE
  )
  matching_fields <- getAIRRFields(clone_file, threads)
  col_read <- names(matching_fields)
  clone_data <- suppressWarnings(
    classes = c("vroom_mismatched_column_name", "vroom_parse_issue"),
    vroom::vroom(clone_file,
      col_types = c(
        vFamilyTies = readr::col_character(),
        jFamilyTies = readr::col_character()
      ),
      na = c("", "NA", "Nan", "NaN", "unresolved"),
      show_col_types = FALSE, progress = FALSE, num_threads = threads,
      .name_repair = ~ stringr::str_replace_all(., matching_fields)
    )
  )
  existing_match <- airr_fields %>%
    colnames() %>%
    intersect(colnames(clone_data))
  if (length(existing_match) == 155) {
    return(clone_data)
  }
  existing_airr_data <- clone_data %>%
    dplyr::select(all_of(existing_match))
  clone_data <- dplyr::bind_rows(airr_fields, existing_airr_data) %>%
    dplyr::slice(-1)
  file_name <- tools::file_path_sans_ext(basename(clone_file))
  clone_data <- clone_data %>%
    dplyr::mutate(
      repertoire_id = file_name,
      d2_call = dplyr::if_else(!is.na(d2_call), stringr::str_split(d2_call, ",")[[1]][2], d2_call),
      cdr3 = stringr::str_sub(junction, 4L, -4L),
      cdr3_aa = stringr::str_sub(junction_aa, 2L, -2L),
      cdr1_end = dplyr::if_else(is.na(cdr1_end), cdr1_end, cdr1_end + cdr1_start),
      cdr2_end = dplyr::if_else(is.na(cdr2_end), cdr2_end, cdr2_end + cdr2_start),
      cdr3_end = dplyr::if_else(is.na(cdr3_end), cdr3_end, cdr3_end + cdr3_start),
      sequence_id = dplyr::row_number(),
      sequence_aa = dplyr::if_else(stringr::str_detect(sequence_aa, "x"), stringr::str_remove_all(sequence_aa, "x"), sequence_aa),
      junction = dplyr::if_else(is.na(junction) & !is.na(sequence), sequence, junction),
      junction_aa = dplyr::if_else(is.na(junction_aa) & !is.na(sequence_aa), sequence_aa, junction_aa),
      junction = dplyr::if_else(stringr::str_detect(junction, "[a-z]+"), toupper(stringr::str_extract(junction, "[a-z]{2,}")), junction),
      junction_aa = dplyr::if_else(stringr::str_detect(junction_aa, "[a-z]+"), toupper(stringr::str_extract(junction_aa, "[a-z]{2,}")), junction_aa),
      junction_length = stringr::str_length(junction),
      junction_aa_length = stringr::str_length(junction_aa),
      rev_comp = FALSE,
      stop_codon = dplyr::if_else(stringr::str_detect(sequence, "\\*") | stringr::str_detect(sequence_aa, "\\*") |
        is.na(sequence) | is.na(sequence_aa), TRUE, FALSE),
      productive = stop_codon,
      v_call = stringr::str_remove(v_call, "/\\w+$"),
      j_call = stringr::str_remove(j_call, "/\\w+$"),
      complete_vdj = dplyr::if_else(is.na(v_call) | is.na(d_call) | is.na(j_call), FALSE, TRUE),
      duplicate_frequency = duplicate_count / sum(duplicate_count),
      reading_frame = dplyr::if_else(stop_codon, "out-of-frame", "in-frame"),
      v_call = stringr::str_c(stringr::str_extract(v_call, "[A-Z]+"),
        as.numeric(stringr::str_extract(v_call, "\\d+")),
        as.numeric(stringr::str_extract(v_call, "-\\d+")),
        sep = ""
      ),
      j_call = stringr::str_c(stringr::str_extract(j_call, "[A-Z]+"),
        as.numeric(stringr::str_extract(j_call, "\\d+")),
        as.numeric(stringr::str_extract(j_call, "-\\d+")),
        sep = ""
      ),
      d_call = stringr::str_c(stringr::str_extract(j_call, "[A-Z]+"),
        as.numeric(stringr::str_extract(j_call, "\\d+")),
        as.numeric(stringr::str_extract(j_call, "-\\d+")),
        sep = ""
      ),
      v_call = stringr::str_replace(v_call, "TCR", "TR"),
      j_call = stringr::str_replace(j_call, "TCR", "TR"),
      d_call = stringr::str_replace(d_call, "TCR", "TR"),
      j_family = stringr::str_extract(j_call, "[A-Z]+\\d+"),
      v_family = stringr::str_extract(v_call, "[A-Z]+\\d+"),
      d_family = stringr::str_extract(d_call, "[A-Z]+\\d+"),
      bio_identity = stringr::str_c(junction_aa, v_call, j_call, sep = "_"),
      sequence_id = stringr::str_c(repertoire_id, dplyr::row_number(), sep = "_"),
      clone_id = bio_identity
    )
  return(clone_data)
}

#' `getAIRRFields()` Given the path to a single AIRRSeq clone file, determine
#' the file type and returns a named vector that can be used to repair headers
#' while reading input.
#' @keywords internal
#' @param clone_file .tsv file containing results from AIRRSeq pipeline
#' @return Named vector of corresponding AIRR fields
#'
#' @import magrittr
#' @noRd
getAIRRFields <- function(clone_file, threads) {
  clone_table <- vroom::vroom(clone_file,
    show_col_types = FALSE,
    n_max = 1, num_threads = threads
  )
  col_names <- base::colnames(clone_table)
  input_type <- getFileType(col_names)
  if (input_type == "immunoSEQ") {
    count_method <- clone_table %>%
      dplyr::pull(counting_method) %>%
      base::unique()
    ## NOTE: These following fields need to be modified when the data is read
    ## 1. d2_call = When it contains a list with the first values being equal to d_call, extract the second element
    ## 2. cdr3 = Remove the first and last codon
    ## 3. cdr3_aa = Remove the first and last amino acid
    ## 4. bioidentity = Reformat in IMGT format
    ## 5. v_call, d_call, d2_call, j_call = Reformat in IMGT format
    ## 6. cd*_end = Add value to cdr*start -3
    ## 7. cd*_start = Subtract 3 from value
    ## 8. junction_aa_length = Divide value by 3
    matching_fields <- c(
      "bio_indetity" = "sequence_id",
      "rearrangement" = "sequence", "amino_acid" = "sequence_aa",
      "frame_type" = "productive", "v_gene" = "v_call",
      "d_gene" = "d_call", "d_gene_ties" = "d2_call", "j_gene" = "j_call",
      "cdr3_rearrangement" = "junction", "cdr3_amino_acid" = "junction_aa",
      "cdr3_sequence" = "junction", "cdr3_sequence_aa" = "junction_aa",
      "cdr1_rearrangement" = "cdr1", "cdr1_amino_acid" = "cdr1_aa",
      "cdr2_rearrangement" = "cdr2", "cdr2_amino_acid" = "cdr2_aa",
      "cdr1_sequence" = "cdr1", "cdr1_sequence_aa" = "cdr1_aa",
      "cdr2_sequence" = "cdr2", "cdr2_sequence_aa" = "cdr2_aa",
      "cdr1_start_index" = "cdr1_start",
      "cdr1_rearrangement_length" = "cdr1_end",
      "cdr2_start_index" = "cdr2_start",
      "cdr2_rearrangement_length" = "cdr2_end",
      "cdr3_start_index" = "cdr3_start",
      "cdr3_rearrangement_length" = "cdr3_end",
      "cdr3_length" = "junction_length",
      "cdr3_length" = "junction_aa_length",
      "n1_insertions" = "n1_length",
      "n2_insertions" = "n2_length"
    )

    if (count_method %in% c("v1")) {
      matching_fields <- c(matching_fields, "seq_reads" = "duplicate_count")
    } else {
      matching_fields <- c(matching_fields, "templates" = "duplicate_count")
    }
  } else if (input_type == "immunoSEQLegacy" | input_type == "BGI") {
    ## NOTE: These following fields need to be modified when the data is read
    ## 1. if input_type == "immunoSEQLegacy", set sequence_aa = junction_aa
    ## 2. if input_type == "immunoSEQLegacy", set sequence = junction
    matching_fields <- c(
      "amino_acid" = "sequence_aa",
      "\\AaminoAcid\\z" = "sequence_aa",
      "\\AaminoAcid.CDR3.in.lowercase.\\z" = "sequence_aa",
      "\\AaminoAcid\\(CDR3 in lowercase\\)\\z" = "sequence_aa",
      "\\ACDR3.stripped.x.a\\z" = "sequence_aa",
      "\\Anucleotide\\z" = "sequence",
      "\\Anucleotide.CDR3.in.lowercase.\\z" = "sequence",
      "\\Anucleotide\\(CDR3 in lowercase\\)\\z" = "sequence",
      "\\Acount \\(templates/reads\\)\\z" = "duplicate_count",
      "\\Acount \\(templates\\)\\z" = "duplicate_count",
      "\\Acount \\(reads\\)\\z" = "duplicate_count",
      "\\Acount\\z" = "duplicate_count",
      "\\AcloneCount\\z" = "duplicate_count",
      "\\Atemplates\\z" = "duplicate_count",
      "frame_type" = "productive", "fuction" = "productive",
      "locus" = "locus", "dGeneName" = "d_call",
      "dGeneNameTies" = "d2_call", "jGeneName" = "j_call",
      "vGeneName" = "v_call", "\\AvGene\\z" = "v_call",
      "\\AdGene\\z" = "d_call", "\\AjGene\\z" = "j_call"
    )
    count_cols <- c(
      "count (template/reads)", "count (templates)",
      "count (reads)", "count", "templates"
    )
    if (length(intersect(col_names, count_cols)) == 0) {
      matching_fields <- c(matching_fields, "estimatedNumberGenomes" = "duplicate_count")
    }
  } else if (input_type == "10X") {
    matching_fields <- col_names
    names(matching_fields) <- col_names
  }
  return(matching_fields)
}
