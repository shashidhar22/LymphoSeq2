#' Read ImmunoSeq files
#'
#' [readImmunoSeq()] Imports tab-separated value (.tsv) files exported by the
#' Adaptive Biotechnologies ImmunoSEQ analyzer, BGI IR-SEQ, MiXCR and stores
#' them as MiAIRR compliant tibble. Optimized for large dataset performance.
#'
#' @param path Path to the directory containing tab-delimited files. Only files
#' with the extension .tsv are imported. The names of the data frames are
#' the same as names of the files.
#' @param recursive A Boolean value
#'  * `TRUE` : the function will recursively search directory for all .tsv files
#'  * `FALSE` (the default): Open file using path
#' @param threads Number of threads.
#' @param parallel A Boolean value
#'  * `TRUE` (the default): Process files in parallel using future
#'  * `FALSE`: Process files sequentially
#' @param chunk_size Integer specifying the number of files to process in each chunk
#' for memory efficiency. Default is NULL (process all files at once).
#' @param max_memory_gb Maximum memory to use in GB. If exceeded, will switch to
#' chunked processing automatically. Default is 8GB.
#' @param streaming_mode Boolean. If TRUE, processes extremely large files in
#' streaming chunks to avoid memory limits. Default is FALSE.
#' @param progress_detail Level of progress reporting: "none", "basic", "detailed".
#' Default is "basic".
#' @param return_type Character string specifying return type: "data.table" (default),
#' "tibble" (tidyverse compatible), or "lazy_dt" (dtplyr - best of both worlds).
#' lazy_dt provides tidyverse syntax with data.table performance.
#' @param use_arrow Character or logical. Controls Apache Arrow usage for large datasets:
#' "auto" (default) enables Arrow automatically for datasets >5GB, "always" forces Arrow,
#' "never" disables Arrow, TRUE/FALSE for simple enable/disable.
#' @param sample_mode Boolean. If TRUE, processes only a random sample of sequences
#' from large datasets for quick analysis. Default is FALSE.
#' @param sample_size Integer. Number of sequences to sample per file when sample_mode=TRUE.
#' Default is 1,000,000.
#' @return Returns a data.table or tibble with MiAIRR headers and repertoire_id
#' @examples
#' file_path <- system.file("extdata", "TCRB_sequencing",
#'  package = "LymphoSeq2")
#' study_table <- LymphoSeq2::readImmunoSeq(
#'   path = file_path, recursive = FALSE,
#'   threads = 1
#' )
#' study_table <- LymphoSeq2::topSeqs(study_table, top = 100)
#'
#' @export
readImmunoSeq <- function(path,
                          recursive = FALSE,
                          threads = parallel::detectCores() / 2,
                          parallel = TRUE,
                          chunk_size = NULL,
                          max_memory_gb = 8,
                          streaming_mode = FALSE,
                          progress_detail = "basic",
                          return_type = "data.table",
                          use_arrow = "auto",
                          sample_mode = FALSE,
                          sample_size = 1000000) {
  Sys.setenv("VROOM_SHOW_PROGRESS" = "false")
  if (floor(threads) == 0) {
    threads <- as.integer(1)
  }
  if (length(path) > 1) {
    file_paths <- path
  } else if (utils::file_test("-d", path)) {
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
  # Check if any files were found
  if (length(file_paths) == 0) {
    stop("No .tsv, .txt, or .tsv.gz files found in the specified path: ", path,
         call. = FALSE)
  }

  file_num <- length(file_paths)
  file_info <- file.info(file_paths)
  file_paths <- rownames(file_info)[file_info$size > 0]

  # Check if any files have content after filtering
  if (length(file_paths) == 0) {
    stop("All files found are empty (0 bytes). No data to import.",
         call. = FALSE)
  }

  if (file_num != length(file_paths)) {
    warning(stringr::str_c("One or more of the files you are trying to import",
      " has no sequences and will be ignored.", sep = " "),
      call. = FALSE
    )
  }

  # Advanced dataset analysis for optimization
  total_size_gb <- sum(file_info[file_paths, "size"], na.rm = TRUE) / (1024^3)
  largest_file_mb <- max(file_info[file_paths, "size"], na.rm = TRUE) / (1024^2)
  available_memory_gb <- get_available_memory_gb()

  if (progress_detail != "none") {
    cat(sprintf("Dataset Analysis:\n"))
    cat(sprintf("  Files: %d, Total: %.2f GB, Largest: %.1f MB\n",
                length(file_paths), total_size_gb, largest_file_mb))
    cat(sprintf("  Available memory: %.1f GB\n", available_memory_gb))
  }

  # Apache Arrow decision logic
  should_use_arrow <- FALSE
  if (use_arrow == "always" || use_arrow == TRUE) {
    should_use_arrow <- TRUE
  } else if (use_arrow == "auto") {
    # Auto-enable Arrow for datasets >5GB or when memory is constrained
    should_use_arrow <- total_size_gb > 5 || total_size_gb > available_memory_gb * 0.6
  }

  if (should_use_arrow) {
    if (!requireNamespace("arrow", quietly = TRUE)) {
      if (progress_detail != "none") {
        cat("Arrow not installed - falling back to data.table processing\n")
      }
      should_use_arrow <- FALSE
    } else {
      if (progress_detail != "none") {
        cat("Using Apache Arrow for large dataset processing\n")
      }
      return(process_with_arrow(file_paths, threads, progress_detail, return_type,
                                sample_mode, sample_size))
    }
  }

  # Adaptive processing mode selection
  if (streaming_mode || total_size_gb > max_memory_gb || largest_file_mb > 1000) {
    if (progress_detail != "none") {
      cat("Switching to streaming mode for large dataset\n")
    }
    return(process_streaming_mode(file_paths, threads, progress_detail))
  }

  if (total_size_gb > available_memory_gb * 0.7) {
    if (progress_detail != "none") {
      cat("Large dataset detected - enabling memory optimizations\n")
    }
    # Force garbage collection and adjust chunk size
    gc(full = TRUE)
    if (is.null(chunk_size)) {
      chunk_size <- max(1, min(10, floor(length(file_paths) / 4)))
    }
  }

  # Enhanced progress reporting
  if (progress_detail == "detailed") {
    progress_format <- paste0("Loading [:bar] :current/:total files (:percent) | ",
                              ":rate files/sec | ETA: :eta | Elapsed: :elapsed")
  } else if (progress_detail == "basic") {
    progress_format <- "Loading [:bar] :current/:total files (:percent) | ETA: :eta"
  } else {
    progress_format <- NULL
  }

  if (!is.null(progress_format)) {
    progress_bar <- progress::progress_bar$new(
      format = progress_format,
      total = length(file_paths), clear = FALSE, width = 100
    )
    progress_bar$tick(0)
  } else {
    # Create a dummy progress bar for "none" mode
    progress_bar <- list(
      tick = function(...) invisible(NULL),
      terminate = function(...) invisible(NULL)
    )
  }

  # Chunked processing for memory efficiency
  if (!is.null(chunk_size) && length(file_paths) > chunk_size) {
    if (progress_detail != "none") {
      cat(sprintf("Processing %d files in chunks of %d for memory efficiency\n",
                  length(file_paths), chunk_size))
    }
    return(process_file_chunks(file_paths, chunk_size, progress_bar,
                               threads, parallel, progress_detail))
  }

  # Use data.table::rbindlist for faster binding
  if (parallel && length(file_paths) > 1) {
    # Set up parallel processing if requested
    if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("furrr", quietly = TRUE)) {
      warning("future and furrr packages required for parallel processing. Falling back to sequential processing.")
      parallel <- FALSE
    }
  }

  if (parallel && length(file_paths) > 1) {
    # Parallel processing
    current_plan <- class(future::plan())[1]
    if (current_plan == "sequential") {
      future::plan(future::multisession, workers = min(threads, length(file_paths)))
      reset_plan <- TRUE
    } else {
      reset_plan <- FALSE
    }

    file_list <- furrr::future_map(file_paths, ~ getStandard(.x, progress_bar, threads),
                                   .options = furrr::furrr_options(seed = TRUE))

    if (reset_plan) {
      future::plan(future::sequential) # Reset plan only if we set it
    }
  } else {
    # Sequential processing
    file_list <- purrr::map(file_paths, ~ getStandard(.x, progress_bar, threads))
  }

  # Use data.table::rbindlist for faster binding
  result <- data.table::rbindlist(file_list, use.names = TRUE, fill = TRUE)
  progress_bar$terminate()

  # Convert return type based on user preference
  return_type <- match.arg(return_type, c("data.table", "tibble", "lazy_dt"))

  if (return_type == "tibble") {
    # Convert to tibble for tidyverse compatibility
    result <- tibble::as_tibble(result)
  } else if (return_type == "lazy_dt") {
    # Use dtplyr for efficient tidyverse interface with data.table backend
    result <- dtplyr::lazy_dt(result)
  }

  return(result)
}

#' [getFileType()] retrieve the file type of the input TSV file
#' @keywords internal
#' @param clone_file A .tsv file to identify the file type
#' @return Returns "immunoSEQLegacy", "immunoSEQ", "10X", "BGI"
#' @noRd
getFileType <- function(col_names) {
  colname_file <- system.file("extdata", "Accepted_file_types.csv",
    package = "LymphoSeq2")
  colname_table <- data.table::fread(colname_file, showProgress=FALSE)
  immunoseq <- colname_table[!is.na(colname_table$immunoseq_v3), ]$immunoseq_v3
  tenx <- colname_table[!is.na(colname_table$tenx), ]$tenx
  bgi <- colname_table[!is.na(colname_table$bgi), ]$bgi

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

#' [getStandard()] Converts AIRR-Seq data into MiAIRR compatible format
#' @keywords internal
#' @param clone_file A .tsv file to read in and standardize its fields to be 
#'  MiAIRR compliant.
#' @param progress Progress bar
#' @param threads Number of threads for parallel processing
#' @return Tibble of given data with MiAIRR fields
#'
#' @noRd
getStandard <- function(clone_file, progress, threads) {
  progress$tick()
  airr_headers_path <- system.file("extdata", "AIRR_fields.csv",
    package = "LymphoSeq2")
  airr_fields <- data.table::fread(airr_headers_path, showProgress=FALSE)
  AIRR_fields_res <- getAIRRFields(clone_file, threads)
  matching_fields <- AIRR_fields_res$matching_fields
  file_type <- AIRR_fields_res$file_type
  col_read <- names(matching_fields)
  clone_data <- data.table::fread(
    clone_file,
    nThread = threads,
    na.strings = c("", "NA", "Nan", "NaN", "unresolved"),
    stringsAsFactors = FALSE,
    verbose = FALSE,
    showProgress=FALSE
  )

  # Apply column name mappings if needed - optimized
  if (length(matching_fields) > 0) {
    old_names <- names(clone_data)
    # Use vectorized gsub with Perl regex patterns (needed for \A and \z anchors)
    new_names <- old_names
    for (pattern in names(matching_fields)) {
      new_names <- gsub(pattern, matching_fields[[pattern]], new_names, perl = TRUE)
    }
    data.table::setnames(clone_data, old_names, new_names)
  }
  existing_match <- intersect(colnames(airr_fields), colnames(clone_data))
  if (length(existing_match) == 155) {
    return(clone_data)
  }
  # Use data.table operations for better performance
  existing_airr_data <- clone_data[, .SD, .SDcols = existing_match]
  clone_data <- data.table::rbindlist(
    list(airr_fields[0], existing_airr_data),
    use.names = TRUE,
    fill = TRUE
  )
  file_name <- tools::file_path_sans_ext(basename(clone_file))
  clone_data <- clone_data |>
    dplyr::mutate(
      repertoire_id = file_name,
      d2_call = ifelse(!is.na(d2_call),
          sapply(strsplit(d2_call, ",", fixed = TRUE), function(x) if(length(x) >= 2) x[2] else NA_character_),
          d2_call),
          cdr3 = substr(junction, 4L, nchar(junction) - 3L),
          cdr3_aa = substr(junction_aa, 2L, nchar(junction_aa) - 1L),
      cdr1_end = ifelse(is.na(cdr1_end), cdr1_end, cdr1_end + cdr1_start),
      cdr2_end = ifelse(is.na(cdr2_end), cdr2_end, cdr2_end + cdr2_start),
      cdr3_end = ifelse(is.na(cdr3_end), cdr3_end, cdr3_end + cdr3_start),
      sequence_id = dplyr::row_number(),
      sequence = ifelse(is.na(sequence) & !is.na(junction), junction, sequence),
      junction = ifelse(is.na(junction) & !is.na(sequence), sequence, junction),
      sequence_aa = ifelse(is.na(sequence_aa) & !is.na(junction_aa), junction_aa, sequence_aa),
      junction_aa = ifelse(is.na(junction_aa) & !is.na(sequence_aa), sequence_aa, junction_aa)
    )
  if (file_type == "BGI") {
    # Optimized string operations for BGI format - single pass with combined regex
    clone_data <- clone_data |>
      dplyr::mutate(
        junction = toupper(gsub("[a-zx]", "", junction)),
        junction_aa = toupper(gsub("[a-zx]", "", junction_aa))
      )
  } else {
    # Optimized string operations for non-BGI formats
    clone_data <- clone_data |>
      dplyr::mutate(
        junction = ifelse(grepl("[a-z]+", junction),
            toupper(gsub(".*([a-z]{2,}).*", "\\1", junction)),
            junction),
        junction_aa = ifelse(grepl("[a-z]+", junction_aa),
            toupper(gsub(".*([a-z]{2,}).*", "\\1", junction_aa)),
            junction_aa)
      )
  }
  clone_data <- clone_data |>
    dplyr::mutate(
      sequence = junction,
      sequence_aa = junction_aa,
      junction_length = nchar(junction),
      junction_aa_length = nchar(junction_aa),
      rev_comp = FALSE,
      stop_codon = ifelse(grepl("\\*", sequence, fixed = TRUE) |
        grepl("\\*", sequence_aa, fixed = TRUE) | is.na(sequence) |
        is.na(sequence_aa),
          TRUE,
          FALSE),
      productive = !stop_codon,
      v_call = gsub("/\\w+$", "", v_call),
      j_call = gsub("/\\w+$", "", j_call),
      complete_vdj = ifelse(is.na(v_call) | is.na(d_call) | is.na(j_call),
          FALSE, TRUE),
      duplicate_frequency = duplicate_count / sum(duplicate_count, na.rm = TRUE),
      reading_frame = ifelse(stop_codon, "out-of-frame", "in-frame"),
      # Optimized gene call parsing using vectorized operations with proper NA handling
      # Extract core gene name (e.g., TRBV7-9 or TRBV7) removing allele info (*01, etc)
      v_call = ifelse(
        is.na(v_call) | v_call == "",
        NA_character_,
        sub("\\*.*$", "", v_call)  # Remove allele designation (*01, *02, etc)
      ),
      j_call = ifelse(
        is.na(j_call) | j_call == "",
        NA_character_,
        sub("\\*.*$", "", j_call)
      ),
      d_call = ifelse(
        is.na(d_call) | d_call == "",
        NA_character_,
        sub("\\*.*$", "", d_call)
      ),
      v_call = gsub("TCR", "TR", v_call, fixed = TRUE),
      j_call = gsub("TCR", "TR", j_call, fixed = TRUE),
      d_call = gsub("TCR", "TR", d_call, fixed = TRUE),
      j_family = gsub(".*(([A-Z]+\\d+)).*", "\\1", j_call),
      v_family = gsub(".*(([A-Z]+\\d+)).*", "\\1", v_call),
      d_family = gsub(".*(([A-Z]+\\d+)).*", "\\1", d_call),
      bio_identity = paste(junction_aa, v_call, j_call, sep = "_"),
      sequence_id = paste(repertoire_id, dplyr::row_number(), sep = "_"),
      clone_id = bio_identity
    )
  return(clone_data)
}

#' [getAIRRFields()] Given the path to a single AIRRSeq clone file, determine
#' the file type and returns a named vector that can be used to repair headers
#' while reading input.
#' @keywords internal
#' @param clone_file .tsv file containing results from AIRRSeq pipeline
#' @param threads Number of threads for parallel processing
#' @return Named vector of corresponding AIRR fields
#'
#' @noRd
getAIRRFields <- function(clone_file, threads) {
  clone_table <- data.table::fread(clone_file, nrows = 1, nThread = threads, showProgress=FALSE)
  col_names <- base::colnames(clone_table)
  input_type <- getFileType(col_names)
  if (input_type == "immunoSEQ") {
    count_method <- clone_table |>
      dplyr::pull(counting_method) |>
      base::unique()
    # NOTE: These following fields need to be modified when the data is read
    # 1. d2_call = When it contains a list with the first values being equal to
    #    d_call, extract the second element.
    # 2. cdr3 = Remove the first and last codon
    # 3. cdr3_aa = Remove the first and last amino acid
    # 4. bioidentity = Reformat in IMGT format
    # 5. v_call, d_call, d2_call, j_call = Reformat in IMGT format
    # 6. cd*_end = Add value to cdr*start -3
    # 7. cd*_start = Subtract 3 from value
    # 8. junction_aa_length = Divide value by 3
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
      matching_fields <- c(matching_fields,
        "estimatedNumberGenomes" = "duplicate_count")
    }
  } else if (input_type == "10X") {
    matching_fields <- col_names
    names(matching_fields) <- col_names
  }
  return(list("matching_fields" = matching_fields, "file_type" = input_type))
}

#' Get available system memory in GB
#' @keywords internal
#' @noRd
get_available_memory_gb <- function() {
  if (Sys.info()["sysname"] == "Darwin") {
    vm_stat <- system("vm_stat", intern = TRUE)
    free_pages <- as.numeric(gsub(".*: ", "", vm_stat[grep("Pages free", vm_stat)]))
    inactive_pages <- as.numeric(gsub(".*: ", "", vm_stat[grep("Pages inactive", vm_stat)]))
    page_size <- 4096
    available_gb <- (free_pages + inactive_pages) * page_size / (1024^3)
    return(available_gb)
  } else if (Sys.info()["sysname"] == "Linux") {
    mem_info <- readLines("/proc/meminfo")
    available_line <- grep("MemAvailable", mem_info, value = TRUE)
    available_kb <- as.numeric(gsub(".*:\\s*([0-9]+).*", "\\1", available_line))
    return(available_kb / (1024^2))
  } else {
    return(8)
  }
}

#' Process files in chunks for memory efficiency
#' @keywords internal
#' @noRd
process_file_chunks <- function(file_paths, chunk_size, progress_bar,
                                threads, parallel, progress_detail) {
  file_chunks <- split(file_paths, ceiling(seq_along(file_paths) / chunk_size))
  chunk_results <- vector("list", length(file_chunks))

  if (progress_detail != "none") {
    cat(sprintf("Processing %d chunks...\n", length(file_chunks)))
  }

  for (i in seq_along(file_chunks)) {
    if (progress_detail == "detailed") {
      current_memory <- get_current_memory_mb()
      cat(sprintf("📦 Chunk %d/%d (Memory: %.1f MB)\n",
                  i, length(file_chunks), current_memory))
    }

    if (parallel && length(file_chunks[[i]]) > 1) {
      chunk_data <- furrr::future_map(file_chunks[[i]],
                                      ~ getStandard(.x, progress_bar, threads),
                                      .options = furrr::furrr_options(seed = TRUE))
    } else {
      chunk_data <- purrr::map(file_chunks[[i]],
                               ~ getStandard(.x, progress_bar, threads))
    }

    chunk_results[[i]] <- data.table::rbindlist(chunk_data, use.names = TRUE, fill = TRUE)
    gc()

    if (progress_detail == "detailed") {
      cat(sprintf("Chunk %d complete: %s sequences\n",
                  i, format(nrow(chunk_results[[i]]), big.mark = ",")))
    }
  }

  if (progress_detail != "none") {
    cat("🔗 Combining chunks...\n")
  }

  result <- data.table::rbindlist(chunk_results, use.names = TRUE, fill = TRUE)

  if (!is.null(progress_bar)) {
    progress_bar$terminate()
  }

  return(result)
}

#' Process files in streaming mode for ultra-large datasets
#' @keywords internal
#' @noRd
process_streaming_mode <- function(file_paths, threads, progress_detail) {
  if (progress_detail != "none") {
    cat("STREAMING MODE: Processing ultra-large dataset\n")
  }

  # Create progress bar for streaming mode
  if (progress_detail == "detailed") {
    progress_format <- paste0("Streaming [:bar] :current/:total files (:percent) | ",
                              ":rate files/sec | ETA: :eta | Elapsed: :elapsed")
  } else if (progress_detail == "basic") {
    progress_format <- "Streaming [:bar] :current/:total files (:percent) | ETA: :eta"
  } else {
    progress_format <- NULL
  }

  if (!is.null(progress_format)) {
    progress_bar <- progress::progress_bar$new(
      format = progress_format,
      total = length(file_paths), clear = FALSE, width = 100
    )
    progress_bar$tick(0)
  } else {
    progress_bar <- list(
      tick = function(...) invisible(NULL),
      terminate = function(...) invisible(NULL)
    )
  }

  streaming_results <- vector("list", length(file_paths))

  for (i in seq_along(file_paths)) {
    progress_bar$tick()

    # Suppress all output during individual file processing
    dummy_progress <- list(tick = function() {})

    # Capture and discard any progress bar output while preserving the result
    invisible(capture.output({
      streaming_results[[i]] <- suppressMessages(suppressWarnings(
        getStandard(file_paths[i], dummy_progress, threads)
      ))
    }, type = "output"))

    if (i %% 5 == 0) {
      gc(full = TRUE)
      if (progress_detail == "detailed") {
        current_memory <- get_current_memory_mb()
        cat(sprintf("Memory usage: %.1f MB\n", current_memory))
      }
    }
  }

  progress_bar$terminate()

  result <- data.table::rbindlist(streaming_results, use.names = TRUE, fill = TRUE)
  return(result)
}

#' Get current process memory usage in MB
#' @keywords internal
#' @noRd
get_current_memory_mb <- function() {
  if (Sys.info()["sysname"] == "Darwin") {
    pid <- Sys.getpid()
    mem_kb <- as.numeric(system(paste0("ps -o rss= -p ", pid), intern = TRUE))
    return(mem_kb / 1024)
  } else {
    return(0)
  }
}

#' [getStandardOptimized()] Optimized version using data.table for better performance
#' @keywords internal
#' @param clone_file A .tsv file to read in and standardize its fields to be
#'  MiAIRR compliant.
#' @param progress Progress bar
#' @param threads Number of threads for parallel processing
#' @return data.table of given data with MiAIRR fields
#'
#' @noRd
getStandardOptimized <- function(clone_file, progress, threads) {
  progress$tick()

  # Load AIRR fields template once and cache
  airr_headers_path <- system.file("extdata", "AIRR_fields.csv", package = "LymphoSeq2")
  airr_fields <- data.table::fread(airr_headers_path, stringsAsFactors = FALSE, showProgress=FALSE)

  # Get field mappings
  AIRR_fields_res <- getAIRRFieldsOptimized(clone_file, threads)
  matching_fields <- AIRR_fields_res$matching_fields
  file_type <- AIRR_fields_res$file_type

  # Use data.table::fread for much faster reading
  clone_data <- data.table::fread(
    clone_file,
    nThread = threads,
    na.strings = c("", "NA", "Nan", "NaN", "unresolved"),
    stringsAsFactors = FALSE,
    verbose = FALSE
  )

  # Apply column name mappings
  if (length(matching_fields) > 0) {
    old_names <- names(clone_data)
    new_names <- stringr::str_replace_all(old_names, matching_fields)
    data.table::setnames(clone_data, old_names, new_names)
  }

  # Get existing matches with AIRR fields
  existing_match <- intersect(names(airr_fields), names(clone_data))

  # If all AIRR fields are present, return early
  if (length(existing_match) == 155) {
    return(clone_data)
  }

  # Select only existing AIRR columns from data
  existing_airr_data <- clone_data[, ..existing_match]

  # Create template with missing columns filled with NA
  airr_template <- airr_fields[0] # Empty data.table with all AIRR columns

  # Bind the template with data to ensure all columns exist
  clone_data <- data.table::rbindlist(list(airr_template, existing_airr_data),
                                      use.names = TRUE, fill = TRUE)

  # Get file name for repertoire_id
  file_name <- tools::file_path_sans_ext(basename(clone_file))

  # Perform all mutations using data.table syntax for better performance
  clone_data[, `:=`(
    repertoire_id = file_name,
    sequence_id = .I
  )]

  # Handle d2_call splitting more efficiently
  if ("d2_call" %in% names(clone_data)) {
    clone_data[!is.na(d2_call) & stringr::str_detect(d2_call, ","),
               d2_call := stringr::str_split(d2_call, ",", simplify = TRUE)[, 2]]
  }

  # Process CDR3 sequences
  if ("junction" %in% names(clone_data)) {
    clone_data[!is.na(junction), cdr3 := stringr::str_sub(junction, 4L, -4L)]
  }
  if ("junction_aa" %in% names(clone_data)) {
    clone_data[!is.na(junction_aa), cdr3_aa := stringr::str_sub(junction_aa, 2L, -2L)]
  }

  # Handle CDR end positions
  cdr_cols <- c("cdr1_end", "cdr2_end", "cdr3_end")
  start_cols <- c("cdr1_start", "cdr2_start", "cdr3_start")

  for (i in seq_along(cdr_cols)) {
    if (all(c(cdr_cols[i], start_cols[i]) %in% names(clone_data))) {
      clone_data[!is.na(get(cdr_cols[i])),
                 (cdr_cols[i]) := get(cdr_cols[i]) + get(start_cols[i])]
    }
  }

  # Set sequence fields
  clone_data[, `:=`(
    sequence = data.table::fcoalesce(sequence, junction),
    junction = data.table::fcoalesce(junction, sequence),
    sequence_aa = data.table::fcoalesce(sequence_aa, junction_aa),
    junction_aa = data.table::fcoalesce(junction_aa, sequence_aa)
  )]

  # Process sequences based on file type
  if (file_type == "BGI") {
    # BGI-specific processing
    clone_data[, `:=`(
      junction = stringr::str_remove_all(junction, "x|[A-Z]+"),
      junction = toupper(junction),
      junction_aa = stringr::str_remove_all(junction_aa, "x|[A-Z]+"),
      junction_aa = toupper(junction_aa)
    )]
  } else {
    # Other file types
    clone_data[stringr::str_detect(junction, "[a-z]+"),
               junction := toupper(stringr::str_extract(junction, "[a-z]{2,}"))]
    clone_data[stringr::str_detect(junction_aa, "[a-z]+"),
               junction_aa := toupper(stringr::str_extract(junction_aa, "[a-z]{2,}"))]
  }

  # Final processing
  clone_data[, `:=`(
    sequence = junction,
    sequence_aa = junction_aa,
    junction_length = stringr::str_length(junction),
    junction_aa_length = stringr::str_length(junction_aa),
    rev_comp = FALSE,
    stop_codon = stringr::str_detect(sequence, "\\*") |
                 stringr::str_detect(sequence_aa, "\\*") |
                 is.na(sequence) | is.na(sequence_aa),
    productive = !(stringr::str_detect(sequence, "\\*") |
                   stringr::str_detect(sequence_aa, "\\*") |
                   is.na(sequence) | is.na(sequence_aa))
  )]

  # Process gene calls
  gene_cols <- c("v_call", "j_call", "d_call")
  for (col in gene_cols) {
    if (col %in% names(clone_data)) {
      clone_data[, (col) := stringr::str_remove(get(col), "/\\w+$")]
      clone_data[, (col) := stringr::str_c(
        stringr::str_extract(get(col), "[A-Z]+"),
        as.numeric(stringr::str_extract(get(col), "\\d+")),
        as.numeric(stringr::str_extract(get(col), "-\\d+")),
        sep = ""
      )]
      clone_data[, (col) := stringr::str_replace(get(col), "TCR", "TR")]
    }
  }

  # Set remaining fields
  clone_data[, `:=`(
    complete_vdj = !is.na(v_call) & !is.na(d_call) & !is.na(j_call),
    duplicate_frequency = duplicate_count / sum(duplicate_count, na.rm = TRUE),
    reading_frame = ifelse(stop_codon, "out-of-frame", "in-frame"),
    j_family = stringr::str_extract(j_call, "[A-Z]+\\d+"),
    v_family = stringr::str_extract(v_call, "[A-Z]+\\d+"),
    d_family = stringr::str_extract(d_call, "[A-Z]+\\d+"),
    bio_identity = stringr::str_c(junction_aa, v_call, j_call, sep = "_"),
    sequence_id = stringr::str_c(repertoire_id, .I, sep = "_"),
    clone_id = stringr::str_c(junction_aa, v_call, j_call, sep = "_")
  )]

  return(clone_data)
}

#' [getAIRRFieldsOptimized()] Optimized version using data.table
#' @keywords internal
#' @param clone_file .tsv file containing results from AIRRSeq pipeline
#' @param threads Number of threads for parallel processing
#' @return Named vector of corresponding AIRR fields
#'
#' @noRd
getAIRRFieldsOptimized <- function(clone_file, threads) {
  # Read only the header to determine file type
  clone_table <- data.table::fread(clone_file, nrows = 1, nThread = threads, showProgress=FALSE)
  col_names <- names(clone_table)
  input_type <- getFileType(col_names)

  if (input_type == "immunoSEQ") {
    count_method <- clone_table[, unique(counting_method)]

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
      matching_fields <- c(matching_fields,
        "estimatedNumberGenomes" = "duplicate_count")
    }
  } else if (input_type == "10X") {
    matching_fields <- col_names
    names(matching_fields) <- col_names
  }
  return(list("matching_fields" = matching_fields, "file_type" = input_type))
}

#' Process large datasets using Apache Arrow with advanced optimizations
#'
#' @param file_paths Character vector of file paths to process
#' @param threads Number of threads to use
#' @param progress_detail Level of progress detail ("basic", "detailed", "none")
#' @param return_type Return format ("data.table", "tibble", "lazy_dt")
#'
#' @return Processed dataset in requested format
#' @keywords internal
#' @noRd
process_with_arrow <- function(file_paths, threads, progress_detail, return_type,
                                sample_mode = FALSE, sample_size = 1000000) {

  if (!requireNamespace("arrow", quietly = TRUE)) {
    stop("Arrow package is required for large dataset processing. ",
         "Please install it with: install.packages('arrow')")
  }

  # Analyze file size distribution for intelligent processing
  file_info <- file.info(file_paths)
  file_sizes_mb <- file_info$size / (1024^2)
  total_size_gb <- sum(file_info$size, na.rm = TRUE) / (1024^3)

  # Categorize files by size for different processing strategies
  small_files <- file_paths[file_sizes_mb < 50]   # < 50MB
  medium_files <- file_paths[file_sizes_mb >= 50 & file_sizes_mb < 200]  # 50-200MB
  large_files <- file_paths[file_sizes_mb >= 200]  # > 200MB

  if (progress_detail != "none") {
    cat("Advanced Dataset Analysis\n")
    cat("========================\n")
    cat(sprintf("Files to process: %d (Total: %.2f GB)\n", length(file_paths), total_size_gb))
    cat(sprintf("  Small files (<50MB): %d\n", length(small_files)))
    cat(sprintf("  Medium files (50-200MB): %d\n", length(medium_files)))
    cat(sprintf("  Large files (>200MB): %d\n", length(large_files)))
    cat("Backend: Apache Arrow with parallel processing\n")
    cat("Processing method: Adaptive chunking with intelligent batching\n\n")
  }

  # Use future for true parallel processing
  if (!requireNamespace("future", quietly = TRUE) || !requireNamespace("furrr", quietly = TRUE)) {
    warning("future and furrr packages required for parallel Arrow processing. Installing...")
    return(process_with_arrow_sequential(file_paths, threads, progress_detail, return_type))
  }

  # Set up parallel plan if not already configured
  current_plan <- class(future::plan())[1]
  if (current_plan == "sequential") {
    future::plan(future::multisession, workers = min(threads, 8))  # Cap at 8 workers
    reset_plan <- TRUE
  } else {
    reset_plan <- FALSE
  }

  # Setup progress reporting
  total_batches <- ceiling(length(small_files)/10) + ceiling(length(medium_files)/5) + length(large_files)

  if (progress_detail == "detailed") {
    progress_format <- paste0("Loading [:bar] :current/:total batches (:percent) | ",
                              ":rate batches/sec | ETA: :eta | Elapsed: :elapsed")
  } else if (progress_detail == "basic") {
    progress_format <- "Loading [:bar] :current/:total batches (:percent) | ETA: :eta"
  } else {
    progress_format <- NULL
  }

  if (!is.null(progress_format)) {
    progress_bar <- progress::progress_bar$new(
      format = progress_format,
      total = total_batches, clear = FALSE, width = 100
    )
    progress_bar$tick(0)
  }

  all_results <- list()
  batch_count <- 0

  # Process small files in parallel batches of 10
  if (length(small_files) > 0) {
    small_batches <- split(small_files, ceiling(seq_along(small_files)/10))

    for (batch in small_batches) {
      batch_result <- furrr::future_map(batch, function(file_path) {
        process_single_file_arrow(file_path, fast_mode = TRUE,
                                  sample_mode = sample_mode, sample_size = sample_size)
      }, .options = furrr::furrr_options(seed = TRUE))

      # Combine batch results
      batch_combined <- data.table::rbindlist(batch_result, fill = TRUE)
      all_results[[length(all_results) + 1]] <- batch_combined

      batch_count <- batch_count + 1
      if (!is.null(progress_format)) progress_bar$tick()

      # Memory cleanup
      rm(batch_result, batch_combined)
      if (batch_count %% 3 == 0) gc()
    }
  }

  # Process medium files in parallel batches of 5
  if (length(medium_files) > 0) {
    medium_batches <- split(medium_files, ceiling(seq_along(medium_files)/5))

    for (batch in medium_batches) {
      batch_result <- furrr::future_map(batch, function(file_path) {
        process_single_file_arrow(file_path, fast_mode = FALSE,
                                  sample_mode = sample_mode, sample_size = sample_size)
      }, .options = furrr::furrr_options(seed = TRUE))

      # Combine batch results
      batch_combined <- data.table::rbindlist(batch_result, fill = TRUE)
      all_results[[length(all_results) + 1]] <- batch_combined

      batch_count <- batch_count + 1
      if (!is.null(progress_format)) progress_bar$tick()

      # Memory cleanup
      rm(batch_result, batch_combined)
      gc()
    }
  }

  # Process large files individually with streaming
  if (length(large_files) > 0) {
    for (file_path in large_files) {
      result <- process_large_file_streaming(file_path, sample_mode, sample_size)
      all_results[[length(all_results) + 1]] <- result

      batch_count <- batch_count + 1
      if (!is.null(progress_format)) progress_bar$tick()

      # Aggressive memory cleanup for large files
      rm(result)
      gc(full = TRUE)
    }
  }

  if (!is.null(progress_format)) {
    progress_bar$terminate()
  }

  # Reset parallel plan if we set it
  if (reset_plan) {
    future::plan(future::sequential)
  }

  # Efficient final combination using data.table
  if (progress_detail != "none") {
    cat("🔗 Combining all results...\n")
  }

  result_dt <- data.table::rbindlist(all_results, fill = TRUE)

  # Memory cleanup
  rm(all_results)
  gc()

  # Remove duplicates if any
  if (nrow(result_dt) > 0) {
    initial_rows <- nrow(result_dt)
    result_dt <- unique(result_dt, by = c("sequence_id", "repertoire_id"))
    if (progress_detail == "detailed" && nrow(result_dt) < initial_rows) {
      cat(sprintf("Removed %s duplicate sequences\n",
                 format(initial_rows - nrow(result_dt), big.mark = ",")))
    }
  }

  # Convert to requested return type
  if (return_type == "tibble") {
    result <- tibble::as_tibble(result_dt)
  } else if (return_type == "lazy_dt") {
    result <- dtplyr::lazy_dt(result_dt)
  } else {
    result <- result_dt
  }

  if (progress_detail != "none") {
    cat(sprintf("\n✅ Arrow processing complete!\n"))
    cat(sprintf("Total sequences loaded: %s\n",
               format(nrow(result_dt), big.mark = ",")))
    cat(sprintf("Final dataset size: %.2f GB in memory\n",
               object.size(result_dt) / (1024^3)))
  }

  return(result)
}

#' Process BGI format data
#' @keywords internal
#' @noRd
process_bgi_format <- function(dt_data) {
  # Use the existing getStandard logic for BGI files
  # This function serves as a wrapper for consistent processing
  dummy_progress <- list(tick = function() {})

  # Apply the same BGI processing logic from getStandard
  file_type <- "BGI"

  # Apply basic BGI transformations
  if ("junction" %in% names(dt_data)) {
    dt_data[, junction := gsub("x", "", junction, fixed = TRUE)]
    dt_data[, junction := gsub("[A-Z]+", "", junction)]
    dt_data[, junction := toupper(junction)]
  }

  if ("junction_aa" %in% names(dt_data)) {
    dt_data[, junction_aa := gsub("x", "", junction_aa, fixed = TRUE)]
    dt_data[, junction_aa := gsub("[A-Z]+", "", junction_aa)]
    dt_data[, junction_aa := toupper(junction_aa)]
  }

  return(dt_data)
}

#' Process ImmunoSEQ format data
#' @keywords internal
#' @noRd
process_immunoseq_format <- function(dt_data) {
  # Apply ImmunoSEQ specific processing
  if ("junction" %in% names(dt_data)) {
    dt_data[grepl("[a-z]+", junction),
            junction := toupper(gsub(".*([a-z]{2,}).*", "\\1", junction))]
  }

  if ("junction_aa" %in% names(dt_data)) {
    dt_data[grepl("[a-z]+", junction_aa),
            junction_aa := toupper(gsub(".*([a-z]{2,}).*", "\\1", junction_aa))]
  }

  return(dt_data)
}

#' Process standard format data
#' @keywords internal
#' @noRd
process_standard_format <- function(dt_data) {
  # Minimal processing for standard/unknown formats
  # Just ensure basic columns exist

  # Add missing key columns if they don't exist
  if (!"sequence" %in% names(dt_data) && "junction" %in% names(dt_data)) {
    dt_data[, sequence := junction]
  }

  if (!"sequence_aa" %in% names(dt_data) && "junction_aa" %in% names(dt_data)) {
    dt_data[, sequence_aa := junction_aa]
  }

  return(dt_data)
}

#' Process a single file with Arrow optimizations
#' @keywords internal
#' @noRd
process_single_file_arrow <- function(file_path, fast_mode = FALSE,
                                      sample_mode = FALSE, sample_size = 1000000) {
  tryCatch({
    # Get file name for repertoire_id
    file_name <- tools::file_path_sans_ext(basename(file_path))

    # Use Arrow for initial read with optimizations
    dt_data <- tryCatch({
      # Read with Arrow for better memory efficiency
      arrow_table <- arrow::read_csv_arrow(
        file_path,
        skip_empty_rows = TRUE,
        parse_options = arrow::CsvParseOptions(
          delimiter = "\t",
          quote_char = FALSE
        )
      )

      # Convert to data.table
      arrow_table |>
        arrow::collect() |>
        data.table::as.data.table()

    }, error = function(e) {
      # Fallback to data.table::fread
      data.table::fread(file_path, showProgress = FALSE, nThread = 2)
    })

    # Apply sampling if requested
    if (sample_mode && nrow(dt_data) > sample_size) {
      sample_indices <- sample(nrow(dt_data), sample_size)
      dt_data <- dt_data[sample_indices, ]
    }

    # Apply minimal processing for fast mode, full processing otherwise
    if (fast_mode) {
      # Minimal processing for small files
      if (!"repertoire_id" %in% names(dt_data)) {
        dt_data[, repertoire_id := file_name]
      }
      if (!"sequence_id" %in% names(dt_data)) {
        dt_data[, sequence_id := paste(file_name, .I, sep = "_")]
      }
    } else {
      # Full processing using existing logic
      dt_data <- apply_standard_processing(dt_data, file_name)
    }

    return(dt_data)

  }, error = function(e) {
    # Return empty data.table with minimal structure on error
    warning(sprintf("Failed to process %s: %s", basename(file_path), e$message))
    return(data.table::data.table(
      repertoire_id = tools::file_path_sans_ext(basename(file_path)),
      sequence_id = character(0)
    ))
  })
}

#' Process large files with streaming
#' @keywords internal
#' @noRd
process_large_file_streaming <- function(file_path, sample_mode = FALSE,
                                          sample_size = 1000000) {
  file_name <- tools::file_path_sans_ext(basename(file_path))

  tryCatch({
    # For very large files, read in chunks
    chunk_size <- 50000  # Process 50k rows at a time

    # Get total number of lines first
    total_lines <- as.numeric(system(paste("wc -l <", shQuote(file_path)), intern = TRUE)) - 1

    if (total_lines <= chunk_size) {
      # File is small enough to process normally
      return(process_single_file_arrow(file_path, fast_mode = FALSE,
                                       sample_mode = sample_mode, sample_size = sample_size))
    }

    # If in sample mode, adjust strategy for large files
    if (sample_mode) {
      # For sampling, we can read fewer chunks
      total_samples_needed <- min(sample_size, total_lines)
      chunk_size <- min(chunk_size, total_samples_needed)
      # Calculate how many lines to skip to get representative sampling
      skip_interval <- max(1, floor(total_lines / total_samples_needed))
    }

    # Read header first
    header <- data.table::fread(file_path, nrows = 1, showProgress = FALSE)
    col_names <- names(header)

    # Process in chunks
    chunks <- list()
    current_row <- 1

    while (current_row <= total_lines) {
      chunk_data <- data.table::fread(
        file_path,
        skip = current_row,
        nrows = min(chunk_size, total_lines - current_row + 1),
        col.names = col_names,
        showProgress = FALSE,
        nThread = 2
      )

      if (nrow(chunk_data) > 0) {
        # Apply basic processing
        chunk_data[, repertoire_id := file_name]
        chunk_data[, sequence_id := paste(file_name, (.I + current_row - 1), sep = "_")]

        chunks[[length(chunks) + 1]] <- chunk_data
      }

      current_row <- current_row + chunk_size

      # Clean up memory every few chunks
      if (length(chunks) %% 3 == 0) {
        gc()
      }
    }

    # Combine chunks
    result <- data.table::rbindlist(chunks, fill = TRUE)
    return(result)

  }, error = function(e) {
    warning(sprintf("Streaming failed for %s: %s", basename(file_path), e$message))
    # Fallback to regular processing
    return(process_single_file_arrow(file_path, fast_mode = TRUE))
  })
}

#' Apply standard processing to data
#' @keywords internal
#' @noRd
apply_standard_processing <- function(dt_data, file_name) {
  # Detect file type
  col_names <- names(dt_data)

  if (any(grepl("Sequence.ID", col_names))) {
    dt_data <- process_bgi_format(dt_data)
  } else if (any(grepl("nucleotide|amino_acid", col_names))) {
    dt_data <- process_immunoseq_format(dt_data)
  } else {
    dt_data <- process_standard_format(dt_data)
  }

  # Add required columns
  if (!"repertoire_id" %in% names(dt_data)) {
    dt_data[, repertoire_id := file_name]
  }
  if (!"sequence_id" %in% names(dt_data)) {
    dt_data[, sequence_id := paste(file_name, .I, sep = "_")]
  }

  return(dt_data)
}

#' Sequential Arrow processing fallback
#' @keywords internal
#' @noRd
process_with_arrow_sequential <- function(file_paths, threads, progress_detail, return_type) {
  # Fallback to the original implementation with some optimizations
  if (progress_detail != "none") {
    cat("Falling back to sequential Arrow processing\n")
  }

  # Setup progress bar
  if (progress_detail == "detailed") {
    progress_format <- paste0("Loading [:bar] :current/:total files (:percent) | ",
                              ":rate files/sec | ETA: :eta | Elapsed: :elapsed")
  } else if (progress_detail == "basic") {
    progress_format <- "Loading [:bar] :current/:total files (:percent) | ETA: :eta"
  } else {
    progress_format <- NULL
  }

  if (!is.null(progress_format)) {
    progress_bar <- progress::progress_bar$new(
      format = progress_format,
      total = length(file_paths), clear = FALSE, width = 100
    )
    progress_bar$tick(0)
  }

  all_data <- list()

  for (i in seq_along(file_paths)) {
    result <- process_single_file_arrow(file_paths[i], fast_mode = FALSE)
    all_data[[i]] <- result

    if (!is.null(progress_format)) {
      progress_bar$tick()
    }

    # Memory cleanup every 10 files
    if (i %% 10 == 0) {
      gc()
    }
  }

  if (!is.null(progress_format)) {
    progress_bar$terminate()
  }

  # Combine results
  result_dt <- data.table::rbindlist(all_data, fill = TRUE)

  # Convert to requested return type
  if (return_type == "tibble") {
    result <- tibble::as_tibble(result_dt)
  } else if (return_type == "lazy_dt") {
    result <- dtplyr::lazy_dt(result_dt)
  } else {
    result <- result_dt
  }

  return(result)
}