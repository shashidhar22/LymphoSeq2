context("Comprehensive readImmunoSeq Tests - All Formats and AIRR Compliance")
library(LymphoSeq2)
library(testthat)
library(data.table)

# Helper function to get AIRR fields
get_expected_airr_fields <- function() {
  airr_headers_path <- system.file("extdata", "AIRR_fields.csv", package = "LymphoSeq2")
  airr_table <- data.table::fread(airr_headers_path, showProgress = FALSE)
  return(names(airr_table))
}

# ============================================================================
# 1. File Format Detection Tests
# ============================================================================

test_that("getFileType correctly identifies BGI format", {
  test_file <- "test_data/sample_BGI.tsv"
  skip_if_not(file.exists(test_file), "BGI test file not found")

  # Read header
  header <- data.table::fread(test_file, nrows = 0, showProgress = FALSE)
  col_names <- names(header)

  # Test internal function
  file_type <- LymphoSeq2:::getFileType(col_names)
  expect_equal(file_type, "BGI")
})

test_that("getFileType correctly identifies ImmunoSEQ Legacy format", {
  test_file <- "test_data/015V12001549_CFAR.tsv"

  header <- data.table::fread(test_file, nrows = 0, showProgress = FALSE)
  col_names <- names(header)

  file_type <- LymphoSeq2:::getFileType(col_names)
  expect_equal(file_type, "immunoSEQLegacy")
})

test_that("getFileType correctly identifies ImmunoSEQ V3 format", {
  test_file <- "test_data/sample_TRB_V3.tsv"
  skip_if_not(file.exists(test_file), "TRB V3 test file not found")

  header <- data.table::fread(test_file, nrows = 0, showProgress = FALSE)
  col_names <- names(header)

  file_type <- LymphoSeq2:::getFileType(col_names)
  expect_equal(file_type, "immunoSEQ")
})

# ============================================================================
# 2. BGI File Processing and AIRR Mapping
# ============================================================================

test_that("readImmunoSeq processes BGI files correctly", {
  test_file <- "test_data/sample_BGI.tsv"
  skip_if_not(file.exists(test_file), "BGI test file not found")

  result <- readImmunoSeq(
    test_file,
    threads = 1,
    parallel = FALSE,
    progress_detail = "none"
  )

  # Basic structure checks
  expect_s3_class(result, "data.table")
  expect_true(nrow(result) > 0)
  expect_true("repertoire_id" %in% names(result))
  expect_true("junction" %in% names(result))
  expect_true("junction_aa" %in% names(result))

  # Check repertoire_id is set
  expect_equal(unique(result$repertoire_id), "sample_BGI")
})

test_that("BGI junction sequences are properly cleaned", {
  test_file <- "test_data/sample_BGI.tsv"
  skip_if_not(file.exists(test_file), "BGI test file not found")

  result <- readImmunoSeq(test_file, threads = 1, progress_detail = "none")

  # BGI-specific processing checks - should remove lowercase and 'x'
  junctions <- result$junction[!is.na(result$junction)]

  # Should be uppercase nucleotides only (N allowed for ambiguous bases per IUPAC)
  expect_false(any(grepl("x", junctions, fixed = TRUE)))
  expect_false(any(grepl("[a-z]", junctions)))
  expect_true(all(grepl("^[ACGTN]*$", junctions)))

  # Check junction_aa
  junctions_aa <- result$junction_aa[!is.na(result$junction_aa)]
  expect_false(any(grepl("x", junctions_aa, fixed = TRUE)))
  expect_true(all(grepl("^[A-Z\\*]*$", junctions_aa)))
})

test_that("BGI files map to essential AIRR fields", {
  test_file <- "test_data/sample_BGI.tsv"
  skip_if_not(file.exists(test_file), "BGI test file not found")

  result <- readImmunoSeq(test_file, threads = 1, progress_detail = "none")

  # Essential AIRR fields
  essential_fields <- c(
    "sequence_id", "sequence", "sequence_aa",
    "v_call", "d_call", "j_call",
    "junction", "junction_aa", "junction_length", "junction_aa_length",
    "duplicate_count", "duplicate_frequency",
    "repertoire_id", "productive", "reading_frame",
    "v_family", "j_family", "d_family",
    "bio_identity", "clone_id"
  )

  result_fields <- names(result)
  missing_fields <- setdiff(essential_fields, result_fields)

  expect_equal(length(missing_fields), 0,
               info = paste("Missing fields:", paste(missing_fields, collapse = ", ")))
})

test_that("BGI gene calls are converted to IMGT format", {
  test_file <- "test_data/sample_BGI.tsv"
  skip_if_not(file.exists(test_file), "BGI test file not found")

  result <- readImmunoSeq(test_file, threads = 1, progress_detail = "none")

  # Check V gene format
  v_calls <- result$v_call[!is.na(result$v_call)]
  if (length(v_calls) > 0) {
    expect_true(all(grepl("^TR[ABGD]V\\d+", v_calls)),
                info = paste("Invalid V calls:", paste(head(v_calls[!grepl("^TR[ABGD]V\\d+", v_calls)]), collapse = ", ")))
  }

  # Check J gene format
  j_calls <- result$j_call[!is.na(result$j_call)]
  if (length(j_calls) > 0) {
    expect_true(all(grepl("^TR[ABGD]J\\d+", j_calls)),
                info = paste("Invalid J calls:", paste(head(j_calls[!grepl("^TR[ABGD]J\\d+", j_calls)]), collapse = ", ")))
  }

  # Check D gene format
  d_calls <- result$d_call[!is.na(result$d_call)]
  if (length(d_calls) > 0) {
    expect_true(all(grepl("^TR[ABGD]D\\d+", d_calls)),
                info = paste("Invalid D calls:", paste(head(d_calls[!grepl("^TR[ABGD]D\\d+", d_calls)]), collapse = ", ")))
  }
})

# ============================================================================
# 3. ImmunoSEQ Legacy File Processing and AIRR Mapping
# ============================================================================

test_that("readImmunoSeq processes ImmunoSEQ Legacy files correctly", {
  test_file <- "test_data/015V12001549_CFAR.tsv"

  result <- readImmunoSeq(test_file, threads = 1, progress_detail = "none")

  expect_s3_class(result, "data.table")
  expect_true(nrow(result) > 0)
  expect_equal(unique(result$repertoire_id), "015V12001549_CFAR")
})

test_that("ImmunoSEQ Legacy duplicate count mapping is correct", {
  test_file <- "test_data/015V12001549_CFAR.tsv"

  result <- readImmunoSeq(test_file, threads = 1, progress_detail = "none")

  # Check duplicate_count exists and is numeric
  expect_true("duplicate_count" %in% names(result))
  expect_true(is.numeric(result$duplicate_count))
  expect_true(all(result$duplicate_count > 0, na.rm = TRUE))
})

test_that("ImmunoSEQ Legacy junction extraction works", {
  test_file <- "test_data/015V12001549_CFAR.tsv"

  result <- readImmunoSeq(test_file, threads = 1, progress_detail = "none")

  # Legacy format uses nucleotide/aminoAcid columns with CDR3 in lowercase
  # Function should extract and uppercase the CDR3
  junctions <- result$junction[!is.na(result$junction)]

  # Should be uppercase
  expect_true(all(junctions == toupper(junctions)))
  expect_true(all(grepl("^[ACGT]+$", junctions)))
})

# ============================================================================
# 4. ImmunoSEQ V3 File Processing and AIRR Mapping
# ============================================================================

test_that("readImmunoSeq processes ImmunoSEQ V3 TRB files correctly", {
  test_file <- "test_data/sample_TRB_V3.tsv"
  skip_if_not(file.exists(test_file), "TRB V3 test file not found")

  result <- readImmunoSeq(test_file, threads = 1, progress_detail = "none")

  expect_s3_class(result, "data.table")
  expect_true(nrow(result) > 0)

  # V3 format should have more detailed annotations
  expect_true("v_call" %in% names(result))
  expect_true("d_call" %in% names(result))
  expect_true("j_call" %in% names(result))
  expect_true("bio_identity" %in% names(result))
})

test_that("ImmunoSEQ V3 CDR3 extraction is correct", {
  test_file <- "test_data/sample_TRB_V3.tsv"
  skip_if_not(file.exists(test_file), "TRB V3 test file not found")

  result <- readImmunoSeq(test_file, threads = 1, progress_detail = "none")

  # CDR3 should be junction minus first and last codon
  valid_rows <- !is.na(result$junction) & !is.na(result$cdr3) &
                nchar(result$junction) >= 7

  if (sum(valid_rows) > 10) {
    subset_data <- result[valid_rows, ]

    # Check that cdr3 is derived from junction
    for (i in 1:min(10, nrow(subset_data))) {
      junction <- subset_data$junction[i]
      cdr3 <- subset_data$cdr3[i]

      # CDR3 should be substring of junction
      if (!is.na(cdr3) && nchar(cdr3) > 0) {
        expect_true(grepl(cdr3, junction, fixed = TRUE))
      }
    }
  }
})

# ============================================================================
# 5. TRA File Processing (No D gene)
# ============================================================================

test_that("readImmunoSeq processes TRA files correctly (no D gene)", {
  test_file <- "test_data/sample_TRA.tsv"
  skip_if_not(file.exists(test_file), "TRA test file not found")

  result <- readImmunoSeq(test_file, threads = 1, progress_detail = "none")

  expect_s3_class(result, "data.table")
  expect_true(nrow(result) > 0)

  # TRA should have V and J but D gene may be NA or absent
  expect_true("v_call" %in% names(result))
  expect_true("j_call" %in% names(result))
  expect_true("d_call" %in% names(result))  # Column should exist

  # complete_vdj should be FALSE for TRA (no D gene)
  # Actually, for TRA, complete_vdj logic may differ
  expect_true("complete_vdj" %in% names(result))
})

# ============================================================================
# 6. IGH File Processing (Heavy chain)
# ============================================================================

test_that("readImmunoSeq processes IGH files correctly", {
  test_file <- "test_data/sample_IGH.tsv"
  skip_if_not(file.exists(test_file), "IGH test file not found")

  result <- readImmunoSeq(test_file, threads = 1, progress_detail = "none")

  expect_s3_class(result, "data.table")
  expect_true(nrow(result) > 0)

  # IGH should have V, D, and J genes
  expect_true("v_call" %in% names(result))
  expect_true("d_call" %in% names(result))
  expect_true("j_call" %in% names(result))
})

# ============================================================================
# 7. AIRR Standard Compliance Tests
# ============================================================================

test_that("All file formats produce AIRR-compliant output structure", {
  test_files <- c(
    "test_data/sample_BGI.tsv",
    "test_data/015V12001549_CFAR.tsv",
    "test_data/sample_TRB_V3.tsv"
  )

  existing_files <- test_files[file.exists(test_files)]
  skip_if(length(existing_files) == 0, "No test files available")

  expected_fields <- get_expected_airr_fields()

  for (test_file in existing_files) {
    result <- readImmunoSeq(test_file, threads = 1, progress_detail = "none")
    result_fields <- names(result)

    # Check that most AIRR fields are present (allow for some optional fields)
    common_fields <- intersect(expected_fields, result_fields)
    coverage <- length(common_fields) / length(expected_fields)

    expect_true(coverage > 0.8,
                info = paste(basename(test_file), "AIRR coverage:",
                           round(coverage * 100, 1), "%"))
  }
})

test_that("duplicate_frequency sums to 1.0 per repertoire", {
  test_files <- c(
    "test_data/sample_BGI.tsv",
    "test_data/015V12001549_CFAR.tsv"
  )

  for (test_file in test_files) {
    if (!file.exists(test_file)) next

    result <- readImmunoSeq(test_file, threads = 1, progress_detail = "none")

    freq_sum <- result[, sum(duplicate_frequency, na.rm = TRUE), by = repertoire_id]

    expect_true(all(abs(freq_sum$V1 - 1.0) < 0.001),
                info = paste(basename(test_file), "frequency sum:",
                           paste(freq_sum$V1, collapse = ", ")))
  }
})

test_that("productive flag is consistent with stop_codon", {
  test_files <- c(
    "test_data/sample_BGI.tsv",
    "test_data/015V12001549_CFAR.tsv"
  )

  for (test_file in test_files) {
    if (!file.exists(test_file)) next

    result <- readImmunoSeq(test_file, threads = 1, progress_detail = "none")

    # productive should be opposite of stop_codon
    inconsistent <- which(result$productive == result$stop_codon &
                         !is.na(result$productive) &
                         !is.na(result$stop_codon))

    expect_equal(length(inconsistent), 0,
                info = paste(basename(test_file), "has", length(inconsistent),
                           "inconsistent productive/stop_codon entries"))
  }
})

test_that("reading_frame is consistent with productive", {
  test_files <- c(
    "test_data/sample_BGI.tsv",
    "test_data/015V12001549_CFAR.tsv"
  )

  for (test_file in test_files) {
    if (!file.exists(test_file)) next

    result <- readImmunoSeq(test_file, threads = 1, progress_detail = "none")

    # in-frame should correspond to productive = TRUE
    inconsistent <- which(
      (result$productive & result$reading_frame != "in-frame") |
      (!result$productive & result$reading_frame != "out-of-frame")
    )

    expect_true(length(inconsistent) == 0,
                info = paste(basename(test_file), "has", length(inconsistent),
                           "inconsistent productive/reading_frame entries"))
  }
})

test_that("bio_identity field is correctly formatted", {
  test_files <- c(
    "test_data/sample_BGI.tsv",
    "test_data/015V12001549_CFAR.tsv"
  )

  for (test_file in test_files) {
    if (!file.exists(test_file)) next

    result <- readImmunoSeq(test_file, threads = 1, progress_detail = "none")

    # bio_identity should be junction_aa_v_call_j_call
    valid_bio <- !is.na(result$bio_identity) &
                 !is.na(result$junction_aa) &
                 !is.na(result$v_call) &
                 !is.na(result$j_call)

    if (sum(valid_bio) > 0) {
      subset_data <- result[valid_bio, ]

      # Check format
      expect_true(all(grepl("_", subset_data$bio_identity)),
                  info = paste(basename(test_file), "bio_identity missing underscores"))
    }
  }
})

# ============================================================================
# 8. Multiple File Processing Tests
# ============================================================================

test_that("readImmunoSeq processes multiple files correctly", {
  test_files <- c(
    "test_data/015V12001549_CFAR.tsv",
    "test_data/015V12001685_CFAR_R.tsv"
  )

  result <- readImmunoSeq(
    test_files,
    threads = 1,
    parallel = FALSE,
    progress_detail = "none"
  )

  expect_s3_class(result, "data.table")

  # Check both repertoires are present
  repertoires <- unique(result$repertoire_id)
  expect_equal(length(repertoires), 2)
  expect_true("015V12001549_CFAR" %in% repertoires)
  expect_true("015V12001685_CFAR_R" %in% repertoires)
})

test_that("readImmunoSeq processes directory correctly", {
  result <- readImmunoSeq(
    "test_data/",
    threads = 1,
    recursive = FALSE,
    progress_detail = "none"
  )

  expect_s3_class(result, "data.table")
  expect_true(nrow(result) > 0)

  # Should have multiple repertoires
  repertoires <- unique(result$repertoire_id)
  expect_true(length(repertoires) >= 4)  # At least the 4 original files
})

# ============================================================================
# 9. Return Type Tests
# ============================================================================

test_that("readImmunoSeq returns data.table when requested", {
  test_file <- "test_data/015V12001549_CFAR.tsv"

  result <- readImmunoSeq(
    test_file,
    return_type = "data.table",
    threads = 1,
    progress_detail = "none"
  )

  expect_s3_class(result, "data.table")
})

test_that("readImmunoSeq returns tibble when requested", {
  test_file <- "test_data/015V12001549_CFAR.tsv"

  result <- readImmunoSeq(
    test_file,
    return_type = "tibble",
    threads = 1,
    progress_detail = "none"
  )

  expect_s3_class(result, "tbl_df")
  expect_s3_class(result, "tbl")
})

test_that("readImmunoSeq returns lazy_dt when requested", {
  test_file <- "test_data/015V12001549_CFAR.tsv"

  result <- readImmunoSeq(
    test_file,
    return_type = "lazy_dt",
    threads = 1,
    progress_detail = "none"
  )

  expect_s3_class(result, "dtplyr_step")
})

# ============================================================================
# 10. Data Quality Tests
# ============================================================================

test_that("sequence_id is unique within each repertoire", {
  test_file <- "test_data/015V12001549_CFAR.tsv"

  result <- readImmunoSeq(test_file, threads = 1, progress_detail = "none")

  # Check for duplicates within repertoire
  dup_check <- result[, .(n = .N), by = .(repertoire_id, sequence_id)]
  duplicates <- dup_check[n > 1]

  expect_equal(nrow(duplicates), 0,
               info = paste("Found", nrow(duplicates), "duplicate sequence_ids"))
})

test_that("junction_length matches actual junction length", {
  test_files <- c(
    "test_data/sample_BGI.tsv",
    "test_data/015V12001549_CFAR.tsv"
  )

  for (test_file in test_files) {
    if (!file.exists(test_file)) next

    result <- readImmunoSeq(test_file, threads = 1, progress_detail = "none")

    valid_rows <- !is.na(result$junction) & !is.na(result$junction_length)

    if (sum(valid_rows) > 0) {
      subset_data <- result[valid_rows, ]
      actual_lengths <- nchar(subset_data$junction)

      mismatches <- which(subset_data$junction_length != actual_lengths)

      expect_equal(length(mismatches), 0,
                   info = paste(basename(test_file), "has", length(mismatches),
                              "junction_length mismatches"))
    }
  }
})

test_that("junction_aa_length matches actual junction_aa length", {
  test_files <- c(
    "test_data/sample_BGI.tsv",
    "test_data/015V12001549_CFAR.tsv"
  )

  for (test_file in test_files) {
    if (!file.exists(test_file)) next

    result <- readImmunoSeq(test_file, threads = 1, progress_detail = "none")

    valid_rows <- !is.na(result$junction_aa) & !is.na(result$junction_aa_length)

    if (sum(valid_rows) > 0) {
      subset_data <- result[valid_rows, ]
      actual_lengths <- nchar(subset_data$junction_aa)

      mismatches <- which(subset_data$junction_aa_length != actual_lengths)

      expect_equal(length(mismatches), 0,
                   info = paste(basename(test_file), "has", length(mismatches),
                              "junction_aa_length mismatches"))
    }
  }
})

test_that("Gene family extraction is consistent", {
  test_file <- "test_data/sample_BGI.tsv"
  skip_if_not(file.exists(test_file), "BGI test file not found")

  result <- readImmunoSeq(test_file, threads = 1, progress_detail = "none")

  # Check v_family is derived from v_call
  valid_v <- !is.na(result$v_call) & !is.na(result$v_family)

  if (sum(valid_v) > 0) {
    subset_data <- result[valid_v, ]

    # v_family should be prefix of v_call (e.g., TRBV7 from TRBV7-9)
    for (i in 1:min(10, nrow(subset_data))) {
      v_call <- subset_data$v_call[i]
      v_family <- subset_data$v_family[i]

      if (!is.na(v_family) && nchar(v_family) > 0) {
        expect_true(grepl(v_family, v_call, fixed = TRUE),
                   info = paste("v_family", v_family, "not found in v_call", v_call))
      }
    }
  }
})

# ============================================================================
# 11. Performance and Threading Tests
# ============================================================================

test_that("readImmunoSeq works with threads = 1", {
  test_file <- "test_data/015V12001549_CFAR.tsv"

  result <- readImmunoSeq(
    test_file,
    threads = 1,
    progress_detail = "none"
  )

  expect_s3_class(result, "data.table")
  expect_true(nrow(result) > 0)
})

test_that("readImmunoSeq works with threads > 1", {
  test_file <- "test_data/015V12001549_CFAR.tsv"

  result <- readImmunoSeq(
    test_file,
    threads = 2,
    progress_detail = "none"
  )

  expect_s3_class(result, "data.table")
  expect_true(nrow(result) > 0)
})

test_that("Parallel and sequential processing produce same row counts", {
  test_files <- c(
    "test_data/015V12001549_CFAR.tsv",
    "test_data/015V12001685_CFAR_R.tsv"
  )

  result_seq <- readImmunoSeq(
    test_files,
    threads = 1,
    parallel = FALSE,
    progress_detail = "none"
  )

  result_par <- readImmunoSeq(
    test_files,
    threads = 2,
    parallel = TRUE,
    progress_detail = "none"
  )

  expect_equal(nrow(result_seq), nrow(result_par))
  expect_equal(length(unique(result_seq$repertoire_id)),
               length(unique(result_par$repertoire_id)))
})

# ============================================================================
# 12. Error Handling Tests
# ============================================================================

test_that("readImmunoSeq handles non-existent files gracefully", {
  expect_error(
    readImmunoSeq("/nonexistent/file.tsv", threads = 1),
    regexp = ".*"
  )
})

test_that("readImmunoSeq handles empty directory gracefully", {
  temp_dir <- tempfile()
  dir.create(temp_dir)

  expect_error(
    readImmunoSeq(temp_dir, threads = 1),
    regexp = ".*"
  )

  unlink(temp_dir, recursive = TRUE)
})

# ============================================================================
# 13. Edge Cases
# ============================================================================

test_that("readImmunoSeq handles files with NA values", {
  test_file <- "test_data/015V12001549_CFAR.tsv"

  result <- readImmunoSeq(test_file, threads = 1, progress_detail = "none")

  # Should handle NA values gracefully
  expect_true("productive" %in% names(result))
  expect_true(is.logical(result$productive) || all(result$productive %in% c(TRUE, FALSE, NA)))
})

test_that("readImmunoSeq handles mixed file formats in directory", {
  # The test_data directory contains both legacy and potentially other formats
  result <- readImmunoSeq(
    "test_data/",
    threads = 1,
    recursive = FALSE,
    progress_detail = "none"
  )

  expect_s3_class(result, "data.table")

  # Should successfully process multiple formats
  repertoires <- unique(result$repertoire_id)
  expect_true(length(repertoires) > 1)
})

# ============================================================================
# 14. Progress Reporting Tests
# ============================================================================

test_that("readImmunoSeq works with progress_detail = none", {
  test_file <- "test_data/015V12001549_CFAR.tsv"

  result <- readImmunoSeq(
    test_file,
    threads = 1,
    progress_detail = "none"
  )

  expect_s3_class(result, "data.table")
})

test_that("readImmunoSeq works with progress_detail = basic", {
  test_file <- "test_data/015V12001549_CFAR.tsv"

  expect_output(
    result <- readImmunoSeq(
      test_file,
      threads = 1,
      progress_detail = "basic"
    ),
    regexp = ".*"
  )

  expect_s3_class(result, "data.table")
})

test_that("readImmunoSeq works with progress_detail = detailed", {
  test_file <- "test_data/015V12001549_CFAR.tsv"

  expect_output(
    result <- readImmunoSeq(
      test_file,
      threads = 1,
      progress_detail = "detailed"
    ),
    regexp = ".*"
  )

  expect_s3_class(result, "data.table")
})
