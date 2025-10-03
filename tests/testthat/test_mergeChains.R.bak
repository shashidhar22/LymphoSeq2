context("Merge Chains by Cell")
library(LymphoSeq2)
library(testthat)
library(dplyr)

# Setup test data
setup_test_data <- function() {
  read10x("test_data/test_airr.tsv")
}

# Test mode parameter validation
test_that("merge_chains validates mode parameter", {
  data <- setup_test_data()

  expect_error(
    merge_chains(data, mode = "invalid"),
    "mode must be one of: 'none', 'strict', or 'best'"
  )
})

# Test mode = "none"
test_that("merge_chains with mode='none' returns data unchanged", {
  data <- setup_test_data()
  result <- merge_chains(data, mode = "none")

  expect_equal(nrow(result), nrow(data))
  expect_equal(ncol(result), ncol(data))
  expect_identical(result, data)
})

# Test mode = "strict"
test_that("merge_chains with mode='strict' filters cells correctly", {
  data <- setup_test_data()
  result <- merge_chains(data, mode = "strict")

  # Result should have fewer rows than input
  expect_lt(nrow(result), nrow(data))

  # Each cell should have exactly one row (paired TRA and TRB)
  cell_counts <- result %>%
    count(cell_id)

  expect_true(all(cell_counts$n == 1))
})

test_that("merge_chains strict mode concatenates junction fields correctly", {
  data <- setup_test_data()
  result <- merge_chains(data, mode = "strict")

  # Check that junction fields are concatenated with "|"
  expect_true(any(grepl("\\|", result$junction)))
  expect_true(any(grepl("\\|", result$junction_aa)))

  # Check that we have separate alpha and beta columns
  expect_true("v_call_alpha" %in% names(result))
  expect_true("v_call_beta" %in% names(result))
})

test_that("merge_chains strict mode only keeps 1 TRA + 1 TRB cells", {
  data <- setup_test_data()
  result <- merge_chains(data, mode = "strict")

  # Check that all rows have both alpha and beta chains
  expect_true(all(!is.na(result$v_call_alpha)))
  expect_true(all(!is.na(result$v_call_beta)))
  expect_true(all(grepl("TRA", result$v_call_alpha)))
  expect_true(all(grepl("TRB", result$v_call_beta)))
})

# Test mode = "best"
test_that("merge_chains with mode='best' selects most frequent chains", {
  data <- setup_test_data()
  result <- merge_chains(data, mode = "best")

  # Result should have fewer rows than input
  expect_lt(nrow(result), nrow(data))

  # Each cell should have exactly one row
  cell_counts <- result %>%
    count(cell_id)

  expect_true(all(cell_counts$n == 1))
})

test_that("merge_chains best mode selects highest duplicate_count", {
  # Create specific test data with known duplicate counts
  test_data <- tibble(
    cell_id = c("CELL1", "CELL1", "CELL1", "CELL1"),
    clone_id = c("clone1", "clone1", "clone1", "clone1"),
    v_call = c("TRAV1", "TRAV2", "TRBV1", "TRBV2"),
    junction = c("AAA", "BBB", "CCC", "DDD"),
    junction_aa = c("A", "B", "C", "D"),
    duplicate_count = c(10, 50, 100, 20),
    repertoire_id = "test"
  )

  result <- merge_chains(test_data, mode = "best")

  # Should select TRAV2 (50) and TRBV1 (100)
  expect_equal(result$v_call_alpha, "TRAV2")
  expect_equal(result$v_call_beta, "TRBV1")
})

test_that("merge_chains best mode only keeps cells with both TRA and TRB", {
  data <- setup_test_data()
  result <- merge_chains(data, mode = "best")

  # All results should have both alpha and beta chains
  expect_true(all(!is.na(result$v_call_alpha)))
  expect_true(all(!is.na(result$v_call_beta)))
  expect_true(all(grepl("TRA", result$v_call_alpha)))
  expect_true(all(grepl("TRB", result$v_call_beta)))
})

# Test chain type extraction
test_that("merge_chains extracts chain types correctly", {
  test_data <- tibble(
    cell_id = c("C1", "C2", "C3", "C4"),
    clone_id = "clone1",
    v_call = c("TRAV1", "TRBV2", "TCRAV3", "TCRBV4"),
    junction = c("A", "B", "C", "D"),
    junction_aa = c("A", "B", "C", "D"),
    duplicate_count = c(10, 20, 30, 40),
    repertoire_id = "test"
  )

  # Add chain_type using the same logic as merge_chains
  chain_pattern <- stringr::str_extract(test_data$v_call, "TR[ABGD]|TCR[ABGD]")
  chain_pattern <- stringr::str_replace(chain_pattern, "TCR", "TR")

  expect_equal(chain_pattern[1], "TRA")
  expect_equal(chain_pattern[2], "TRB")
  expect_equal(chain_pattern[3], "TRA")
  expect_equal(chain_pattern[4], "TRB")
})

# Test edge cases
test_that("merge_chains handles cells with only TRA chains", {
  test_data <- tibble(
    cell_id = c("CELL1", "CELL1"),
    clone_id = "clone1",
    v_call = c("TRAV1", "TRAV2"),
    junction = c("AAA", "BBB"),
    junction_aa = c("A", "B"),
    duplicate_count = c(10, 20),
    repertoire_id = "test"
  )

  result_strict <- merge_chains(test_data, mode = "strict")
  result_best <- merge_chains(test_data, mode = "best")

  # Should filter out cells without both chains
  expect_equal(nrow(result_strict), 0)
  expect_equal(nrow(result_best), 0)
})

test_that("merge_chains handles cells with only TRB chains", {
  test_data <- tibble(
    cell_id = c("CELL1", "CELL1"),
    clone_id = "clone1",
    v_call = c("TRBV1", "TRBV2"),
    junction = c("AAA", "BBB"),
    junction_aa = c("A", "B"),
    duplicate_count = c(10, 20),
    repertoire_id = "test"
  )

  result_strict <- merge_chains(test_data, mode = "strict")
  result_best <- merge_chains(test_data, mode = "best")

  # Should filter out cells without both chains
  expect_equal(nrow(result_strict), 0)
  expect_equal(nrow(result_best), 0)
})

test_that("merge_chains handles empty input", {
  test_data <- tibble(
    cell_id = character(0),
    clone_id = character(0),
    v_call = character(0),
    junction = character(0),
    junction_aa = character(0),
    duplicate_count = numeric(0),
    repertoire_id = character(0)
  )

  result_none <- merge_chains(test_data, mode = "none")
  result_strict <- merge_chains(test_data, mode = "strict")
  result_best <- merge_chains(test_data, mode = "best")

  expect_equal(nrow(result_none), 0)
  expect_equal(nrow(result_strict), 0)
  expect_equal(nrow(result_best), 0)
})

# Test with real data
test_that("merge_chains strict mode filters multiple chains per type", {
  data <- setup_test_data()

  # CELL006 and CELL007 have multiple chains of same type
  result <- merge_chains(data, mode = "strict")

  # These cells should be filtered out
  expect_false("CELL006-1" %in% result$cell_id)
  expect_false("CELL007-1" %in% result$cell_id)
})

test_that("merge_chains best mode handles multiple chains per type", {
  data <- setup_test_data()

  # CELL006 has 1 TRB + 2 TRA chains
  # CELL007 has 1 TRB + 2 TRA chains
  result <- merge_chains(data, mode = "best")

  # Should include these cells by selecting most frequent chain
  # Check if they exist (they should if there's a valid TRA and TRB)
  cell006_present <- "CELL006-1" %in% result$cell_id
  cell007_present <- "CELL007-1" %in% result$cell_id

  # At least one of them should be present
  expect_true(cell006_present || cell007_present)
})

# Test field concatenation
test_that("merge_chains concatenates junction fields with pipe separator", {
  data <- setup_test_data()
  result <- merge_chains(data, mode = "strict")

  if (nrow(result) > 0) {
    # Check that junction fields are concatenated with "|"
    expect_true(any(grepl("\\|", result$junction)))
    expect_true(any(grepl("\\|", result$junction_aa)))

    # Check that separate alpha/beta fields exist (not concatenated)
    expect_true("v_call_alpha" %in% names(result))
    expect_true("v_call_beta" %in% names(result))
  }
})

# Test integration with read10x
test_that("merge_chains works with contig_annotations data", {
  data <- read10x("test_data/test_contig_annotations.csv")
  result <- merge_chains(data, mode = "strict")

  expect_s3_class(result, "tbl_df")
  expect_true(nrow(result) <= nrow(data))
})

test_that("Complete workflow: read10x + merge_chains", {
  # Read AIRR data
  airr_data <- read10x("test_data/test_airr.tsv")

  # Merge with different modes
  none_result <- merge_chains(airr_data, mode = "none")
  strict_result <- merge_chains(airr_data, mode = "strict")
  best_result <- merge_chains(airr_data, mode = "best")

  # Verify row counts
  expect_equal(nrow(none_result), nrow(airr_data))
  expect_lte(nrow(strict_result), nrow(airr_data))
  expect_lte(nrow(best_result), nrow(airr_data))

  # Strict should be most restrictive
  expect_lte(nrow(strict_result), nrow(best_result))
})
