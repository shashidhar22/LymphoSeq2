context("Read 10x Genomics Files")
library(LymphoSeq2)
library(testthat)
library(dplyr)

# Test by reading files with read10x and checking results
# (Internal functions are tested indirectly through read10x)

# Test read10x main function
test_that("read10x reads single AIRR file", {
  data <- read10x("test_data/test_airr.tsv")

  expect_s3_class(data, "tbl_df")
  expect_equal(nrow(data), 16)
  expect_equal(length(unique(data$cell_id)), 7)
})

test_that("read10x reads single contig_annotations file", {
  data <- read10x("test_data/test_contig_annotations.csv")

  expect_s3_class(data, "tbl_df")
  expect_equal(nrow(data), 6)
  expect_equal(length(unique(data$cell_id)), 3)
})

test_that("read10x reads multiple files", {
  data <- read10x(c(
    "test_data/test_airr.tsv",
    "test_data/test_contig_annotations.csv"
  ))

  expect_s3_class(data, "tbl_df")
  expect_equal(nrow(data), 22)
  # Both files may have same base name, so just check we get data
  expect_true(length(unique(data$repertoire_id)) >= 1)
})

test_that("read10x errors on empty file", {
  # Create empty file for test
  temp_file <- tempfile(fileext = ".tsv")
  file.create(temp_file)

  expect_error(
    read10x(temp_file),
    "No valid files"
  )

  unlink(temp_file)
})

# Test chain counting
test_that("AIRR file has correct chain distribution", {
  data <- read10x("test_data/test_airr.tsv")

  chain_counts <- data |>
    mutate(chain = ifelse(grepl("TRA", v_call), "TRA", "TRB")) |>
    count(cell_id, chain) |>
    tidyr::pivot_wider(names_from = chain, values_from = n, values_fill = 0)

  # Check cells with exactly 1 TRA and 1 TRB
  single_pair_cells <- chain_counts |>
    filter(TRA == 1, TRB == 1) |>
    nrow()

  expect_true(single_pair_cells >= 1)
})

test_that("Contig annotations file has correct chain types", {
  data <- read10x("test_data/test_contig_annotations.csv")

  chain_types <- unique(stringr::str_extract(data$v_call, "TR[AB]"))
  expect_true(all(chain_types %in% c("TRA", "TRB")))
})

# Test junction lengths
test_that("Junction lengths are calculated correctly", {
  data <- read10x("test_data/test_contig_annotations.csv")

  calculated_lengths <- stringr::str_length(data$junction)
  expect_equal(data$junction_length, calculated_lengths)

  calculated_aa_lengths <- stringr::str_length(data$junction_aa)
  expect_equal(data$junction_aa_length, calculated_aa_lengths)
})

# Test repertoire_id extraction
test_that("Repertoire ID is extracted from AIRR filename", {
  data <- read10x("test_data/test_airr.tsv")
  # Check that repertoire_id exists and is not empty
  expect_true(nchar(unique(data$repertoire_id)) > 0)
  expect_true("test" %in% unique(data$repertoire_id))
})

test_that("Repertoire ID is extracted from contig_annotations filename", {
  data <- read10x("test_data/test_contig_annotations.csv")
  # Check that repertoire_id exists and is not empty
  expect_true(nchar(unique(data$repertoire_id)) > 0)
  expect_true("test" %in% unique(data$repertoire_id))
})
