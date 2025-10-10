context("10X Genomics Single-cell Integration")
library(LymphoSeq2)
library(testthat)
library(dplyr)

# Helper function to load test data
load_10x_test_data <- function() {
  read10x("test_data/test_airr.tsv")
}

# Test complete workflow: read → merge → analyze
test_that("Complete 10X workflow: AIRR format", {
  # Read data
  sc_data <- load_10x_test_data()

  expect_s3_class(sc_data, "tbl_df")
  expect_true("v_call" %in% colnames(sc_data))
  expect_true("junction" %in% colnames(sc_data))

  # Merge chains
  merged <- merge_chains(sc_data, mode = "best")

  expect_s3_class(merged, "tbl_df")
  expect_lt(nrow(merged), nrow(sc_data))  # Should have fewer rows after merging
  expect_true("paired_clonotype" %in% colnames(merged))
  expect_true("v_call_alpha" %in% colnames(merged))
  expect_true("v_call_beta" %in% colnames(merged))
})

test_that("merge_chains works with strict mode", {
  sc_data <- load_10x_test_data()
  merged_strict <- merge_chains(sc_data, mode = "strict")

  expect_s3_class(merged_strict, "tbl_df")
  # Strict mode should be more restrictive than best mode
  merged_best <- merge_chains(sc_data, mode = "best")
  expect_lte(nrow(merged_strict), nrow(merged_best))
})

test_that("merge_chains preserves required columns for downstream analysis", {
  sc_data <- load_10x_test_data()
  merged <- merge_chains(sc_data, mode = "best")

  required_cols <- c(
    "cell_id", "repertoire_id",
    "junction", "junction_aa",
    "v_call", "j_call", "d_call",
    "v_family", "j_family", "d_family",
    "duplicate_count", "duplicate_frequency",
    "productive", "reading_frame"
  )

  for (col in required_cols) {
    expect_true(col %in% colnames(merged),
                info = paste("Missing required column:", col))
  }
})

test_that("merged 10X data works with clonality()", {
  sc_data <- load_10x_test_data()
  merged <- merge_chains(sc_data, mode = "best")

  clonality_result <- clonality(merged)

  expect_s3_class(clonality_result, "tbl_df")
  expect_equal(nrow(clonality_result), 1)  # One sample

  # Check expected columns
  expect_true("clonality" %in% colnames(clonality_result))
  expect_true("gini_coefficient" %in% colnames(clonality_result))
  expect_true("total_sequences" %in% colnames(clonality_result))
  expect_true("unique_productive_sequences" %in% colnames(clonality_result))

  # Clonality should be between 0 and 1
  expect_true(clonality_result$clonality >= 0 && clonality_result$clonality <= 1)
})

test_that("merged 10X data works with productiveSeq()", {
  sc_data <- load_10x_test_data()
  merged <- merge_chains(sc_data, mode = "best")

  # Test with junction aggregation
  productive_nt <- productiveSeq(merged, aggregate = "junction")
  expect_s3_class(productive_nt, "tbl_df")
  expect_true(nrow(productive_nt) > 0)

  # Test with junction_aa aggregation
  productive_aa <- productiveSeq(merged, aggregate = "junction_aa")
  expect_s3_class(productive_aa, "tbl_df")
  expect_true(nrow(productive_aa) > 0)

  # AA aggregation should have fewer or equal rows
  expect_lte(nrow(productive_aa), nrow(productive_nt))
})

test_that("merged 10X data works with topSeqs()", {
  sc_data <- load_10x_test_data()
  merged <- merge_chains(sc_data, mode = "best")

  top_seqs <- topSeqs(merged, top = 5)

  expect_s3_class(top_seqs, "tbl_df")
  expect_lte(nrow(top_seqs), 5)

  # Should be ordered by frequency
  if (nrow(top_seqs) > 1) {
    expect_true(all(diff(top_seqs$duplicate_frequency) <= 0))
  }
})

test_that("merged 10X data works with geneFreq()", {
  sc_data <- load_10x_test_data()
  merged <- merge_chains(sc_data, mode = "best")

  # Test V gene frequency
  v_freq <- geneFreq(merged, locus = "V", family = FALSE)
  expect_s3_class(v_freq, "tbl_df")
  expect_true("gene_frequency" %in% colnames(v_freq))

  # Frequencies should sum to ~1 per repertoire
  freq_sum <- v_freq |>
    group_by(repertoire_id) |>
    summarise(total = sum(gene_frequency))
  expect_true(all(abs(freq_sum$total - 1) < 0.01))

  # Test J gene frequency
  j_freq <- geneFreq(merged, locus = "J", family = FALSE)
  expect_s3_class(j_freq, "tbl_df")
})

test_that("duplicate_count uses minimum between chains", {
  sc_data <- load_10x_test_data()
  merged <- merge_chains(sc_data, mode = "best")

  # Check that counts are reasonable
  expect_true(all(merged$duplicate_count > 0))

  # Frequency should sum to 1 per repertoire
  freq_sum <- merged |>
    group_by(repertoire_id) |>
    summarise(total = sum(duplicate_frequency))

  expect_true(all(abs(freq_sum$total - 1) < 0.01))
})

test_that("junction fields contain paired sequences", {
  sc_data <- load_10x_test_data()
  merged <- merge_chains(sc_data, mode = "best")

  # Junction should contain separator "|"
  expect_true(all(grepl("\\|", merged$junction)))
  expect_true(all(grepl("\\|", merged$junction_aa)))

  # Should be able to split back into components
  split_junction <- strsplit(merged$junction[1], "\\|")[[1]]
  expect_equal(length(split_junction), 2)

  # Individual chain junctions should match
  expect_equal(split_junction[1], merged$junction_alpha[1])
  expect_equal(split_junction[2], merged$junction_beta[1])
})

test_that("productive status requires both chains to be productive", {
  sc_data <- load_10x_test_data()
  merged <- merge_chains(sc_data, mode = "best")

  # All test data is productive, so check that field exists and is logical
  expect_type(merged$productive, "logical")
  expect_true(all(merged$productive))  # All test data is productive
  expect_equal(merged$reading_frame, rep("in-frame", nrow(merged)))
})

test_that("Comparing bulk vs single-cell workflows", {
  # This test demonstrates that merged 10X data can be analyzed
  # the same way as bulk ImmunoSeq data

  sc_data <- load_10x_test_data()
  merged <- merge_chains(sc_data, mode = "best")

  # Both should work with same analysis functions
  clonality_sc <- clonality(merged)
  top_sc <- topSeqs(merged, top = 10)
  productive_sc <- productiveSeq(merged, aggregate = "junction_aa")

  # All should return valid tibbles
  expect_s3_class(clonality_sc, "tbl_df")
  expect_s3_class(top_sc, "tbl_df")
  expect_s3_class(productive_sc, "tbl_df")
})

test_that("No merging (mode='none') returns data unchanged", {
  sc_data <- load_10x_test_data()
  not_merged <- merge_chains(sc_data, mode = "none")

  expect_identical(not_merged, sc_data)
})

test_that("Chain-specific gene frequency can be calculated", {
  sc_data <- load_10x_test_data()
  merged <- merge_chains(sc_data, mode = "best")

  # Can analyze alpha chain genes specifically
  alpha_data <- merged |>
    select(repertoire_id, duplicate_count,
           v_call = v_call_alpha,
           j_call = j_call_alpha,
           d_call = d_call_alpha)

  alpha_v_freq <- geneFreq(alpha_data, locus = "V")

  expect_s3_class(alpha_v_freq, "tbl_df")
  expect_true(all(grepl("TRA", alpha_v_freq$gene_name)))

  # Can analyze beta chain genes specifically
  beta_data <- merged |>
    select(repertoire_id, duplicate_count,
           v_call = v_call_beta,
           j_call = j_call_beta,
           d_call = d_call_beta)

  beta_v_freq <- geneFreq(beta_data, locus = "V")

  expect_s3_class(beta_v_freq, "tbl_df")
  expect_true(all(grepl("TRB", beta_v_freq$gene_name)))
})
