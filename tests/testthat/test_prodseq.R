context("Get productive sequences")
library(LymphoSeq2)

test_that("Gather productive aminoacid sequences", {
  stable <- LymphoSeq2::readImmunoSeq(
    c(
      "test_data/015V12001549_CFAR.tsv",
      "test_data/015V12001685_CFAR_R.tsv"
    ),
    threads = 1
  )
  atable <- LymphoSeq2::productiveSeq(stable, aggregate = "junction_aa")

  # Test 1: All amino acid sequences should be unique
  expect_equal(nrow(atable), length(unique(atable$junction_aa)))

  # Test 2: No stop codons should be present
  expect_true(all(!stringr::str_detect(atable$junction_aa, "\\*")))

  # Test 3: All sequences should be in-frame
  expect_true(all(atable$reading_frame == "in-frame"))

  # Test 4: Frequencies should sum to 1 per repertoire
  freq_sums <- atable %>%
    dplyr::group_by(repertoire_id) %>%
    dplyr::summarize(total_freq = sum(duplicate_frequency), .groups = "drop")
  expect_true(all(abs(freq_sums$total_freq - 1) < 1e-10))
})

test_that("Gather productive nucleotide sequences", {
  stable <- LymphoSeq2::readImmunoSeq(
    c(
      "test_data/015V12001549_CFAR.tsv",
      "test_data/015V12001685_CFAR_R.tsv"
    ),
    threads = 1
  )
  atable <- LymphoSeq2::productiveSeq(stable, aggregate = "junction")

  # Test 1: All nucleotide sequences should be unique
  expect_equal(nrow(atable), length(unique(atable$junction)))

  # Test 2: No stop codons should be present in nucleotide sequences
  expect_true(all(!stringr::str_detect(atable$junction, "\\*")))

  # Test 3: All sequences should be in-frame
  expect_true(all(atable$reading_frame == "in-frame"))

  # Test 4: Frequencies should sum to 1 per repertoire
  freq_sums <- atable %>%
    dplyr::group_by(repertoire_id) %>%
    dplyr::summarize(total_freq = sum(duplicate_frequency), .groups = "drop")
  expect_true(all(abs(freq_sums$total_freq - 1) < 1e-10))
})


test_that("Count of collapse amino acid sequences match", {
  stable <- LymphoSeq2::readImmunoSeq("test_data/015V06013979_CFAR.tsv", threads = 1)
  atable <- LymphoSeq2::productiveSeq(stable, aggregate = "junction_aa")
  atable_filtered <- atable %>%
    dplyr::filter(junction_aa %in% c("CASSIASAGGPDTQYF", "CASSMGQGATVGYTF")) %>%
    dplyr::select(repertoire_id, junction_aa, duplicate_count)

  # Check that expected sequences exist with correct counts
  expect_equal(nrow(atable_filtered), 2)
  expect_equal(sort(atable_filtered$duplicate_count), c(184, 608))
  expect_true(all(c("CASSIASAGGPDTQYF", "CASSMGQGATVGYTF") %in% atable_filtered$junction_aa))
})


test_that("Prevalence of amino acid sequences is correct", {
  skip_if_not_installed("LymphoSeqDB")

  stable <- LymphoSeq2::readImmunoSeq("test_data", threads = 1) |>
    dplyr::filter(stringr::str_starts(repertoire_id, "015V"))
  ntable <- LymphoSeq2::productiveSeq(stable, aggregate = "junction_aa", prevalence = TRUE) %>%
    dplyr::select(prevalence, junction_aa) %>%
    dplyr::filter(prevalence != 0) %>%
    dplyr::arrange(junction_aa) %>%
    dplyr::distinct()
  junction_list <- ntable %>%
    dplyr::pull(junction_aa) %>%
    base::unique()
  prevalenceTRB <- LymphoSeqDB::prevalenceTRB %>%
    dplyr::rename(junction_aa = "aminoAcid") %>%
    dplyr::filter(junction_aa %in% junction_list) %>%
    dplyr::arrange(junction_aa)
  expect_true(isTRUE(all.equal(prevalenceTRB, ntable, check.attributes = FALSE)))
})
