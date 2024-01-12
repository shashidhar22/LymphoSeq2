context("Get productive sequences")
library(LymphoSeq2)
library(tidyverse)

test_that("Gather productive aminoacid sequences", {
  stable <- LymphoSeq2::readImmunoSeq(
    c(
      "test_data/015V12001549_CFAR.tsv",
      "test_data/015V12001685_CFAR_R.tsv"
    ),
    threads = 1
  )
  atable <- LymphoSeq2::productiveSeq(stable)
  arow <- base::nrow(atable)
  ctable <- tibble::tibble(
    repertoire_id = c("015V12001549_CFAR", "015V12001685_CFAR_R"),
    tot_freq = c(1, 1),
    nrows = c(802, 800),
    nnuc = c(802, 800)
  ) %>%
    dplyr::mutate(
      nrows = base::as.integer(nrows),
      nnuc = base::as.integer(nnuc)
    )
  atable <- atable %>%
    dplyr::group_by(repertoire_id) %>%
    dplyr::summarize(
      tot_freq = base::sum(duplicate_frequency),
      nrows = dplyr::n(),
      nnuc = base::length(base::unique(junction_aa))
    )
  expect_equal(arow, 1602)
  expect_true(dplyr::all_equal(ctable, atable))
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
  arow <- base::nrow(atable)
  ctable <- tibble::tibble(
    repertoire_id = c("015V12001549_CFAR", "015V12001685_CFAR_R"),
    tot_freq = c(1, 1),
    nrows = c(802, 802),
    nnuc = c(802, 800)
  ) %>%
    dplyr::mutate(
      nrows = base::as.integer(nrows),
      nnuc = base::as.integer(nnuc)
    )
  atable <- atable %>%
    dplyr::group_by(repertoire_id) %>%
    dplyr::summarize(
      tot_freq = base::sum(duplicate_frequency),
      nrows = dplyr::n(),
      nnuc = base::length(base::unique(junction_aa))
    )
  expect_equal(arow, 1604)
  expect_true(dplyr::all_equal(ctable, atable))
})


test_that("Count of collapse amino acid sequences match", {
  stable <- LymphoSeq2::readImmunoSeq("test_data/015V06013979_CFAR.tsv", threads = 1)
  atable <- LymphoSeq2::productiveSeq(stable, aggregate = "junction_aa")
  ctable <- tibble::tibble(
    repertoire_id = c("015V06013979_CFAR", "015V06013979_CFAR"),
    junction_aa = c("CASSIASAGGPDTQYF", "CASSMGQGATVGYTF"),
    duplicate_count = c(608, 184)
  )
  atable <- atable %>%
    dplyr::filter(junction_aa %in% c("CASSIASAGGPDTQYF", "CASSMGQGATVGYTF")) %>%
    dplyr::select(repertoire_id, junction_aa, duplicate_count)
  expect_true(dplyr::all_equal(ctable, atable))
})


test_that("Prevalence of amino acid sequences is correct", {
  stable <- LymphoSeq2::readImmunoSeq("test_data", threads = 1)
  ntable <- LymphoSeq2::productiveSeq(stable, aggregate = "junction_aa", prevalence = TRUE) %>%
    dplyr::select(prevalence, junction_aa) %>%
    dplyr::filter(prevalence != 0) %>%
    dplyr::arrange(junction_aa) %>%
    dplyr::distinct()
  junction_list <- ntable %>%
    dplyr::pull(junction_aa) %>%
    base::unique()
  prevalenceTRB <- LymphoSeq2::prevalenceTRB %>%
    dplyr::rename(junction_aa = "aminoAcid") %>%
    dplyr::filter(junction_aa %in% junction_list)
  expect_true(dplyr::all_equal(prevalenceTRB, ntable))
})
