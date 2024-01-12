context("Check if summary statistics for dataset are correct")
library(LymphoSeq2)
library(tidyverse)

test_that("Check if summary statistics for test data are correct", {
  stable <- LymphoSeq2::readImmunoSeq("test_data/015V06013979_CFAR.tsv", threads = 1)
  ctable <- LymphoSeq2::clonality(stable)
  ttseq <- ctable %>%
    dplyr::pull(total_sequences)
  tupseq <- ctable %>%
    dplyr::pull(unique_productive_sequences)
  ttcount <- ctable %>%
    dplyr::pull(total_count)
  tclonality <- ctable %>%
    dplyr::pull(clonality)
  tcon <- ctable %>%
    dplyr::pull(convergence)
  tgc <- ctable %>%
    dplyr::pull(gini_coefficient)
  ttps <- ctable %>%
    dplyr::pull(top_productive_sequence)
  expect_equal(ttseq, 1000)
  expect_equal(tupseq, 836)
  expect_equal(ttcount, 2404)
  expect_equal(base::round(tclonality, 3), 0.325)
  expect_equal(base::round(tgc, 3), 0.599)
  expect_equal(base::round(ttps, 3), 28.639)
  expect_equal(tcon, 1)
})
