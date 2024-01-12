context("Calculate clonal realtedness of samples")
library(LymphoSeq2)
library(tidyverse)

test_that("Calculate clonal relatedness of all sequences", {
  stable <- LymphoSeq2::readImmunoSeq("test_data/", threads = 1)
  ttable <- LymphoSeq2::clonalRelatedness(stable)
  ctable <- tibble::tibble(
    repertoire_id = c("015V06013979_CFAR", "015V12001549_CFAR", "015V12001685_CFAR_R", "015V12003105_CFAR"),
    relatedness = c(0.001, 0.001, 0.002, 0.001)
  )
  expect_true(base::all.equal(ttable, ctable))
})

test_that("Calculate clonal relatedness of all productive nucleotide sequences", {
  stable <- LymphoSeq2::readImmunoSeq("test_data/", threads = 1)
  ntable <- LymphoSeq2::productiveSeq(stable, aggregate = "junction")
  ttable <- LymphoSeq2::clonalRelatedness(ntable) %>%
    dplyr::mutate(relatedness = base::round(relatedness, 3))
  ctable <- tibble::tibble(
    repertoire_id = c("015V06013979_CFAR", "015V12001549_CFAR", "015V12001685_CFAR_R", "015V12003105_CFAR"),
    relatedness = c(0.001, 0.001, 0.002, 0.001)
  )
  expect_true(base::all.equal(ttable, ctable))
})
