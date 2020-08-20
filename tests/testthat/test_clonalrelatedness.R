context("Calculate clonal realtedness of samples")
library(LymphoSeq2)

test_that("Calculate clonal relatedness of all sequences", {
  stable <- LymphoSeq2::readImmunoSeq("test_data/")
  ttable <- LymphoSeq2::clonalRelatedness(stable)
  ctable <- tibble::tibble(repertoire_id = c("015V06013979_CFAR", "015V12001549_CFAR", "015V12001685_CFAR_R", "015V12003105_CFAR"),
                           clonalRelatedness = c(0.02247191, 0.02000000, 0.02000000, 0.02000000))
  expect_true(base::all.equal(ttable, ctable))
})

test_that("Calculate clonal relatedness of all productive nucleotide sequences", {
  stable <- LymphoSeq2::readImmunoSeq("test_data/")
  ntable <- LymphoSeq2::productiveSeq(stable, aggregate = "junction")
  ttable <- LymphoSeq2::clonalRelatedness(ntable) %>%
            dplyr::mutate(clonalRelatedness = base::round(clonalRelatedness, 3))
  ctable <- tibble::tibble(repertoire_id = c("015V06013979_CFAR", "015V12001549_CFAR", "015V12001685_CFAR_R", "015V12003105_CFAR"),
    clonalRelatedness = c(0.022, 0.025, 0.026, 0.027))
  expect_true(base::all.equal(ttable, ctable))
})
