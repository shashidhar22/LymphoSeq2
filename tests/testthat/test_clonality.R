context("Check if summary statistics for dataset are correct")
library(LymphoSeq2)

test_that("Check if summary statistics for test data are correct", {
  stable <- LymphoSeq2::readImmunoSeq("test_data/015V06013979_CFAR.tsv", threads = 1)
  ctable <- LymphoSeq2::clonality(stable)
  ttseq <- ctable |>
    dplyr::pull(total_sequences)
  tupseq <- ctable |>
    dplyr::pull(unique_productive_sequences)
  ttcount <- ctable |>
    dplyr::pull(total_count)
  tclonality <- ctable |>
    dplyr::pull(clonality)
  tcon <- ctable |>
    dplyr::pull(convergence)
  tgc <- ctable |>
    dplyr::pull(gini_coefficient)
  tsimpson <- ctable |>
    dplyr::pull(simpson_index)
  tinv_simpson <- ctable |>
    dplyr::pull(inverse_simpson)
  ttps <- ctable |>
    dplyr::pull(top_productive_sequence)

  # Basic counts
  expect_equal(ttseq, 1000)
  expect_equal(tupseq, 846)
  expect_equal(ttcount, 2404)

  # Diversity metrics
  expect_equal(base::round(tclonality, 3), 0.323)
  expect_equal(base::round(tgc, 3), 0.605)
  expect_equal(base::round(tsimpson, 3), 0.091)
  expect_equal(base::round(tinv_simpson, 1), 11.0)

  # Other metrics
  expect_equal(base::round(ttps, 3), 27.813)
  expect_equal(tcon, 1)
})
