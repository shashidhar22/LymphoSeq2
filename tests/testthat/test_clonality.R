context("Check if summary statistics for dataset are correct")
library(LymphoSeq2)


test_that("Check if summary statistics for test data are correct", {
  stable <- LymphoSeq2::readImmunoSeq("test_data/015V06013979_CFAR.tsv")
  ctable <- LymphoSeq2::clonality(stable)
  ttseq <- ctable %>%
           dplyr::pull(total_sequences)
  tupseq <- ctable %>%
            dplyr::pull(unique_productive_sequences)
  ttcount <- ctable %>% 
             dplyr::pull(total_count)
  tclonality <- ctable %>%
                dplyr::pull(clonality)
  tsi <- ctable %>%
         dplyr::pull(simpson_index)
  tis <- ctable %>%
         dplyr::pull(inverse_simpson)
  tgc <- ctable %>% 
         dplyr::pull(gini_coefficient)
  ttps <- ctable %>%
          dplyr::pull(top_productive_sequence)
  tce <- ctable %>%
         dplyr::pull(chao_estimate)
  tke <- ctable %>% 
         dplyr::pull(kemp_estimate)
  the <- ctable %>%
         dplyr::pull(hill_estimate)
  expect_equal(ttseq, 267)
  expect_equal(tupseq, 267)
  expect_equal(ttcount, 16645)
  expect_equal(base::round(tclonality, 3), 0.388)
  expect_equal(base::round(tsi, 3), 0.937)
  expect_equal(base::round(tis, 3), 15.912)
  expect_equal(base::round(tgc, 3), 0.901)
  expect_equal(base::round(ttps, 3), 18.522)
  expect_equal(base::round(tce, 3), 316.390)
  expect_equal(base::round(tke, 3), 411.088)
  expect_equal(the, 267)
})
