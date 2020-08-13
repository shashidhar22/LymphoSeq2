context("Read Immuno Seq")
library(LymphoSeq2)

test_that("Reads a single AIRR file correctly", {
  stable <- LymphoSeq2::readImmunoSeq("test_data/015V12001549_CFAR.tsv")
  snuc <- base::length(base::unique(stable$junction))
  srow <- base::nrow(stable)
  tot_freq <- stable %>% dplyr::pull(duplicate_frequency) %>% base::sum()
  expect_equal(base::nrow(stable), 50)
  expect_equal(snuc, srow)
  expect_equal(tot_freq, 1)
})

test_that("Reads a list of AIRR file correctly", {
  stable <- LymphoSeq2::readImmunoSeq(c("test_data/015V12001549_CFAR.tsv", "test_data/015V12001685_CFAR_R.tsv"))
  srow <- base::nrow(stable)
  ctable <- tibble::tibble(repertoire_id = c("015V12001549_CFAR", "015V12001685_CFAR_R"),
                           tot_freq = c(1, 1),
                           nrows = c(50, 50),
                           nnuc = c(50, 50)) %>%
            dplyr::mutate(nrows = base::as.integer(nrows),
                   nnuc = base::as.integer(nnuc))
  stable <- stable %>% 
            dplyr::group_by(repertoire_id) %>% 
            dplyr::summarize(tot_freq = base::sum(duplicate_frequency), 
                             nrows = dplyr::n(), 
                             nnuc = base::length(base::unique(junction)))
  expect_equal(srow, 100)
  expect_true(dplyr::all_equal(ctable, stable))
})

test_that("Reads AIRR files from a path correctly", {
  stable <- LymphoSeq2::readImmunoSeq("test_data/")
  srow <- base::nrow(stable)
  ctable <- tibble::tibble(repertoire_id = c("015V12001549_CFAR", "015V12001685_CFAR_R", "015V12003105_CFAR", "015V06013979_CFAR"),
                           tot_freq = c(1, 1, 1, 1),
                           nrows = c(50, 50, 50, 267),
                           nnuc = c(50, 50, 50, 267)) %>%
             dplyr::mutate(nrows = base::as.integer(nrows),
                           nnuc = base::as.integer(nnuc))
  stable <- stable %>% 
            dplyr::group_by(repertoire_id) %>% 
            dplyr::summarize(tot_freq = base::sum(duplicate_frequency), 
                             nrows = dplyr::n(), 
                             nnuc = base::length(base::unique(junction)))
  expect_equal(srow, 417)
  expect_true(dplyr::all_equal(ctable, stable))
})
