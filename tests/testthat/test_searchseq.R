context("Search of similar sequences from table")
library(LymphoSeq2)
library(tidyverse)

test_that("Find all nucleotide sequences with edit distance ten", {
  stable <- LymphoSeq2::readImmunoSeq("test_data/015V06013979_CFAR.tsv", threads = 1)
  ntable <- LymphoSeq2::productiveSeq(stable, aggregate = "junction")
  seq <- "TTGGAGCTGGGGGACTCGGCCCTTTATCTTTGCGCCAGCAGCTCCGGGACAGGGGGCTCGGGCAATCAGCCCCAGCATTTTGGTGAT"
  nlist <- ntable %>%
           dplyr::pull(junction) %>%
           base::unique()
  edist <- utils::adist(seq, nlist, partial = FALSE)
  elist <- nlist[(edist <= 10)]
  slist <- LymphoSeq2::searchSeq(ntable, seq, edit_distance = 10) %>%
           dplyr::pull(junction) %>%
           base::unique()
  expect_equal(slist, elist)
})

test_that("Find all amino acid sequences with edit distance ten from list of reference", {
  stable <- LymphoSeq2::readImmunoSeq("test_data/015V06013979_CFAR.tsv", threads = 1)
  ntable <- LymphoSeq2::productiveSeq(stable, aggregate = "junction_aa")
  seq <- "CASSSGTGGSGNQPQHF"
  nlist <- ntable %>%
    dplyr::pull(junction_aa) %>%
    base::unique()
  edist <- utils::adist(seq, nlist, partial = FALSE)
  elist <- nlist[(edist <= 10)]
  slist <- LymphoSeq2::searchSeq(ntable, seq, edit_distance = 10, seq_type = "junction_aa") %>%
    dplyr::pull(junction_aa) %>%
    base::unique()
  expect_equal(slist, elist)
})
