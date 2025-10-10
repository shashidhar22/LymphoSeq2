context("Get VDJ gene frequency")
library(LymphoSeq2)

test_that("Gene frequency charcterized by gene name sums to 1", {
  stable <- LymphoSeq2::readImmunoSeq("test_data", threads = 1) |>
    dplyr::filter(stringr::str_starts(repertoire_id, "015V"))
  ntable <- LymphoSeq2::productiveSeq(stable, aggregate = "junction")
  gtable <- LymphoSeq2::geneFreq(ntable) |>
    dplyr::group_by(repertoire_id, gene_type) |>
    dplyr::summarise(freq_tot = sum(gene_frequency), .groups = "drop") |>
    dplyr::pull(freq_tot)
  expect_true(base::all(abs(gtable - 1) < 1e-10))
})

test_that("Gene frequency charcterized by gene name sums to 1", {
  stable <- LymphoSeq2::readImmunoSeq("test_data", threads = 1) |>
    dplyr::filter(stringr::str_starts(repertoire_id, "015V"))
  ntable <- LymphoSeq2::productiveSeq(stable, aggregate = "junction")
  gtable <- LymphoSeq2::geneFreq(ntable, family = TRUE) |>
    dplyr::group_by(repertoire_id, gene_type) |>
    dplyr::summarise(freq_tot = sum(gene_frequency), .groups = "drop") |>
    dplyr::pull(freq_tot)
  expect_true(base::all(abs(gtable - 1) < 1e-10))
})
