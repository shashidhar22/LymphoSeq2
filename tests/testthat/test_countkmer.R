context("Count k-mers in nucleotide sequence")
library(LymphoSeq2)

test_that("Number of k-mers are counted correctly", {
    junction <- "ATCGATCAC"
    study_table <- tibble(junction)
    ktable <- LymphoSeq2::countKmer(study_table = study_table, k = 3)
    num_rows <- base::nrow(ktable)
    Kmer <- c("ATC", "TCG", "CGA", "GAT", "TCA", "CAC")
    Count <- c(2, 1, 1, 1, 1, 1)
    kmer_table <- tibble(Kmer, Count)
    expect_equal(num_rows, 6)
    expect_equal(ktable, kmer_table)
})