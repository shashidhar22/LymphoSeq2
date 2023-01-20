context("Count k-mers in nucleotide sequence")
library(LymphoSeq2)
library(tidyverse)

test_that("Number of k-mers are counted correctly", {
    junction <- "ATCGATCAC"
    study_table <- LymphoSeq2::readImmunoSeq("test_data/015V06013979_CFAR.tsv")
    ktable <- LymphoSeq2::countKmer(study_table = study_table, k = 3)
    num_rows <- base::nrow(ktable)
    rep_id <- study_table %>% 
      dplyr::pull(repertoire_id) %>% 
      unique()
    Kmer <- c('AAA', 'AAC', 'AAG', 'AAT', 'ACA', 'ACC', 'ACG', 'ACT', 'AGA', 
              'AGC', 'AGG', 'AGT', 'ATA', 'ATC', 'ATG', 'ATT', 'CAA', 'CAC', 
              'CAG', 'CAT', 'CCA', 'CCC', 'CCG', 'CCT', 'CGA', 'CGC', 'CGG', 
              'CGT', 'CTA', 'CTC', 'CTG', 'CTT', 'GAA', 'GAC', 'GAG', 'GAT', 
              'GCA', 'GCC', 'GCG', 'GCT', 'GGA', 'GGC', 'GGG', 'GGT', 'GTA', 
              'GTC', 'GTG', 'GTT', 'TAA', 'TAC', 'TAG', 'TAT', 'TCA', 'TCC', 
              'TCG', 'TCT', 'TGA', 'TGC', 'TGG', 'TGT', 'TTA', 'TTC', 'TTG', 
              'TTT')
    Count <- c(558, 761, 762, 328, 1750, 1090, 503, 1626, 1300, 3260, 1186, 
               1316, 362, 1207, 457, 592, 755, 1078, 4335, 794, 2623, 1890, 
               796, 1402, 529, 612, 1734, 231, 868, 1986, 2630, 1957, 1036, 
               1802, 1514, 489, 2034, 2685, 691, 1877, 2030, 1784, 3380, 661, 
               1289, 441, 1848, 665, 194, 1275, 424, 1018, 694, 1037, 1379, 
               2411, 902, 1566, 1534, 2006, 390, 1816, 1044, 1806)
    kmer_table <- tibble::tibble(Kmer, Count) 
    kmer_table <- magrittr::set_colnames(kmer_table, c("Kmer", rep_id))
    expect_equal(num_rows, 64)
    expect_equal(ktable, kmer_table)
})
