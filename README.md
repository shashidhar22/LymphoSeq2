
<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- Logo was obtained from this public domain link https://www.publicdomainpictures.net/en/view-image.php?image=149252&picture=cat-stretch-black-silhouette -->

# LymphoSeq2 <a href="https://shashidhar22.github.io/LymphoSeq2"><img src="man/figures/logo.png" align="right" height="138" alt="LymphoSeq2 website" /></a>

<!-- badges: start -->

[![R-CMD-check](https://github.com/shashidhar22/LymphoSeq2/workflows/R-CMD-check/badge.svg)](https://github.com/shashidhar22/LymphoSeq2/actions)
<!-- badges: end -->

## Overview

Adaptive Immune Receptor Repertoire Sequencing (AIRR-seq) provides a
unique opportunity to interrogate the adaptive immune repertoire under
various clinical conditions. The utility offered by this technology has
quickly garnered interest from a community of clinicians and researchers
investigating the immunological landscapes of a large spectrum of health
and disease states. LymphoSeq2 is a toolkit that allows users to import,
manipulate and visualize AIRR-Seq data from various AIRR-Seq assays such
as Adaptive ImmunoSEQ, BGI-IRSeq, and 10X VDJ sequencing. The platform also supports the importing of AIRR-seq data
processed using the MiXCR pipeline. The vignette highlights some of the
key features of LymphoSeq2.

## Installation

To install the latest version of LymphoSeq2 you can use the `devtools`
package and install LymphoSeq2 from GitHub

``` r
# install.packages("devtools")
devtools::install_github("shashidhar22/LymphoSeq2", build_vignettes = TRUE)
```

## Getting started

To import AIRR-Seq data using `LymphoSeq2` we can use the
`readImmunoSeq` function. Currently the function can import data from
MiXCR, Adaptive ImmunoSEQ, BGI IR-SEQ, and 10X Genomic single cell VDJ
rearrangements.

``` r
l> library(LymphoSeq2)

Attaching package: ‘LymphoSeq2’

The following object is masked _by_ ‘.GlobalEnv’:

    merge_chains

Warning message:
replacing previous import ‘data.table:::=’ by ‘rlang:::=’ when loading ‘LymphoSeq2’ 
> study_files <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq2")
> study_table <- LymphoSeq2::readImmunoSeq(study_files)
Dataset Analysis:
  Files: 10, Total: 0.00 GB, Largest: 0.3 MB
  Available memory: 1.0 GB
Loading [============================================================] 10/10 files (100%) | ETA:  0s
```

To get a quick summary of repertoire characteristics, use the
`clonality` function. This will calculate many standard repertoire
diversity metrics such `clonality`, `gini coefficient`, `convergence`,
and `unique productive sequence` for each of the repertoires in the
input dataset.

``` r
> summary_table <- LymphoSeq2::clonality(study_table)
> summary_table
# A tibble: 10 × 10
   repertoire_id    total_sequences unique_productive_se…¹ total_count clonality
   <chr>                      <int>                  <int>       <int>     <dbl>
 1 TRB_CD4_949                 1000                    851       25769     0.440
 2 TRB_CD8_949                 1000                    811       26239     0.430
 3 TRB_CD8_CMV_369              414                    290        1794     0.328
 4 TRB_Unsorted_0              1000                    850       18161     0.280
 5 TRB_Unsorted_13…            1000                    848      178190     0.420
 6 TRB_Unsorted_14…            1000                    845       33669     0.387
 7 TRB_Unsorted_32              920                    783       31078     0.133
 8 TRB_Unsorted_369            1000                    835      339413     0.426
 9 TRB_Unsorted_83             1000                    835      236732     0.335
10 TRB_Unsorted_949            1000                    845        6549     0.303
# ℹ abbreviated name: ¹​unique_productive_sequences
# ℹ 5 more variables: gini_coefficient <dbl>, simpson_index <dbl>,
#   inverse_simpson <dbl>, top_productive_sequence <dbl>, convergence <dbl>
```

To compare samples with varying depth of sequencing, you can use the
`clonality` function to sample down all repertoires to a minimum number
of sequences. Since we randomly sample sequences from each repertoire,
in this mode the `clonality` function will repeat this operation for a
user specified number of `iterations` and caculate the average value for
all the diversity metrics.

``` r
> sampled_summary <- LymphoSeq2::clonality(study_table, rarefy = TRUE, iterations = 5, min_count = 1000)
> sampled_summary
# A tibble: 10 × 10
   repertoire_id    total_sequences unique_productive_se…¹ total_count clonality
   <chr>                      <dbl>                  <dbl>       <dbl>     <dbl>
 1 TRB_CD4_949                 163.                   137.        1000    0.303 
 2 TRB_CD8_949                 199                    162.        1000    0.292 
 3 TRB_CD8_CMV_369             272                    194.        1000    0.295 
 4 TRB_Unsorted_0              249                    210         1000    0.163 
 5 TRB_Unsorted_13…            178.                   145         1000    0.272 
 6 TRB_Unsorted_14…            207.                   173.        1000    0.255 
 7 TRB_Unsorted_32             413.                   356         1000    0.0851
 8 TRB_Unsorted_369            244.                   200.        1000    0.327 
 9 TRB_Unsorted_83             307.                   256         1000    0.270 
10 TRB_Unsorted_949            299.                   248.        1000    0.230 
# ℹ abbreviated name: ¹​unique_productive_sequences
# ℹ 5 more variables: gini_coefficient <dbl>, simpson_index <dbl>,
#   inverse_simpson <dbl>, top_productive_sequence <dbl>, convergence <dbl>
```
