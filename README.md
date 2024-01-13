
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
as Adaptive ImmunoSEQ and BGI-IRSeq, with support for 10X VDJ sequencing
coming soon. The platform also supports the importing of AIRR-seq data
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
library(LymphoSeq2)
#> Loading required package: data.table
study_files <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq2")
study_table <- LymphoSeq2::readImmunoSeq(study_files)
#> Registered S3 methods overwritten by 'readr':
#>   method                    from 
#>   as.data.frame.spec_tbl_df vroom
#>   as_tibble.spec_tbl_df     vroom
#>   format.col_spec           vroom
#>   print.col_spec            vroom
#>   print.collector           vroom
#>   print.date_names          vroom
#>   print.locale              vroom
#>   str.col_spec              vroom
```

To get a quick summary of repertoire characteristics, use the
`clonality` function. This will calculate many standard repertoire
diversity metrics such `clonality`, `gini coefficient`, `convergence`,
and `unique productive sequence` for each of the repertoires in the
input dataset.

``` r
summary_table <- LymphoSeq2::clonality(study_table)
summary_table
#> # A tibble: 10 × 8
#>    repertoire_id    total_sequences unique_productive_se…¹ total_count clonality
#>    <chr>                      <int>                  <int>       <dbl>     <dbl>
#>  1 TRB_CD4_949                 1000                    845       25769     0.443
#>  2 TRB_CD8_949                 1000                    794       26239     0.431
#>  3 TRB_CD8_CMV_369              414                    281        1794     0.332
#>  4 TRB_Unsorted_0              1000                    838       18161     0.281
#>  5 TRB_Unsorted_13…            1000                    838      178190     0.422
#>  6 TRB_Unsorted_14…            1000                    832       33669     0.389
#>  7 TRB_Unsorted_32              920                    767       31078     0.134
#>  8 TRB_Unsorted_369            1000                    830      339413     0.426
#>  9 TRB_Unsorted_83             1000                    823      236732     0.338
#> 10 TRB_Unsorted_949            1000                    831        6549     0.306
#> # ℹ abbreviated name: ¹​unique_productive_sequences
#> # ℹ 3 more variables: gini_coefficient <dbl>, top_productive_sequence <dbl>,
#> #   convergence <dbl>
```

To compare samples with varying depth of sequencing, you can use the
`clonality` function to sample down all repertoires to a minimum number
of sequences. Since we randomly sample sequences from each repertoire,
in this mode the `clonality` function will repeat this operation for a
user specified number of `iterations` and caculate the average value for
all the diversity metrics.

``` r
sampled_summary <- LymphoSeq2::clonality(study_table, rarefy = TRUE, iterations = 5, min_count = 1000)
sampled_summary
#> # A tibble: 10 × 8
#>    repertoire_id    total_sequences unique_productive_se…¹ total_count clonality
#>    <chr>                      <dbl>                  <dbl>       <dbl>     <dbl>
#>  1 TRB_CD4_949                 160.                   134         1000    0.311 
#>  2 TRB_CD8_949                 190.                   150.        1000    0.301 
#>  3 TRB_CD8_CMV_369             276                    188.        1000    0.298 
#>  4 TRB_Unsorted_0              248.                   209         1000    0.157 
#>  5 TRB_Unsorted_13…            184.                   149         1000    0.282 
#>  6 TRB_Unsorted_14…            217                    182.        1000    0.262 
#>  7 TRB_Unsorted_32             416.                   350.        1000    0.0928
#>  8 TRB_Unsorted_369            244.                   200.        1000    0.338 
#>  9 TRB_Unsorted_83             307.                   250         1000    0.272 
#> 10 TRB_Unsorted_949            307                    255.        1000    0.221 
#> # ℹ abbreviated name: ¹​unique_productive_sequences
#> # ℹ 3 more variables: gini_coefficient <dbl>, top_productive_sequence <dbl>,
#> #   convergence <dbl>
```
