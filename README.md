
<!-- README.md is generated from README.Rmd. Please edit that file -->

# LymphoSeq2

<!-- badges: start -->

[![R-CMD-check](https://github.com/shashidhar22/LymphoSeq2/workflows/R-CMD-check/badge.svg)](https://github.com/shashidhar22/LymphoSeq2/actions)
<!-- badges: end -->

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
#> Registered S3 methods overwritten by 'ggalt':
#>   method                  from   
#>   grid.draw.absoluteGrob  ggplot2
#>   grobHeight.absoluteGrob ggplot2
#>   grobWidth.absoluteGrob  ggplot2
#>   grobX.absoluteGrob      ggplot2
#>   grobY.absoluteGrob      ggplot2
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
#>    repertoire_id     total_seq…¹ uniqu…² total…³ clona…⁴ gini_…⁵ top_p…⁶ conve…⁷
#>    <chr>                   <int>   <int>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
#>  1 TRB_CD4_949              1000     845   25769   0.443   0.867   30.1     1   
#>  2 TRB_CD8_949              1000     794   26239   0.431   0.903   19.3     1.01
#>  3 TRB_CD8_CMV_369           414     281    1794   0.332   0.761   16.5     1.12
#>  4 TRB_Unsorted_0           1000     838   18161   0.281   0.818    5.77    1   
#>  5 TRB_Unsorted_1320        1000     838  178190   0.422   0.902   14.6     1   
#>  6 TRB_Unsorted_1496        1000     832   33669   0.389   0.881   14.2     1   
#>  7 TRB_Unsorted_32           920     767   31078   0.134   0.601    4.87    1.01
#>  8 TRB_Unsorted_369         1000     830  339413   0.426   0.845   29.7     1   
#>  9 TRB_Unsorted_83          1000     823  236732   0.338   0.777   23.6     1   
#> 10 TRB_Unsorted_949         1000     831    6549   0.306   0.765   13.8     1   
#> # … with abbreviated variable names ¹​total_sequences,
#> #   ²​unique_productive_sequences, ³​total_count, ⁴​clonality, ⁵​gini_coefficient,
#> #   ⁶​top_productive_sequence, ⁷​convergence
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
#>    repertoire_id     total_seq…¹ uniqu…² total…³ clona…⁴ gini_…⁵ top_p…⁶ conve…⁷
#>    <chr>                   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>
#>  1 TRB_CD4_949              156.    129     1000  0.319    0.735   31.3     1   
#>  2 TRB_CD8_949              192.    154.    1000  0.293    0.734   18.6     1.01
#>  3 TRB_CD8_CMV_369          266     183.    1000  0.303    0.723   17.1     1.09
#>  4 TRB_Unsorted_0           254.    212.    1000  0.158    0.616    6.11    1.00
#>  5 TRB_Unsorted_1320        191.    156.    1000  0.277    0.726   14.6     1.02
#>  6 TRB_Unsorted_1496        206     168.    1000  0.260    0.710   14.0     1.00
#>  7 TRB_Unsorted_32          417.    349.    1000  0.0952   0.466    5.54    1.01
#>  8 TRB_Unsorted_369         245.    203     1000  0.336    0.711   30.2     1   
#>  9 TRB_Unsorted_83          320.    264.    1000  0.261    0.639   23.4     1   
#> 10 TRB_Unsorted_949         301.    249.    1000  0.222    0.630   13.0     1.00
#> # … with abbreviated variable names ¹​total_sequences,
#> #   ²​unique_productive_sequences, ³​total_count, ⁴​clonality, ⁵​gini_coefficient,
#> #   ⁶​top_productive_sequence, ⁷​convergence
```
