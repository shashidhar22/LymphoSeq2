Analysis of Adaptive Immune Receptor Repertoire Sequencing data with
LymphoSeq2
================
Shashidhar Ravishankar
2020-08-25

Adaptive Immune Receptor Repertoire Sequencing (AIRR-seq) provides a
unique opportunity to interrogate the adaptive immune repertoire under
various clinical conditions. The utility offered by this technology has
quickly garnered interest from a community of clinicians and researchers
investigating the immunological landscapes of a large spectrum of health
and disease states. LymphoSeq2 is a toolikt that allows users to import,
manipulate and visualize AIRR-Seq dta from various AIRR-Seq assays such
as Adaptive ImmunoSEQ and BGI-IRSeq, with support for 10X VDJ sequencing
coming soon. The platform also supports the importing of AIRR-seq data
processed using the MiXCR pipeline. The vignette highlights some of the
key features of LymphoSeq2.

# Importing data

The function `readImmunoSeq` imports AIRR-seq receptor files from
Adaptive ImmunoSEQ assay as well well as BGI-IRSeq assay. The sequences
can be (.tsv) files processed using one of the three following
platforms: Adaptive Biotechnologies ImmunoSEQ analyzer, BGI IR-SEQ
iMonitor platform, and the MiXCR pipeline for AIRR-seq data analysis.
The function has the ability to identify file type based on the headers
provided in the (.tsv) file, accordingly the data is tranformed into a
format that is compatible AIRR-Community guidelines
(<https://github.com/airr-community/airr-standards>).

To explore the features of LymphoSeq, this package includes 2 example
data sets. The first is a data set of T cell receptor beta (TCRB)
sequencing from 10 blood samples acquired serially from a single patient
who underwent a bone marrow transplant (Kanakry, C.G., et al. JCI
Insight 2016;1(5):pii: e86252). The second, is a data set of B cell
receptor immunoglobulin heavy (IGH) chain sequencing from Burkitt
lymphoma tumor biopsies acquired from 10 different individuals
(Lombardo, K.A., et al. Blood Advances 2017 1:535-544). To improve
performance, both data sets contain only the top 1,000 most frequent
sequences. The complete data sets are publicly available through
Adapatives’ immuneACCESS portal. As shown in the example below, you can
specify the path to the example data sets using the command

``` r
system.file("extdata", "TCRB_sequencing", package = "LymphoSeq") #For the TCRB files
#> [1] ""
system.file("extdata", "IGH_sequencing", package = "LymphoSeq") #For the IGH files.
#> [1] ""
```

`readImmunoSeq` can take as input, a single file name, a list of files
or the path to a directory containing AIRR-seq data. The columns are
renamed to follow AIRR-community guidelines based on the input file
type. The function returns a tibble with individual file names set as
the `repertoire_id`. The CDR3 nucleotide and amino acid sequences are
denoted by the `junction` and `junction_aa` fields respectively. The
counts of the CDR3 sequeneces observed, and their frequency in each
individual repertoire is denoted by the `duplicate_count` and
`duplicate_frequency` field respectively.

``` r
library(LymphoSeq2)
library(tidyverse)
#> ── Attaching packages ──────────────────────────────────────────────── tidyverse 1.3.0 ──
#> ✓ ggplot2 3.3.2     ✓ purrr   0.3.4
#> ✓ tibble  3.0.3     ✓ dplyr   1.0.2
#> ✓ tidyr   1.1.1     ✓ stringr 1.4.0
#> ✓ readr   1.3.1     ✓ forcats 0.5.0
#> ── Conflicts ─────────────────────────────────────────────────── tidyverse_conflicts() ──
#> x dplyr::filter() masks stats::filter()
#> x dplyr::lag()    masks stats::lag()
study_files <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq2") 
study_table <- readImmunoSeq(study_files)
```

Looking at the `study_table` we see a tibble with 12 columns and 9334
rows

``` r
study_table
#> # A tibble: 9,334 x 12
#>    repertoire_id junction junction_aa v_call j_call d_call v_family j_family
#>    <chr>         <chr>    <chr>       <chr>  <chr>  <chr>  <chr>    <chr>   
#>  1 TRB_CD4_949   AAATCTT… <NA>        unres… TCRBJ… unres… unresol… TCRBJ01 
#>  2 TRB_CD4_949   AACATGA… CASSLGGGYG… TCRBV… TCRBJ… TCRBD… TCRBV13  TCRBJ01 
#>  3 TRB_CD4_949   AACCCGA… CASSINPVQG… TCRBV… TCRBJ… TCRBD… TCRBV19  TCRBJ02 
#>  4 TRB_CD4_949   AACCCGA… CASSIPTGLA… TCRBV… TCRBJ… TCRBD… TCRBV19  TCRBJ02 
#>  5 TRB_CD4_949   AACCTGA… CASSVEVGSA… TCRBV… TCRBJ… unres… TCRBV09  TCRBJ01 
#>  6 TRB_CD4_949   AACCTGA… CASSVEVAGW… TCRBV… TCRBJ… TCRBD… TCRBV09  TCRBJ01 
#>  7 TRB_CD4_949   AACCTGA… CASSVMVGTE… TCRBV… TCRBJ… unres… TCRBV09  TCRBJ01 
#>  8 TRB_CD4_949   AACGCCT… CASSPGGTGM… TCRBV… TCRBJ… TCRBD… TCRBV05  TCRBJ02 
#>  9 TRB_CD4_949   AACGCCT… CASSDSGVPG… TCRBV… TCRBJ… TCRBD… TCRBV05  TCRBJ01 
#> 10 TRB_CD4_949   AACGCCT… CASSFLRGGS… TCRBV… TCRBJ… TCRBD… TCRBV05  TCRBJ01 
#> # … with 9,324 more rows, and 4 more variables: d_family <chr>,
#> #   reading_frame <chr>, duplicate_count <int>, duplicate_frequency <dbl>
```

Since the study table is a tibble, we can use tidyverse syntax to
extract a list of sample names

``` r
study_table %>% pull(repertoire_id) %>% unique()
#>  [1] "TRB_CD4_949"       "TRB_CD8_949"       "TRB_CD8_CMV_369"  
#>  [4] "TRB_Unsorted_0"    "TRB_Unsorted_1320" "TRB_Unsorted_1496"
#>  [7] "TRB_Unsorted_32"   "TRB_Unsorted_369"  "TRB_Unsorted_83"  
#> [10] "TRB_Unsorted_949"
```

# Subsetting Data

The tibble structure of the TCR data allows for easy subsampling of
data. To select the TCR sequences from any given samples in the dataset,
the `filter` function from the `dplyr` package can be used.

``` r
TRB_Unsorted_0 <- study_table %>% filter(repertoire_id == "TRB_Unsorted_0")
TRB_Unsorted_0
#> # A tibble: 1,000 x 12
#>    repertoire_id junction junction_aa v_call j_call d_call v_family j_family
#>    <chr>         <chr>    <chr>       <chr>  <chr>  <chr>  <chr>    <chr>   
#>  1 TRB_Unsorted… AAAAGAA… <NA>        TCRBV… TCRBJ… TCRBD… TCRBV19  TCRBJ02 
#>  2 TRB_Unsorted… AAACCTT… <NA>        TCRBV… TCRBJ… TCRBD… TCRBV04  TCRBJ01 
#>  3 TRB_Unsorted… AACATGA… CSVRMLKTGV… TCRBV… TCRBJ… TCRBD… TCRBV29  TCRBJ01 
#>  4 TRB_Unsorted… AACATGA… CSVEERTVSG… TCRBV… TCRBJ… TCRBD… TCRBV29  TCRBJ02 
#>  5 TRB_Unsorted… AACGCCT… CASSTGTGVY… TCRBV… TCRBJ… TCRBD… TCRBV05  TCRBJ01 
#>  6 TRB_Unsorted… AACGCCT… CASSYNSGTS… TCRBV… TCRBJ… TCRBD… TCRBV05  TCRBJ02 
#>  7 TRB_Unsorted… AACGCCT… CASSSDRGVD… TCRBV… TCRBJ… unres… TCRBV05  TCRBJ01 
#>  8 TRB_Unsorted… AACGCCT… CASSLLGLTN… TCRBV… TCRBJ… TCRBD… TCRBV05  TCRBJ01 
#>  9 TRB_Unsorted… AACGCCT… CASSLDQGRN… TCRBV… TCRBJ… unres… TCRBV05  TCRBJ01 
#> 10 TRB_Unsorted… AACGCCT… CASSLGSSGA… TCRBV… TCRBJ… TCRBD… TCRBV05  TCRBJ02 
#> # … with 990 more rows, and 4 more variables: d_family <chr>,
#> #   reading_frame <chr>, duplicate_count <int>, duplicate_frequency <dbl>
```

The `str_detect` function from `stringr` package can be used in
conjunction with the `filter` to find samples using a pattern

``` r
CMV <- study_table %>% filter(str_detect(repertoire_id, "CMV"))
CMV
#> # A tibble: 414 x 12
#>    repertoire_id junction junction_aa v_call j_call d_call v_family j_family
#>    <chr>         <chr>    <chr>       <chr>  <chr>  <chr>  <chr>    <chr>   
#>  1 TRB_CD8_CMV_… AAAAGAA… <NA>        TCRBV… TCRBJ… TCRBD… TCRBV19  TCRBJ02 
#>  2 TRB_CD8_CMV_… AACCTGA… CASSAGRTMP… TCRBV… TCRBJ… TCRBD… TCRBV09  TCRBJ01 
#>  3 TRB_CD8_CMV_… AACTGGA… <NA>        TCRBV… TCRBJ… TCRBD… TCRBV14  TCRBJ02 
#>  4 TRB_CD8_CMV_… AAGAACC… CASSIDVVRT… TCRBV… TCRBJ… TCRBD… TCRBV19  TCRBJ01 
#>  5 TRB_CD8_CMV_… AAGAAGC… CAWKAPRSTN… TCRBV… TCRBJ… unres… TCRBV30  TCRBJ01 
#>  6 TRB_CD8_CMV_… AAGAAGC… CAWSGGQGAS… TCRBV… TCRBJ… TCRBD… TCRBV30  TCRBJ02 
#>  7 TRB_CD8_CMV_… AAGAAGC… CAWSLKGAMN… TCRBV… TCRBJ… TCRBD… TCRBV30  TCRBJ01 
#>  8 TRB_CD8_CMV_… AAGATCC… CASSDSFNQP… unres… TCRBJ… TCRBD… unresol… TCRBJ01 
#>  9 TRB_CD8_CMV_… AAGATCC… CASSFGGLTE… TCRBV… TCRBJ… TCRBD… TCRBV07  TCRBJ01 
#> 10 TRB_CD8_CMV_… AAGATCC… CATAGLEGDE… TCRBV… TCRBJ… unres… TCRBV02  TCRBJ02 
#> # … with 404 more rows, and 4 more variables: d_family <chr>,
#> #   reading_frame <chr>, duplicate_count <int>, duplicate_frequency <dbl>
```

A metadata file for the TCR sequencing samples can easily be combined
with the `study_table` by reading in the metadata file as a tibble and
using the `dplyr::left_join` function to merge the two tables. In the
example below, a metadata file is imported for the example TCRB data set
which contains information on the number of days post bone marrow
transplant the sample was collected and the cellular phenotype the blood
sample was sorted for prior to sequencing.

``` r
TCRB_metadata <- read.csv(system.file("extdata", "TCRB_metadata.csv", package = "LymphoSeq2"))
TCRB_metadata
#>              samples  day phenotype
#> 1     TRB_Unsorted_0    0  Unsorted
#> 2    TRB_Unsorted_32   32  Unsorted
#> 3    TRB_Unsorted_83   82  Unsorted
#> 4    TRB_CD8_CMV_369  369  CD8+CMV+
#> 5   TRB_Unsorted_369  369  Unsorted
#> 6        TRB_CD4_949  949      CD4+
#> 7        TRB_CD8_949  949      CD8+
#> 8   TRB_Unsorted_949  949  Unsorted
#> 9  TRB_Unsorted_1320 1320  Unsorted
#> 10 TRB_Unsorted_1496 1496  Unsorted
```

``` r
study_table <- left_join(study_table, TCRB_metadata, by = c("repertoire_id" = "samples"))
study_table
#> # A tibble: 9,334 x 14
#>    repertoire_id junction junction_aa v_call j_call d_call v_family j_family
#>    <chr>         <chr>    <chr>       <chr>  <chr>  <chr>  <chr>    <chr>   
#>  1 TRB_CD4_949   AAATCTT… <NA>        unres… TCRBJ… unres… unresol… TCRBJ01 
#>  2 TRB_CD4_949   AACATGA… CASSLGGGYG… TCRBV… TCRBJ… TCRBD… TCRBV13  TCRBJ01 
#>  3 TRB_CD4_949   AACCCGA… CASSINPVQG… TCRBV… TCRBJ… TCRBD… TCRBV19  TCRBJ02 
#>  4 TRB_CD4_949   AACCCGA… CASSIPTGLA… TCRBV… TCRBJ… TCRBD… TCRBV19  TCRBJ02 
#>  5 TRB_CD4_949   AACCTGA… CASSVEVGSA… TCRBV… TCRBJ… unres… TCRBV09  TCRBJ01 
#>  6 TRB_CD4_949   AACCTGA… CASSVEVAGW… TCRBV… TCRBJ… TCRBD… TCRBV09  TCRBJ01 
#>  7 TRB_CD4_949   AACCTGA… CASSVMVGTE… TCRBV… TCRBJ… unres… TCRBV09  TCRBJ01 
#>  8 TRB_CD4_949   AACGCCT… CASSPGGTGM… TCRBV… TCRBJ… TCRBD… TCRBV05  TCRBJ02 
#>  9 TRB_CD4_949   AACGCCT… CASSDSGVPG… TCRBV… TCRBJ… TCRBD… TCRBV05  TCRBJ01 
#> 10 TRB_CD4_949   AACGCCT… CASSFLRGGS… TCRBV… TCRBJ… TCRBD… TCRBV05  TCRBJ01 
#> # … with 9,324 more rows, and 6 more variables: d_family <chr>,
#> #   reading_frame <chr>, duplicate_count <int>, duplicate_frequency <dbl>,
#> #   day <int>, phenotype <chr>
```

Now the metadata information can be used to further subset the data. For
instance to select all “Unsorted” samples collected more than 300 days
after bone marrow transplant, we would use the following code

``` r
unsorted_300 <- study_table %>% filter(day > 300 & phenotype == "Unsorted")
unsorted_300
#> # A tibble: 4,000 x 14
#>    repertoire_id junction junction_aa v_call j_call d_call v_family j_family
#>    <chr>         <chr>    <chr>       <chr>  <chr>  <chr>  <chr>    <chr>   
#>  1 TRB_Unsorted… AACATGA… CSAYGTATAT… TCRBV… TCRBJ… TCRBD… TCRBV29  TCRBJ01 
#>  2 TRB_Unsorted… AACATGA… CASSSGRTYE… TCRBV… TCRBJ… TCRBD… TCRBV13  TCRBJ02 
#>  3 TRB_Unsorted… AACCTGA… CASTSRGGGE… TCRBV… TCRBJ… TCRBD… TCRBV09  TCRBJ02 
#>  4 TRB_Unsorted… AACCTGA… CASSAGRTMP… TCRBV… TCRBJ… TCRBD… TCRBV09  TCRBJ01 
#>  5 TRB_Unsorted… AACCTGA… CASSVEVGSA… TCRBV… TCRBJ… unres… TCRBV09  TCRBJ01 
#>  6 TRB_Unsorted… AACCTGA… CASSVAGGQP… TCRBV… TCRBJ… TCRBD… TCRBV09  TCRBJ01 
#>  7 TRB_Unsorted… AACCTGA… CASSVMVGTE… TCRBV… TCRBJ… unres… TCRBV09  TCRBJ01 
#>  8 TRB_Unsorted… AACCTGA… CASSVLQGTE… TCRBV… TCRBJ… TCRBD… TCRBV09  TCRBJ01 
#>  9 TRB_Unsorted… AACCTGA… CASSVGMNSP… TCRBV… TCRBJ… TCRBD… TCRBV09  TCRBJ01 
#> 10 TRB_Unsorted… AACGCCT… CASSQEGGRD… TCRBV… TCRBJ… TCRBD… TCRBV05  TCRBJ02 
#> # … with 3,990 more rows, and 6 more variables: d_family <chr>,
#> #   reading_frame <chr>, duplicate_count <int>, duplicate_frequency <dbl>,
#> #   day <int>, phenotype <chr>
```

# Extracting productive sequences

When AIRR-seq samples are derived from genomic DNA rather than
complimentary DNA made from RNA, then you will find productive and
unproductive sequences. Productive sequences are defined as in-frame
sequences without any early stop codons. To filter out these productive
sequences, you can use the `productiveSeq` to remove unproductive
sequences and recompute the duplicate\_frequency to reflect the
productive amino acid or nucleotide sequence frequencies.

If you are interested in just the complementarity determining region 3
(CDR3) amino acid sequences, then set aggregate to `junction_aa` and the
`duplicate_count` for duplicate amino acid sequences will be summed. The
resulting tibble will have `junction_aa`, `duplicate_count`,
`duplicate_frequency`, `reading_frame`, and the most frequent VDJ gene
combinations for each of the duplicated amino acid sequences and the
corresponding gene family names. These gene names are only kept for
consistency of the tibble structure, but since a single amino acid
sequence can be generated from different VDJ combinations, it is
inadvisable to use these values for downstream analysis

``` r
aa_table <- productiveSeq(study_table = study_table, aggregate = "junction_aa", 
                          prevalence = FALSE)
aa_table
#> # A tibble: 7,533 x 11
#>    repertoire_id junction_aa v_call d_call j_call v_family d_family j_family
#>    <chr>         <chr>       <chr>  <chr>  <chr>  <chr>    <chr>    <chr>   
#>  1 TRB_CD4_949   CAIKPGQGAS… TCRBV… TCRBD… TCRBJ… TCRBV10  TCRBD01  TCRBJ01 
#>  2 TRB_CD4_949   CAIRAGTSTD… TCRBV… TCRBD… TCRBJ… TCRBV10  TCRBD02  TCRBJ02 
#>  3 TRB_CD4_949   CAISDETPGE… TCRBV… unres… TCRBJ… TCRBV10  unresol… TCRBJ02 
#>  4 TRB_CD4_949   CAISDLGRGD… TCRBV… TCRBD… TCRBJ… TCRBV10  TCRBD01  TCRBJ01 
#>  5 TRB_CD4_949   CAISDLKEQP… TCRBV… TCRBD… TCRBJ… TCRBV10  TCRBD02  TCRBJ01 
#>  6 TRB_CD4_949   CAISDQGGDQ… TCRBV… TCRBD… TCRBJ… TCRBV10  TCRBD02  TCRBJ01 
#>  7 TRB_CD4_949   CAISEREQGA… TCRBV… TCRBD… TCRBJ… TCRBV10  TCRBD01  TCRBJ01 
#>  8 TRB_CD4_949   CAISEWSGSS… TCRBV… unres… TCRBJ… TCRBV10  unresol… TCRBJ02 
#>  9 TRB_CD4_949   CAISGQGSTE… TCRBV… TCRBD… TCRBJ… TCRBV10  TCRBD01  TCRBJ01 
#> 10 TRB_CD4_949   CAISLNSGGA… TCRBV… TCRBD… TCRBJ… TCRBV10  TCRBD02  TCRBJ02 
#> # … with 7,523 more rows, and 3 more variables: reading_frame <chr>,
#> #   duplicate_count <int>, duplicate_frequency <dbl>
```

Alternatively you can aggregate by `junction` to group sequences by CDR3
nucleotide sequences. This option will produce a tibble similar to the
output of `readImmunoSeq`. Many of the functions within LymphoSeq2 use
the results from `productiveSeq` function. Please be sure to check the
function documentation.

``` r
nuc_table <- productiveSeq(study_table = study_table, aggregate = "junction", 
                          prevalence = FALSE)
nuc_table
#> # A tibble: 7,679 x 12
#>    repertoire_id junction junction_aa v_call d_call j_call v_family d_family
#>    <chr>         <chr>    <chr>       <chr>  <chr>  <chr>  <chr>    <chr>   
#>  1 TRB_CD4_949   AACATGA… CASSLGGGYG… TCRBV… TCRBD… TCRBJ… TCRBV13  TCRBD02 
#>  2 TRB_CD4_949   AACCCGA… CASSINPVQG… TCRBV… TCRBD… TCRBJ… TCRBV19  TCRBD01 
#>  3 TRB_CD4_949   AACCCGA… CASSIPTGLA… TCRBV… TCRBD… TCRBJ… TCRBV19  TCRBD02 
#>  4 TRB_CD4_949   AACCTGA… CASSVEVGSA… TCRBV… unres… TCRBJ… TCRBV09  unresol…
#>  5 TRB_CD4_949   AACCTGA… CASSVEVAGW… TCRBV… TCRBD… TCRBJ… TCRBV09  TCRBD01 
#>  6 TRB_CD4_949   AACCTGA… CASSVMVGTE… TCRBV… unres… TCRBJ… TCRBV09  unresol…
#>  7 TRB_CD4_949   AACGCCT… CASSPGGTGM… TCRBV… TCRBD… TCRBJ… TCRBV05  TCRBD01 
#>  8 TRB_CD4_949   AACGCCT… CASSDSGVPG… TCRBV… TCRBD… TCRBJ… TCRBV05  TCRBD01 
#>  9 TRB_CD4_949   AACGCCT… CASSFLRGGS… TCRBV… TCRBD… TCRBJ… TCRBV05  TCRBD02 
#> 10 TRB_CD4_949   AACGCCT… CASSLATKNT… TCRBV… TCRBD… TCRBJ… TCRBV05  TCRBD01 
#> # … with 7,669 more rows, and 4 more variables: j_family <chr>,
#> #   reading_frame <chr>, duplicate_count <int>, duplicate_frequency <dbl>
```

If the parameter prevalence is set to TRUE, then a new column is added
to each of the data frames giving the prevalence (%) of each TCR beta
CDR3 amino acid sequence in 55 healthy donor peripheral blood samples.
Values range from 0 to 100% where 100% means the sequence appeared in
the blood of all 55 individuals.

Notice in the example below that there are no amino acid sequences given
in the first and fourth row of the study\_table table for sample
“TRB\_Unsorted\_949”. This is because the nucleotide sequence is out
of frame and does not produce a productively transcribed amino acid
sequence. If an asterisk (\*) appears in the amino acid sequences, this
would indicate an early stop codon.

``` r
study_table %>% filter(repertoire_id == "TRB_Unsorted_0") 
#> # A tibble: 1,000 x 14
#>    repertoire_id junction junction_aa v_call j_call d_call v_family j_family
#>    <chr>         <chr>    <chr>       <chr>  <chr>  <chr>  <chr>    <chr>   
#>  1 TRB_Unsorted… AAAAGAA… <NA>        TCRBV… TCRBJ… TCRBD… TCRBV19  TCRBJ02 
#>  2 TRB_Unsorted… AAACCTT… <NA>        TCRBV… TCRBJ… TCRBD… TCRBV04  TCRBJ01 
#>  3 TRB_Unsorted… AACATGA… CSVRMLKTGV… TCRBV… TCRBJ… TCRBD… TCRBV29  TCRBJ01 
#>  4 TRB_Unsorted… AACATGA… CSVEERTVSG… TCRBV… TCRBJ… TCRBD… TCRBV29  TCRBJ02 
#>  5 TRB_Unsorted… AACGCCT… CASSTGTGVY… TCRBV… TCRBJ… TCRBD… TCRBV05  TCRBJ01 
#>  6 TRB_Unsorted… AACGCCT… CASSYNSGTS… TCRBV… TCRBJ… TCRBD… TCRBV05  TCRBJ02 
#>  7 TRB_Unsorted… AACGCCT… CASSSDRGVD… TCRBV… TCRBJ… unres… TCRBV05  TCRBJ01 
#>  8 TRB_Unsorted… AACGCCT… CASSLLGLTN… TCRBV… TCRBJ… TCRBD… TCRBV05  TCRBJ01 
#>  9 TRB_Unsorted… AACGCCT… CASSLDQGRN… TCRBV… TCRBJ… unres… TCRBV05  TCRBJ01 
#> 10 TRB_Unsorted… AACGCCT… CASSLGSSGA… TCRBV… TCRBJ… TCRBD… TCRBV05  TCRBJ02 
#> # … with 990 more rows, and 6 more variables: d_family <chr>,
#> #   reading_frame <chr>, duplicate_count <int>, duplicate_frequency <dbl>,
#> #   day <int>, phenotype <chr>
```

After `productiveSeq` is run, the unproductive sequences are removed and
the `duplicate_frequency` is recalculated for each sequence. If there
were two identical amino acid sequences that differed in their
nucleotide sequence, they would be combined and their counts added
together.

``` r
aa_table %>% filter(repertoire_id == "TRB_Unsorted_0")
#> # A tibble: 833 x 11
#>    repertoire_id junction_aa v_call d_call j_call v_family d_family j_family
#>    <chr>         <chr>       <chr>  <chr>  <chr>  <chr>    <chr>    <chr>   
#>  1 TRB_Unsorted… CAAGDTTLYE… TCRBV… TCRBD… TCRBJ… TCRBV06  TCRBD02  TCRBJ02 
#>  2 TRB_Unsorted… CAARGGGESY… TCRBV… TCRBD… TCRBJ… TCRBV06  TCRBD02  TCRBJ02 
#>  3 TRB_Unsorted… CAGGRGLEEA… TCRBV… TCRBD… TCRBJ… TCRBV28  TCRBD01  TCRBJ01 
#>  4 TRB_Unsorted… CAGQGVGYTE… TCRBV… unres… TCRBJ… TCRBV02  unresol… TCRBJ01 
#>  5 TRB_Unsorted… CAGRVGYTF   TCRBV… unres… TCRBJ… TCRBV02  unresol… TCRBJ01 
#>  6 TRB_Unsorted… CAIAQPNEKL… TCRBV… unres… TCRBJ… TCRBV10  unresol… TCRBJ01 
#>  7 TRB_Unsorted… CAIGAAYEQYF TCRBV… TCRBD… TCRBJ… TCRBV10  TCRBD02  TCRBJ02 
#>  8 TRB_Unsorted… CAIGYNQPQHF TCRBV… unres… TCRBJ… TCRBV10  unresol… TCRBJ01 
#>  9 TRB_Unsorted… CAINPGGTGN… TCRBV… TCRBD… TCRBJ… TCRBV10  TCRBD01  TCRBJ01 
#> 10 TRB_Unsorted… CAIRDGYNSP… TCRBV… unres… TCRBJ… TCRBV10  unresol… TCRBJ01 
#> # … with 823 more rows, and 3 more variables: reading_frame <chr>,
#> #   duplicate_count <int>, duplicate_frequency <dbl>
```

# Create a table of summary statistics

To create a table summarizing the total number of sequences, number of
unique productive sequences, number of genomes, clonality, Gini
coefficient, and the frequency (%) of the top productive sequence,
Simpson index, Inverse Simpson index, Hill diversity index, Chao1 index
and Kemp index in each imported file, use the function `clonality`.

``` r
clonality(study_table = study_table)
#> # A tibble: 10 x 12
#>    repertoire_id total_sequences unique_producti… total_count clonality
#>    <chr>                   <int>            <int>       <int>     <dbl>
#>  1 TRB_CD4_949              1000              845       25769     0.443
#>  2 TRB_CD8_949              1000              794       26239     0.431
#>  3 TRB_CD8_CMV_…             414              281        1794     0.332
#>  4 TRB_Unsorted…            1000              838       18161     0.281
#>  5 TRB_Unsorted…            1000              838      178190     0.422
#>  6 TRB_Unsorted…            1000              832       33669     0.389
#>  7 TRB_Unsorted…             920              767       31078     0.134
#>  8 TRB_Unsorted…            1000              830      339413     0.426
#>  9 TRB_Unsorted…            1000              823      236732     0.338
#> 10 TRB_Unsorted…            1000              831        6549     0.306
#> # … with 7 more variables: simpson_index <dbl>, inverse_simpson <dbl>,
#> #   gini_coefficient <dbl>, top_productive_sequence <dbl>, chao_estimate <dbl>,
#> #   kemp_estimate <dbl>, hill_estimate <dbl>
```

The clonality score is derived from the Shannon entropy, which is
calculated from the frequencies of all productive sequences divided by
the logarithm of the total number of unique productive sequences. This
normalized entropy value is then inverted (1 - normalized entropy) to
produce the clonality metric.

The Gini coefficient, Chao1 estimate, Kemp estimate, Hill estimate,
Simpson index and Inverse Simpson index are alternative metric to
measure sequence diversity within the immune repertoire.

The Gini coefficient is an alternative metric used to calculate
repertoire diversity and is derived from the Lorenz curve. The Lorenz
curve is drawn such that x-axis represents the cumulative percentage of
unique sequences and the y-axis represents the cumulative percentage of
reads. A line passing through the origin with a slope of 1 reflects
equal frequencies of all clones. The Gini coefficient is the ratio of
the area between the line of equality and the observed Lorenz curve over
the total area under the line of equality.

# Calculate clonal relatedness

One of the drawbacks of the clonality metric is that it does not take
into account sequence similarity. This is particularly important when
studying affinity maturation or B cell malignancies(Lombardo, K.A., et
al. Blood Advances 2017 1:535-544). Clonal relatedness is a useful
metric that takes into account sequence similarity without regard for
clonal frequency. It is defined as the proportion of nucleotide
sequences that are related by a defined edit distance threshold. The
value ranges from 0 to 1 where 0 indicates no sequences are related and
1 indicates all sequences are related. Edit distance is a way of
quantifying how dissimilar two sequences are to one another by counting
the minimum number of operations required to transform one sequence into
the other. For example, an edit distance of 0 means the sequences are
identical and an edit distance of 1 indicates that the sequences
different by a single amino acid or nucleotide.

``` r
IGH_path <- system.file("extdata", "IGH_sequencing", package = "LymphoSeq2")
IGH_table <- readImmunoSeq(path = IGH_path)
clonalRelatedness(study_table = IGH_table, editDistance = 10)
#> # A tibble: 10 x 2
#>    repertoire_id     clonalRelatedness
#>    <chr>                         <dbl>
#>  1 IGH_MVQ108911A_BL           0.631  
#>  2 IGH_MVQ194745A_BL           0.847  
#>  3 IGH_MVQ81231A_BL            0.631  
#>  4 IGH_MVQ89037A_BL            0.290  
#>  5 IGH_MVQ90143A_BL            0.00610
#>  6 IGH_MVQ92552A_BL            0.273  
#>  7 IGH_MVQ93505A_BL            0.408  
#>  8 IGH_MVQ93631A_BL            0.757  
#>  9 IGH_MVQ94865A_BL            0.007  
#> 10 IGH_MVQ95413A_BL            0.00364
```

# Draw a phylogenetic tree

A phylogenetic tree is a useful way to visualize the similarity between
sequences. The `phyloTree` function create a phylogenetic tree of a
single sample using neighbor joining tree estimation for amino acid or
nucleotide CDR3 sequences. Each leaf in the tree represents a sequence
color coded by the V, D, and J gene usage. The number next to each leaf
refers to the sequence count. A triangle shaped leaf indicates the most
frequent sequence. The distance between leaves on the horizontal axis
corresponds to the sequence similarity (i.e. the further apart the
leaves are horizontally, the less similar the sequences are to one
another).

``` r
nuc_IGH_table <- productiveSeq(study_table = IGH_table, aggregate = "junction")
#phyloTree(study_table, nuc_IGH_table, repertoire_id = "IGH_MVQ92552A_BL", type = "junction", 
#         layout = "rectangular")
```
