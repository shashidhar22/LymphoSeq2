#' Antigen DB
#'
#' @description
#' The database is created from flat file versions of IEDB,
#' McPAS-TCR, and VdjDB slim in that order and generates an RDA file of all
#' T-cells and B-cells recorded to show antigenic specificity in these three
#' databases. This is a script for internal use and will be run once every six
#' months to update the databases.
#' @format ## `antigen_db`
#' A tibble with 320,786 rows and 16 columns:
#' \describe{
#'   \item{tra_cdr3_aa}{T-cell receptor alpha chain amino acid sequence}
#'   \item{gene}{Antigen gene/protein name}
#'   \item{epitope}{Antigen epitope sequence/target}
#'   \item{pathology}{Pathology associated with epitope}
#'   \item{antigen}{Name of antigen}
#'   \item{tra_v_call}{T-cell receptor alpha chain V gene}
#'   \item{tra_j_call}{T-cell receptor alpha chain J gene}
#'   \item{mhc_allele}{MHC allele associated with the epitope}
#'   \item{reference}{Reference for known antigenic specificity}
#'   \item{score}{Confidence of the antigenic specificity}
#'   \item{cell_type}{Immune cell type}
#'   \item{source}{Database name}
#'   \item{trb_cdr3_aa}{T-cell receptor beta chain amino acid sequence}
#'   \item{trb_v_call}{T-cell receptor beta chain V gene}
#'   \item{trb_j_call}{T-cell receptor beta chain J gene}
#'   \item{Species}{Species}
#' }
#' @source
#' The flat files for this function can be
#' downloaded at the following links:
#' 1. IEDB: https://www.iedb.org/downloader.php?file_name=doc/receptor_full_v3.zip
#' 2. McPAS-TCR: http://friedmanlab.weizmann.ac.il/McPAS-TCR/
#' 3. VdjDB: https://github.com/antigenomics/vdjdb-db/releases
#' Note: Here the vdjdb.slim.txt is used as the input
#' Note: The order of the input matters, please provide the flat file paths in
#' this order IEDB, McPAS-TCR, VdjDB.
"antigen_db"
#' Prevalence TRB
#'
#' @description
#' The database describes the frequency at which a CDR3 amino acid sequence was
#' found in cohort of 55 PBMC samples from healthy individuals
#' @format ## `prevalenceTRB`
#' A tibble with 11,724,294 rows and 2 columns:
#' \describe{
#'   \item{prevalence}{Frequency of sequences in 55 healthy PBMC samples}
#'   \item{aminoAcid}{T-cell receptor beta chain amino acid sequence}
#' }
#' @source
#' TCR beta sequencing data from PBMCs of 55 healthy individuals
"prevalenceTRB"
#' Published TRB
#'
#' @description
#' The database describes T-cell beta chain amino acid sequences found to be
#' associated with antigenic specificity
#' @format ## `publishedTRB`
#' A tibble with 11,724,294 rows and 2 columns:
#' \describe{
#'   \item{prevalence}{Frequency of sequences in 55 healthy PBMC samples}
#'   \item{aminoAcid}{T-cell receptor beta chain amino acid sequence}
#'   \item{PMID}{Pubmed ID}
#'   \item{HLA}{Known MHC restriction}
#'   \item{antigen}{Antigen name}
#'   \item{epitope}{Known epitope sequence}
#' }
#' @source
#' T-cells previously published T-cell beta chain CDR3 amino acid sequences
#' with known antigenic specificity
"publishedTRB"
