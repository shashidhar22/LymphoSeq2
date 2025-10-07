#' Generate rarefaction curves for diversity estimation
#'
#' Estimate repertoire diversity across different sequencing depths using
#' rarefaction (interpolation) and extrapolation. This uses the iNEXT algorithm
#' to help determine if sequencing depth is sufficient to capture repertoire diversity.
#'
#' @param sample_table A tibble from [readImmunoSeq()] or [productiveSeq()] containing
#' "junction_aa", "duplicate_count", and "repertoire_id" columns. Can contain one or
#' multiple repertoires.
#' @param q Diversity order to calculate:
#'   * 0 (default): Species richness (number of unique clones)
#'   * 1: Shannon diversity (accounts for evenness)
#'   * 2: Simpson diversity (emphasizes abundant clones)
#' @param endpoint Maximum sequencing depth for extrapolation. Default is 100000.
#' Set higher to predict diversity at deeper sequencing.
#' @param nboot Number of bootstrap iterations for confidence intervals (default 10).
#' Higher values give more precise estimates but take longer.
#' @param conf Confidence level for intervals (default 0.95 for 95% CI)
#'
#' @return A tibble with rarefaction/extrapolation results:
#' * `m`: Sample size (sequencing depth)
#' * `Method`: "Rarefaction", "Observed", or "Extrapolation"
#' * `Order.q`: Diversity order (same as `q` parameter)
#' * `qD`: Estimated diversity at depth `m`
#' * `qD.LCL`: Lower confidence limit
#' * `qD.UCL`: Upper confidence limit
#' * `SC`: Standard error from bootstrap
#' * `repertoire_id`: Sample identifier
#'
#' @details
#' This function wraps the iNEXT package (Chao et al. 2014) for rarefaction and
#' extrapolation analysis. It converts the sample table into a matrix format where
#' rows are unique sequences and columns are repertoires, then runs iNEXT on all
#' samples simultaneously.
#'
#' Rarefaction vs Extrapolation:
#'
#' Rarefaction (m less than observed): Subsample sequences to depth m and count
#' unique clones. Shows how diversity increases with sequencing depth. Useful to
#' compare samples at equal depth.
#'
#' Extrapolation (m greater than observed): Predict diversity at deeper sequencing using
#' Chao1 estimator for unseen species. Shows whether sequencing is complete
#' (plateau) or more diversity remains (still increasing).
#'
#' How to interpret the curve:
#' Plateau reached = Sequencing depth is sufficient, most clones captured.
#' Still increasing steeply = Need deeper sequencing to capture full diversity.
#' Comparing samples = Use rarefied diversity at same depth, not raw counts.
#'
#' @references
#' Chao, A., et al. (2014). Rarefaction and extrapolation with Hill numbers:
#' a framework for sampling and estimation in species diversity studies.
#' Ecological Monographs, 84(1), 45-67.
#'
#' @examples
#' \dontrun{
#' file_path <- system.file("extdata", "TCRB_sequencing",
#'  package = "LymphoSeq2")
#' study_table <- LymphoSeq2::readImmunoSeq(path = file_path, threads = 1)
#' amino_table <- LymphoSeq2::productiveSeq(study_table,
#'   aggregate = "junction_aa",
#'   prevalence = TRUE
#' )
#' # Run on all samples at once
#' rarefaction_table <- LymphoSeq2::runINext(amino_table)
#' }
#' @export
runINext <- function(sample_table, q = 0, endpoint = 100000, nboot = 10, conf = 0.95) {
  # Check if iNEXT is available
  if (!requireNamespace("iNEXT", quietly = TRUE)) {
    stop("Package 'iNEXT' is required for rarefaction analysis. ",
         "Please install it with: install.packages('iNEXT')",
         call. = FALSE)
  }

  # Convert sample table to matrix format for iNEXT
  # Group by junction_aa and repertoire_id, sum duplicate counts
  rarefaction_matrix <- sample_table |>
    dplyr::group_by(junction_aa, repertoire_id) |>
    dplyr::summarise(duplicate_count = sum(duplicate_count), .groups = "drop") |>
    tidyr::pivot_wider(
      names_from = repertoire_id,
      id_cols = junction_aa,
      values_from = duplicate_count,
      values_fill = list(duplicate_count = 0)
    ) |>
    dplyr::ungroup() |>
    dplyr::select(-junction_aa) |>
    as.matrix()

  # Run iNEXT on the matrix
  rarefaction_result <- iNEXT::iNEXT(
    x = rarefaction_matrix,
    q = q,
    datatype = "abundance",
    endpoint = endpoint,
    se = TRUE,
    conf = conf,
    nboot = nboot
  )

  # Extract results and convert to tibble
  rarefaction_tables <- rarefaction_result$iNextEst$size_based |>
    dplyr::as_tibble() |>
    dplyr::rename(repertoire_id = Assemblage)

  return(rarefaction_tables)
}
