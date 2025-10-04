#' Calculate repertoire diversity metrics
#'
#' Quantify immune repertoire diversity using multiple complementary metrics
#' including Shannon clonality, Gini coefficient, Simpson indices, and TCR/BCR
#' convergence. Optionally correct for sequencing depth bias using rarefaction.
#'
#' @param study_table A tibble of antigen receptor sequences from [readImmunoSeq()].
#' Must contain "junction_aa", "duplicate_count", and "duplicate_frequency" columns.
#' Use productive junction sequences (not aggregated by amino acid) for accurate
#' clonality estimates.
#' @param rarefy Logical. Should diversity be normalized for sequencing depth?
#'    * `TRUE`: Apply rarefaction by subsampling all repertoires to `min_count`
#'    depth, repeating for `iterations`, and averaging the results. Use this when
#'    comparing samples with different sequencing depths.
#'    * `FALSE` (default): Calculate raw diversity metrics without normalization.
#' @param iterations Number of bootstrap iterations for rarefaction (default 100).
#' Higher values increase precision but take longer to compute.
#' @param min_count Target sequencing depth for rarefaction (default 1000).
#' Repertoires with fewer sequences than this will be excluded with a warning.
#'
#' @return A tibble with one row per repertoire containing:
#' * `total_sequences`: Number of total sequences
#' * `unique_productive_sequences`: Number of unique clones
#' * `total_count`: Sum of UMI/read counts
#' * `clonality`: Shannon clonality (0 = diverse, 1 = monoclonal)
#' * `gini_coefficient`: Gini coefficient (0 = even, 1 = skewed)
#' * `simpson_index`: Simpson's D (0 = diverse, 1 = monoclonal)
#' * `inverse_simpson`: Effective number of dominant clones
#' * `top_productive_sequence`: Frequency (%) of most abundant clone
#' * `convergence`: Average nucleotide sequences per amino acid
#'
#' @details
#' Diversity Metrics:
#'
#' Shannon Clonality - Measures evenness of clone distribution. Calculated as
#' 1 - (entropy / log(unique clones)). Values near 0 indicate diverse repertoires;
#' near 1 indicate oligoclonal expansion.
#'
#' Gini Coefficient - Borrowed from economics to measure inequality. Based on
#' the Lorenz curve of cumulative clone frequencies. Ranges 0-1 where 0 is perfect
#' equality and 1 is maximal inequality (single dominant clone).
#'
#' Simpson Index - Probability that two randomly selected sequences belong to
#' the same clone. Higher values indicate lower diversity.
#'
#' Inverse Simpson - Number of equally-abundant clones needed to achieve the
#' observed diversity. More intuitive than Simpson's D (higher = more diverse).
#'
#' Rarefaction (rarefy = TRUE):
#'
#' When samples have different sequencing depths, raw diversity metrics are not
#' comparable. Rarefaction corrects this by: (1) Subsampling all repertoires to
#' the same depth (min_count), (2) Calculating diversity on the subsampled data,
#' (3) Repeating steps 1-2 for iterations, and (4) Averaging the results.
#'
#' This allows fair comparison between a deeply-sequenced blood sample and a
#' shallow tumor sample. Samples with fewer than min_count sequences are excluded.
#' @examples
#' file_path <- system.file("extdata", "TCRB_sequencing",
#'  package = "LymphoSeq2")
#' study_table <- LymphoSeq2::readImmunoSeq(path = file_path, threads = 1)
#' study_table <- LymphoSeq2::topSeqs(study_table, top = 100)
#' raw_clonality <- LymphoSeq2::clonality(study_table)
#' sampled_clonality <- LymphoSeq2::clonality(study_table,
#'   rarefy = TRUE,
#'   iterations = 100,
#'   min_count = 100
#' )
#' @seealso [LymphoSeq2::lorenzCurve()]
#' @export
clonality <- function(study_table,
                      rarefy = FALSE,
                      iterations = 100,
                      min_count = 1000) {
  if (rarefy) {
    low_count <- study_table |>
      dplyr::group_by(repertoire_id) |>
      dplyr::summarize(total = sum(duplicate_count)) |>
      dplyr::filter(total < min_count) |>
      dplyr::pull(repertoire_id)
    if (length(low_count) >= 1) {
      warning(stringr::str_c("Dropping the following samples since they have",
        "less than ", min_count, "sequences \n",
        stringr::str_c(low_count, sep = ""), sep = ""
      ))
    }
    study_table <- study_table |>
      dplyr::filter(!(repertoire_id %in% low_count)) |>
      dplyr::group_by(repertoire_id) |>
      dplyr::group_split() |>
      purrr::map(~ iterativeSummary(.x, iterations, min_count)) |>
      dplyr::bind_rows()
  } else {
    study_table <- study_table |>
      dplyr::group_by(repertoire_id) |>
      dplyr::group_split() |>
      purrr::map(summarySeq) |>
      dplyr::bind_rows()
  }

  return(study_table)
}

#' Get summary statistics for each repertoire_id in the analysis
#'
#' @inheritParams clonality
#' @return Tibble summarizing the sequence information for each repertoire_id
#'
#' @export

summarySeq <- function(study_table) {
  productive <- LymphoSeq2::productiveSeq(study_table, aggregate = "junction")
  frequency <- productive |>
    dplyr::pull(duplicate_frequency)
  counts <- productive |>
    dplyr::pull(duplicate_count)

  # Shannon entropy and clonality
  entropy <- -base::sum(frequency * base::log2(frequency), na.rm = TRUE)
  clonality <- 1 - base::round(entropy / base::log2(base::nrow(productive)),
    digits = 6)

  # Simpson index and inverse Simpson
  simpson <- calculate_simpson(frequency)
  inv_simpson <- calculate_inverse_simpson(frequency)

  convergence <- productive |>
    LymphoSeq2::topSeqs(top = 100) |>
    dplyr::group_by(junction_aa) |>
    dplyr::summarise(convergence = length(unique(junction))) |>
    dplyr::pull(convergence) |>
    mean()

  study_summary <- tibble::tibble(
    repertoire_id = study_table$repertoire_id[1],
    total_sequences = base::nrow(study_table),
    unique_productive_sequences = base::nrow(productive),
    total_count = base::sum(study_table$duplicate_count),
    clonality = clonality,
    gini_coefficient = calculate_gini(productive$duplicate_frequency),
    simpson_index = simpson,
    inverse_simpson = inv_simpson,
    top_productive_sequence = base::max((productive$duplicate_frequency) * 100),
    convergence = convergence
  )
  return(study_summary)
}

#' Calculate rarefied diversity metrics through repeated subsampling
#'
#' Normalize diversity estimates across samples by subsampling to a common depth.
#' This corrects for sequencing depth bias when comparing repertoires.
#'
#' @inheritParams clonality
#' @return Tibble with diversity metrics averaged across bootstrap iterations
#'
#' @export
iterativeSummary <- function(study_table, iterations, min_count = 1000) {
  # Convert to data.table for performance
  dt_table <- data.table::as.data.table(study_table)

  # Expand rows by duplicate_count (uncount operation)
  dt_expanded <- dt_table[rep(seq_len(.N), duplicate_count)]

  # Pre-allocate list for results
  summary_list <- vector("list", iterations)

  # Run iterations
  for (i in seq_len(iterations)) {
    summary_list[[i]] <- sampledSummary(dt_expanded, min_count)
  }

  # Combine and average results using data.table
  dt_combined <- data.table::rbindlist(summary_list)

  # Calculate means by repertoire_id
  numeric_cols <- setdiff(names(dt_combined), "repertoire_id")
  summary_table <- dt_combined[, lapply(.SD, mean), by = repertoire_id, .SDcols = numeric_cols]

  # Convert to tibble for user-facing output
  return(tibble::as_tibble(summary_table))
}


#' Sample repertoire to fixed depth and calculate diversity
#'
#' Randomly sample sequences to a specified depth and compute diversity metrics.
#' Used within bootstrap iterations for rarefied diversity estimation.
#'
#' @inheritParams clonality
#' @return Tibble with diversity metrics for the subsampled repertoire
#'
#' @export
sampledSummary <- function(study_table, min_count) {
  # Convert to data.table for performance
  dt_table <- data.table::as.data.table(study_table)

  # Check if table has enough rows
  n_rows <- nrow(dt_table)
  if (n_rows < min_count) {
    stop("Sample has ", n_rows, " sequences, but min_count is ", min_count,
         ". Cannot sample more sequences than available.", call. = FALSE)
  }

  # Sample using data.table (faster than dplyr::sample_n)
  sampled_idx <- sample(n_rows, min_count, replace = FALSE)
  dt_sampled <- dt_table[sampled_idx]

  # Aggregate by all columns except duplicate_frequency and duplicate_count
  # (we'll recompute duplicate_count from the sampled data)
  group_cols <- setdiff(names(dt_sampled), c("duplicate_frequency", "duplicate_count"))
  dt_aggregated <- dt_sampled[, .(duplicate_count = .N), by = group_cols]

  # Calculate frequency
  dt_aggregated[, duplicate_frequency := duplicate_count / sum(duplicate_count)]

  # Convert to tibble for summarySeq (which expects tibble)
  summary_table <- summarySeq(tibble::as_tibble(dt_aggregated))

  return(summary_table)
}

#' Calculate Gini coefficient
#'
#' Native implementation of the Gini coefficient calculation
#' (previously used ineq::Gini). The Gini coefficient measures
#' inequality in a distribution and ranges from 0 (perfect equality)
#' to 1 (maximal inequality).
#'
#' @param x Numeric vector of values (e.g., clone frequencies)
#' @return Numeric value between 0 and 1
#' @keywords internal
calculate_gini <- function(x) {
  # Remove NA values and sort
  x <- sort(x[!is.na(x)])
  n <- length(x)

  # Handle edge cases
  if (n == 0) return(NA_real_)
  if (n == 1) return(0)

  # Calculate Gini coefficient
  # G = (2 * sum(i * x_i)) / (n * sum(x_i)) - (n + 1) / n
  index <- seq_len(n)
  gini <- (2 * sum(index * x)) / (n * sum(x)) - (n + 1) / n

  return(gini)
}

#' Calculate Lorenz curve coordinates
#'
#' Native implementation of Lorenz curve calculation
#' (previously used ineq::Lc). Returns cumulative proportions
#' for creating Lorenz curves.
#'
#' @param x Numeric vector of values (e.g., clone frequencies)
#' @return List with components p (cumulative proportion of population)
#'   and L (cumulative proportion of values)
#' @keywords internal
calculate_lorenz <- function(x) {
  # Remove NA values and sort
  x <- sort(x[!is.na(x)])
  n <- length(x)

  # Handle edge case
  if (n == 0) {
    return(list(p = numeric(0), L = numeric(0)))
  }

  # Calculate cumulative proportions
  cumsum_x <- cumsum(x)
  total_x <- sum(x)

  # p: cumulative proportion of population (0 to 1)
  p <- c(0, seq_len(n) / n)

  # L: cumulative proportion of total values (0 to 1)
  L <- c(0, cumsum_x / total_x)

  return(list(p = p, L = L))
}

#' Calculate Simpson diversity index
#'
#' Compute the Simpson index (D), which represents the probability that two
#' randomly selected sequences from a repertoire belong to the same clone.
#'
#' @param x Numeric vector of clone frequencies (must sum to 1)
#' @return Numeric value between 0 and 1
#' @details
#' The Simpson index is calculated as D = sum(p_i^2) where p_i is the relative
#' frequency of clone i. This metric quantifies the probability of sampling
#' the same clone twice with replacement.
#'
#' Interpretation:
#' 0 = Infinite diversity (every clone equally frequent);
#' 1 = No diversity (single clone dominates);
#' 0.1 = Low clonality (highly diverse);
#' 0.5 = Moderate clonality (some dominant clones);
#' 0.9 = High clonality (oligoclonal repertoire).
#'
#' The Simpson index is less sensitive to rare clones than Shannon
#' entropy, making it useful for detecting oligoclonal expansions.
#'
#' @keywords internal
calculate_simpson <- function(x) {
  # Remove NA values
  x <- x[!is.na(x)]
  n <- length(x)

  # Handle edge cases
  if (n == 0) return(NA_real_)
  if (n == 1) return(1)

  # Calculate Simpson index: sum of squared frequencies
  simpson <- sum(x^2)

  return(simpson)
}

#' Calculate inverse Simpson diversity index
#'
#' Compute the effective number of dominant clones in a repertoire. This metric
#' transforms Simpson's D into an intuitive measure where higher values mean
#' greater diversity.
#'
#' @param x Numeric vector of clone frequencies (must sum to 1)
#' @return Numeric value >= 1 representing effective number of clones
#' @details
#' The inverse Simpson index is calculated as 1/D = 1 / sum(p_i^2), also known as
#' the Hill number of order 2. This represents the number of equally-abundant
#' clones that would produce the same diversity as the observed repertoire.
#'
#' Interpretation:
#' 1 = Monoclonal (single clone dominates);
#' 5-10 = Oligoclonal (few dominant clones);
#' 50-100 = Polyclonal (many clones contributing);
#' Greater than 100 = Highly diverse repertoire.
#'
#' Example: An inverse Simpson of 10 means the repertoire has the same
#' diversity as if it contained exactly 10 equally-abundant clones.
#'
#' Advantages over Simpson's D include: more intuitive interpretation
#' (higher = more diverse), interpretable as effective clone count,
#' and linear scale that is easier to compare across samples.
#'
#' @keywords internal
calculate_inverse_simpson <- function(x) {
  # Calculate Simpson index first
  simpson <- calculate_simpson(x)

  # Handle edge cases
  if (is.na(simpson) || simpson == 0) return(NA_real_)

  # Return inverse
  return(1 / simpson)
}
