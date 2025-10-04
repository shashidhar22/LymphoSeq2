#' Generate rarefaction curves for diversity estimation
#'
#' Estimate repertoire diversity across different sequencing depths using
#' rarefaction (interpolation) and extrapolation. This helps determine if
#' sequencing depth is sufficient to capture repertoire diversity.
#'
#' @param sample_table A tibble from [readImmunoSeq()] containing "junction_aa",
#' "duplicate_count", and "duplicate_frequency" columns. Use a single repertoire
#' (filter to one repertoire_id before calling this function).
#' @param q Diversity order to calculate:
#'   * 0 (default): Species richness (number of unique clones)
#'   * 1: Shannon diversity (accounts for evenness)
#'   * 2: Simpson diversity (emphasizes abundant clones)
#' @param endpoint Maximum sequencing depth for extrapolation. Default is 2x
#' the observed sample size. Set higher to predict diversity at deeper sequencing.
#' @param nboot Number of bootstrap iterations for confidence intervals (default 50).
#' Higher values give more precise estimates but take longer.
#' @param conf Confidence level for intervals (default 0.95 for 95% CI)
#'
#' @return A tibble with rarefaction/extrapolation results:
#' * `m`: Sample size (sequencing depth)
#' * `method`: "interpolated", "observed", or "extrapolated"
#' * `order`: Diversity order (same as `q` parameter)
#' * `qD`: Estimated diversity at depth `m`
#' * `qD.LCL`: Lower confidence limit
#' * `qD.UCL`: Upper confidence limit
#' * `SC`: Standard error from bootstrap
#' * `repertoire_id`: Sample identifier
#'
#' @details
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
#' Confidence intervals: Computed via bootstrap for rarefaction points only.
#' Wider intervals indicate more uncertainty in the estimate.
#' @examples
#' file_path <- system.file("extdata", "TCRB_sequencing",
#'  package = "LymphoSeq2")
#' study_table <- LymphoSeq2::readImmunoSeq(path = file_path, threads = 1) |>
#'   LymphoSeq2::topSeqs(top = 100)
#' amino_table <- LymphoSeq2::productiveSeq(study_table,
#'   aggregate = "junction_aa",
#'   prevalence = TRUE
#' )
#' amino_table <- amino_table |>
#'   dplyr::filter(repertoire_id == "TRB_Unsorted_1320")
#' rarefaction_table <- LymphoSeq2::runINext(amino_table)
#' @export
runINext <- function(sample_table, q = 0, endpoint = NULL, nboot = 50, conf = 0.95) {
  repertoire_id <- sample_table$repertoire_id[1]

  # Get abundance vector
  abundance <- sample_table$duplicate_count
  n <- sum(abundance)  # Total sample size

  # Set endpoint if not provided (2x observed size)
  if (is.null(endpoint)) {
    endpoint <- n * 2
  }

  # Generate sample sizes for rarefaction/extrapolation
  sample_sizes <- unique(c(
    seq(1, n, length.out = 40),  # Interpolation
    seq(n, endpoint, length.out = 40)  # Extrapolation
  ))
  sample_sizes <- sort(sample_sizes)

  # Calculate diversity at each sample size
  result_list <- lapply(sample_sizes, function(m) {
    if (m <= n) {
      # Rarefaction (interpolation)
      est <- rarefaction_diversity(abundance, m, q)
    } else {
      # Extrapolation
      est <- extrapolation_diversity(abundance, m, q)
    }

    # Bootstrap CI if requested
    if (nboot > 0 && m <= n) {
      boot_vals <- replicate(nboot, {
        boot_sample <- sample(rep(1:length(abundance), abundance),
                               size = n, replace = TRUE)
        boot_abundance <- table(factor(boot_sample, levels = 1:length(abundance)))
        rarefaction_diversity(as.numeric(boot_abundance), m, q)
      })
      se <- sd(boot_vals)
      ci_lower <- est - qnorm(1 - (1 - conf)/2) * se
      ci_upper <- est + qnorm(1 - (1 - conf)/2) * se
    } else {
      se <- NA
      ci_lower <- NA
      ci_upper <- NA
    }

    data.frame(
      m = m,
      method = ifelse(m < n, "interpolated",
                      ifelse(m == n, "observed", "extrapolated")),
      order = q,
      qD = est,
      qD.LCL = ci_lower,
      qD.UCL = ci_upper,
      SC = se,
      repertoire_id = repertoire_id
    )
  })

  result <- dplyr::bind_rows(result_list)
  return(tibble::as_tibble(result))
}

#' Calculate rarefaction diversity (interpolation)
#'
#' @param abundance Vector of species abundances
#' @param m Sample size
#' @param q Diversity order
#' @return Diversity estimate
#' @keywords internal
rarefaction_diversity <- function(abundance, m, q) {
  n <- sum(abundance)

  if (m == n) {
    # Observed diversity
    if (q == 0) {
      return(sum(abundance > 0))  # Species richness
    } else if (q == 1) {
      p <- abundance[abundance > 0] / n
      return(exp(-sum(p * log(p))))  # Shannon
    } else if (q == 2) {
      p <- abundance[abundance > 0] / n
      return(1 / sum(p^2))  # Inverse Simpson
    }
  }

  if (q == 0) {
    # Species richness rarefaction
    S_rare <- sum(1 - sapply(abundance[abundance > 0], function(x) {
      if (n - x < m) return(0)
      exp(lchoose(n - x, m) - lchoose(n, m))
    }))
    return(S_rare)
  } else {
    # For q > 0, use sampling approach
    # This is a simplified version
    p_rare <- sapply(abundance[abundance > 0], function(x) {
      x * (1 - (1 - x/n)^m)
    }) / m

    if (q == 1) {
      p_rare <- p_rare[p_rare > 0]
      return(exp(-sum(p_rare * log(p_rare))))
    } else if (q == 2) {
      return(1 / sum(p_rare^2))
    }
  }
}

#' Calculate extrapolation diversity
#'
#' @param abundance Vector of species abundances
#' @param m Sample size (> observed)
#' @param q Diversity order
#' @return Diversity estimate
#' @keywords internal
extrapolation_diversity <- function(abundance, m, q) {
  n <- sum(abundance)
  f1 <- sum(abundance == 1)  # Singletons
  f2 <- sum(abundance == 2)  # Doubletons

  if (q == 0) {
    # Chao1-based extrapolation for species richness
    S_obs <- sum(abundance > 0)

    # Chao1 estimator
    if (f2 > 0) {
      S_chao <- S_obs + f1^2 / (2 * f2)
    } else if (f1 > 0) {
      S_chao <- S_obs + f1 * (f1 - 1) / 2
    } else {
      S_chao <- S_obs
    }

    # Linear extrapolation
    slope <- (S_chao - S_obs) / n
    S_extrap <- S_obs + slope * (m - n)

    return(max(S_obs, S_extrap))
  } else {
    # For q > 0, use asymptotic approach
    p <- abundance[abundance > 0] / n

    if (q == 1) {
      return(exp(-sum(p * log(p))))
    } else if (q == 2) {
      return(1 / sum(p^2))
    }
  }
}
