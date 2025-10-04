#' Merge paired TCR chains by cell
#'
#' Combines alpha and beta TCR chains from the same cell into a single observation,
#' creating paired-chain clonotypes for single-cell VDJ data. Uses conservative
#' counting to avoid inflating clone frequencies.
#'
#' @param airr_data A tibble with AIRR-formatted data from [read10x()]
#' @param mode Merging mode:
#'  * `"none"` (default): Keep all chains separate (no merging)
#'  * `"strict"`: Only keep cells with exactly one alpha and one beta chain
#'  * `"best"`: Select the most frequent alpha and beta chains per cell
#' @param chain_separator Character to use between chain sequences. Default: "|"
#'
#' @return A tibble with one row per cell containing both chains. The output
#' includes:
#'  * `cell_id`: Cell barcode
#'  * `repertoire_id`: Sample identifier
#'  * `junction`: Combined nucleotide CDR3 (alpha|beta)
#'  * `junction_aa`: Combined amino acid CDR3 (alpha|beta)
#'  * `v_call_alpha`, `v_call_beta`: Separate V gene calls
#'  * `j_call_alpha`, `j_call_beta`: Separate J gene calls
#'  * `duplicate_count`: Minimum count between the two chains
#'  * `duplicate_frequency`: Recalculated frequency
#'  * `paired_clonotype`: Combined identifier for the pair
#'
#' @details
#' This function pairs TRA and TRB chains from the same cell while maintaining
#' compatibility with downstream LymphoSeq2 functions. Key features:
#'
#' **Conservative Counting**: Uses the minimum `duplicate_count` between paired
#' chains to avoid overestimating clone frequency. This is important because
#' UMI counts can differ between chains.
#'
#' **Chain Identification**: Stores individual chain information in separate
#' columns (alpha/beta) so gene usage and other metrics can still be calculated.
#'
#' **Combined Sequences**: The `junction` and `junction_aa` fields contain both
#' chains separated by "|" for use in diversity and similarity analyses.
#'
#' **Mode Selection**:
#' - `strict`: Most stringent - requires exactly 1 alpha + 1 beta
#' - `best`: More permissive - selects top chain when multiples exist
#' - `none`: No merging (returns input unchanged)
#'
#' @examples
#' \dontrun{
#' # Read 10X AIRR data
#' sc_data <- read10x("path/to/airr.tsv")
#'
#' # Strict pairing
#' paired <- merge_chains(sc_data, mode = "strict")
#'
#' # Use with standard analyses
#' diversity <- clonality(paired)
#' top_clones <- topSeqs(paired, top = 100)
#'
#' # Gene frequency on alpha chains
#' alpha_genes <- paired %>%
#'   select(repertoire_id, duplicate_count, v_call = v_call_alpha) %>%
#'   geneFreq(locus = "V")
#' }
#'
#' @seealso [read10x()], [clonality()], [productiveSeq()]
#' @export
merge_chains <- function(airr_data, mode = "none", chain_separator = "|") {
  if (!mode %in% c("none", "strict", "best")) {
    stop("mode must be one of: 'none', 'strict', or 'best'")
  }

  if (mode == "none") {
    return(airr_data)
  }

  # Extract chain type from v_call
  data_with_chain <- airr_data |>
    dplyr::mutate(
      chain_type = stringr::str_extract(v_call, "TR[ABGD]|TCR[ABGD]"),
      chain_type = stringr::str_replace(chain_type, "TCR", "TR"),
      chain_type = dplyr::coalesce(chain_type, "Unknown")
    )

  # Use data.table for performance
  dt_data <- data_with_chain |>
    dtplyr::lazy_dt()
  # Filter and select chains based on mode
  if (mode == "strict") {
    # Only cells with exactly 1 TRA and 1 TRB
    valid_cells <- dt_data |>
      tibble::as_tibble() |>
      dplyr::filter(chain_type %in% c("TRA", "TRB")) |>
      dplyr::group_by(cell_id, chain_type) |>
      dplyr::summarise(n = dplyr::n(), .groups = "drop_last") |>
      dplyr::filter(n == 1) |>
      dplyr::summarise(n_chains = dplyr::n(), .groups = "drop") |>
      dplyr::filter(n_chains == 2) |>
      dplyr::pull(cell_id)
    
    paired_data <- dt_data |>
      tibble::as_tibble() |>
      dplyr::filter(cell_id %in% valid_cells, chain_type %in% c("TRA", "TRB")) |>
      dplyr::as_tibble()

  } else if (mode == "best") {
    # Select most frequent chain of each type using dplyr
    paired_data <- data_with_chain |>
      tibble::as_tibble() |>
      dplyr::filter(chain_type %in% c("TRA", "TRB")) |>
      dplyr::group_by(cell_id, chain_type) |>
      dplyr::slice_max(order_by = duplicate_count, n = 1, with_ties = FALSE) |>
      dplyr::ungroup() |>
      dplyr::group_by(cell_id) |>
      dplyr::filter(dplyr::n() == 2) |>
      dplyr::ungroup()
  }

  # Convert to tibble first to ensure consistent data structure
  paired_data <- tibble::as_tibble(paired_data)

  # Ensure required columns exist with defaults
  if (!"d_call" %in% names(paired_data)) {
    paired_data$d_call <- NA_character_
  }
  if (!"j_call" %in% names(paired_data)) {
    paired_data$j_call <- NA_character_
  }
  if (!"c_call" %in% names(paired_data)) {
    paired_data$c_call <- NA_character_
  }
  if (!"productive" %in% names(paired_data)) {
    paired_data$productive <- TRUE
  }

  # Separate alpha and beta chains
  alpha_chains <- paired_data |>
    dplyr::filter(chain_type == "TRA") |>
    dplyr::select(
      cell_id, repertoire_id,
      v_call_alpha = v_call,
      d_call_alpha = d_call,
      j_call_alpha = j_call,
      c_call_alpha = c_call,
      junction_alpha = junction,
      junction_aa_alpha = junction_aa,
      count_alpha = duplicate_count,
      productive_alpha = productive
    )

  beta_chains <- paired_data |>
    dplyr::filter(chain_type == "TRB") |>
    dplyr::select(
      cell_id,
      v_call_beta = v_call,
      d_call_beta = d_call,
      j_call_beta = j_call,
      c_call_beta = c_call,
      junction_beta = junction,
      junction_aa_beta = junction_aa,
      count_beta = duplicate_count,
      productive_beta = productive
    )

  # Join chains by cell_id
  merged <- dplyr::inner_join(alpha_chains, beta_chains, by = "cell_id")

  # Create combined fields
  result <- merged |>
    tibble::as_tibble() |>
    dplyr::mutate(
      # Combined sequences using separator
      junction = paste(junction_alpha, junction_beta, sep = chain_separator),
      junction_aa = paste(junction_aa_alpha, junction_aa_beta, sep = chain_separator),

      # Paired clonotype identifier
      paired_clonotype = junction_aa,

      # Use minimum count (conservative)
      duplicate_count = pmin(count_alpha, count_beta),

      # Productive only if both chains are productive
      productive = productive_alpha & productive_beta,

      # Junction lengths (combined)
      junction_length = stringr::str_length(junction),
      junction_aa_length = stringr::str_length(junction_aa),

      # Reading frame - in-frame if both chains are productive
      reading_frame = dplyr::if_else(productive, "in-frame", "out-of-frame"),

      # V and J calls - keep both for downstream analysis
      # For functions expecting single v_call, they can use v_call_alpha or v_call_beta
      v_call = v_call_alpha,  # Default to alpha for compatibility
      j_call = j_call_alpha,
      d_call = NA_character_,  # Not applicable for paired

      # Extract gene families (needed for some downstream functions)
      v_family = dplyr::if_else(
        stringr::str_detect(v_call_alpha, "(TRA|TCRA)V"),
        stringr::str_extract(v_call_alpha, "(TRA|TCRA)V\\d+"),
        "unrecognized"
      ),
      j_family = dplyr::if_else(
        stringr::str_detect(j_call_alpha, "(TRA|TCRA)J"),
        stringr::str_extract(j_call_alpha, "(TRA|TCRA)J\\d+"),
        "unrecognized"
      ),
      d_family = NA_character_  # Not applicable for TRA
    ) |>
    dplyr::group_by(repertoire_id) |>
    dplyr::mutate(
      duplicate_frequency = duplicate_count / sum(duplicate_count)
    ) |>
    dplyr::ungroup() |>
    dplyr::select(
      cell_id, repertoire_id,
      junction, junction_aa,
      v_call, j_call, d_call,  # For compatibility
      v_family, j_family, d_family,  # Gene families
      v_call_alpha, j_call_alpha, d_call_alpha, c_call_alpha,
      v_call_beta, j_call_beta, d_call_beta, c_call_beta,
      junction_alpha, junction_aa_alpha,
      junction_beta, junction_aa_beta,
      junction_length, junction_aa_length,
      duplicate_count, duplicate_frequency,
      productive, reading_frame,
      paired_clonotype
    )

  return(result)
}
