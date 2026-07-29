# ------------------------------------------------------------
# HERVariumR compact user-facing feature tables
# ------------------------------------------------------------

.as_herv_bool <- function(x) {
  if (is.null(x)) {
    return(logical(0))
  }

  if (is.logical(x)) {
    x[is.na(x)] <- FALSE
    return(x)
  }

  z <- tolower(as.character(x))
  out <- z %in% c("true", "t", "1", "yes")
  out[is.na(out)] <- FALSE
  out
}

.numeric_or_zero <- function(x, n = NULL) {
  if (is.null(x)) {
    if (is.null(n)) return(numeric(0))
    return(rep(0, n))
  }

  out <- suppressWarnings(as.numeric(x))
  out[is.na(out)] <- 0
  out
}

.display_or_dot <- function(x, n = NULL) {
  if (is.null(x)) {
    if (is.null(n)) return(character(0))
    return(rep(".", n))
  }

  out <- as.character(x)
  out[is.na(out) | out == "" | out == "NA"] <- "."
  out
}

.first_existing_col <- function(features, cols) {
  cols <- cols[cols %in% colnames(features)]
  if (length(cols) == 0) return(NULL)
  cols[[1]]
}

.collapse_detected_encode_layers <- function(features, prefix) {
  n <- nrow(features)

  layer_map <- c(
    DNase = paste0(prefix, "_dnase_detected"),
    H3K27ac = paste0(prefix, "_h3k27ac_detected"),
    H3K4me3 = paste0(prefix, "_h3k4me3_detected"),
    transcription = paste0(prefix, "_transcription_detected")
  )

  out <- rep(".", n)

  for (i in seq_len(n)) {
    layers <- names(layer_map)[
      vapply(layer_map, function(col) {
        col %in% colnames(features) && .as_herv_bool(features[[col]])[i]
      }, logical(1))
    ]

    if (length(layers) > 0) {
      out[i] <- paste(layers, collapse = ";")
    }
  }

  out
}

.add_internal_transcription_direction_columns <- function(features) {
  n <- nrow(features)

  same <- if ("internal_same_strand_transcription_detected" %in% colnames(features)) {
    .as_herv_bool(features$internal_same_strand_transcription_detected)
  } else {
    rep(FALSE, n)
  }

  opposite <- if ("internal_opposite_strand_transcription_detected" %in% colnames(features)) {
    .as_herv_bool(features$internal_opposite_strand_transcription_detected)
  } else {
    rep(FALSE, n)
  }

  bidirectional <- if ("internal_bidirectional_transcription_detected" %in% colnames(features)) {
    .as_herv_bool(features$internal_bidirectional_transcription_detected)
  } else {
    same & opposite
  }

  if ("internal_same_strand_transcription_detected" %in% colnames(features) ||
      "internal_opposite_strand_transcription_detected" %in% colnames(features) ||
      "internal_bidirectional_transcription_detected" %in% colnames(features)) {
    features$internal_encode_same_strand_transcription_detected <- same
    features$internal_encode_opposite_strand_transcription_detected <- opposite
    features$internal_encode_bidirectional_transcription_detected <- bidirectional
  }

  features
}

.add_ltr_transcription_direction_columns <- function(features, prefix) {
  n <- nrow(features)

  same_col <- paste0(prefix, "_same_strand_transcription_detected")
  opposite_col <- paste0(prefix, "_opposite_strand_transcription_detected")
  bidirectional_col <- paste0(prefix, "_bidirectional_transcription_detected")
  class_col <- paste0(prefix, "_transcription_direction_class")

  if (!any(c(same_col, opposite_col, class_col) %in% colnames(features))) {
    return(features)
  }

  same <- if (same_col %in% colnames(features)) {
    .as_herv_bool(features[[same_col]])
  } else {
    rep(FALSE, n)
  }

  opposite <- if (opposite_col %in% colnames(features)) {
    .as_herv_bool(features[[opposite_col]])
  } else {
    rep(FALSE, n)
  }

  bidirectional <- same & opposite

  if (class_col %in% colnames(features)) {
    class <- tolower(as.character(features[[class_col]]))
    bidirectional <- bidirectional | grepl("bidirectional|both", class)
    bidirectional[is.na(bidirectional)] <- FALSE
  }

  features[[same_col]] <- same
  features[[opposite_col]] <- opposite
  features[[bidirectional_col]] <- bidirectional

  features
}

.add_compact_derived_columns <- function(features) {
  if (!is.data.frame(features) || nrow(features) == 0) {
    return(features)
  }

  n <- nrow(features)

  if (any(c(
    "internal_dnase_detected",
    "internal_h3k27ac_detected",
    "internal_h3k4me3_detected",
    "internal_transcription_detected"
  ) %in% colnames(features))) {
    features$internal_encode_detected_layers <-
      .collapse_detected_encode_layers(features, "internal")
  }

  if (any(c(
    "ltr5_encode_dnase_detected",
    "ltr5_encode_h3k27ac_detected",
    "ltr5_encode_h3k4me3_detected",
    "ltr5_encode_transcription_detected"
  ) %in% colnames(features))) {
    features$ltr5_encode_detected_layers <-
      .collapse_detected_encode_layers(features, "ltr5_encode")
  }

  if (any(c(
    "ltr3_encode_dnase_detected",
    "ltr3_encode_h3k27ac_detected",
    "ltr3_encode_h3k4me3_detected",
    "ltr3_encode_transcription_detected"
  ) %in% colnames(features))) {
    features$ltr3_encode_detected_layers <-
      .collapse_detected_encode_layers(features, "ltr3_encode")
  }

  features <- .add_internal_transcription_direction_columns(features)
  features <- .add_ltr_transcription_direction_columns(features, "ltr5_encode")
  features <- .add_ltr_transcription_direction_columns(features, "ltr3_encode")

  if (all(c("ltr5_encode_detected_layer_count", "ltr5_encode_ccre_overlap") %in% colnames(features))) {
    features$ltr5_encode_evidence_including_ccre <-
      .numeric_or_zero(features$ltr5_encode_detected_layer_count, n) > 0 |
      .as_herv_bool(features$ltr5_encode_ccre_overlap)
  }

  if (all(c("ltr3_encode_detected_layer_count", "ltr3_encode_ccre_overlap") %in% colnames(features))) {
    features$ltr3_encode_evidence_including_ccre <-
      .numeric_or_zero(features$ltr3_encode_detected_layer_count, n) > 0 |
      .as_herv_bool(features$ltr3_encode_ccre_overlap)
  }

  # Normalize optional list columns when they exist. Detailed TFBM lists are
  # only exposed when their optional layer was enabled.
  list_cols <- c(
    "ltr5_tfbm_tf_names_all",
    "ltr3_tfbm_tf_names_all"
  )

  for (col in list_cols[list_cols %in% colnames(features)]) {
    features[[col]] <- .display_or_dot(features[[col]], n)
  }

  # Create one canonical, non-truncated TF ChIP-seq TF list per LTR. Prefer the
  # explicit full column from the compact TF ChIP-seq asset, then the canonical
  # list column, and use the display column only as a last-resort fallback for
  # older/custom assets.
  ltr5_tfchip_source <- .first_existing_col(
    features,
    c("ltr5_tfchip_tf_list_full", "ltr5_tfchip_tf_list", "ltr5_tfchip_tf_list_display")
  )
  if (!is.null(ltr5_tfchip_source)) {
    features$ltr5_tfchip_tf_list <- .display_or_dot(features[[ltr5_tfchip_source]], n)
  }

  ltr3_tfchip_source <- .first_existing_col(
    features,
    c("ltr3_tfchip_tf_list_full", "ltr3_tfchip_tf_list", "ltr3_tfchip_tf_list_display")
  )
  if (!is.null(ltr3_tfchip_source)) {
    features$ltr3_tfchip_tf_list <- .display_or_dot(features[[ltr3_tfchip_source]], n)
  }

  features
}

select_compact_herv_features <- function(features, include_comparison_group = FALSE) {
  if (!is.data.frame(features)) {
    stop("features must be a data.frame.")
  }

  features <- .add_compact_derived_columns(features)

  compact_cols <- c(
    "HERV_id",
    "subfamily",
    "locid",
    "chrom",
    "start",
    "end",
    "strand",

    "has_domain",
    "domain_count",
    "max_domain_cov",
    "has_gag",
    "has_pol",
    "has_env",
    "has_accessory",
    "has_complete_gag_pol_env",
    "domains_type",

    "has_ltr5",
    "has_ltr3",
    "has_both_ltrs",
    "ltr5_name",
    "ltr3_name",
    "ltr5_tfbm_burden",
    "ltr3_tfbm_burden",

    # Optional detailed TFBM columns. These appear only when
    # add_tfbm_details = TRUE and a detailed FIMO file was supplied.
    "ltr5_n_unique_tfbm_tf_names_detailed",
    "ltr5_tfbm_tf_names_all",
    "ltr3_n_unique_tfbm_tf_names_detailed",
    "ltr3_tfbm_tf_names_all",

    "rbp_burden",
    "rbp_unique",

    "has_ltr5_stat1_motif",
    "has_ltr3_stat1_motif",
    "has_ltr5_stat1stat2_motif",
    "has_ltr3_stat1stat2_motif",

    "has_any_terminal_exon_domain",
    "has_terminal_exon_domain_protein_coding",
    "has_terminal_exon_domain_lncRNA",
    "max_terminal_domain_coverage_protein_coding",
    "max_terminal_domain_coverage_lncRNA",

    "feature_overlap",
    "ov_gene_types",
    "is_intergenic",
    "overlaps_gene",
    "overlaps_exon",
    "overlaps_intron",
    "overlaps_cds",
    "overlaps_utr",
    "overlaps_lncRNA",
    "overlaps_protein_coding",

    # Keep one umbrella internal-evidence flag. The detected-layer count and
    # readable layer list retain the more informative multilayer detail.
    "has_internal_encode_evidence",
    "internal_encode_detected_layer_count",
    "internal_encode_detected_layers",
    "internal_encode_same_strand_transcription_detected",
    "internal_encode_opposite_strand_transcription_detected",
    "internal_encode_bidirectional_transcription_detected",

    "ltr5_encode_evidence_including_ccre",
    "ltr3_encode_evidence_including_ccre",
    "ltr5_encode_detected_layer_count",
    "ltr3_encode_detected_layer_count",
    "ltr5_encode_detected_layers",
    "ltr3_encode_detected_layers",
    "ltr5_encode_same_strand_transcription_detected",
    "ltr5_encode_opposite_strand_transcription_detected",
    "ltr5_encode_bidirectional_transcription_detected",
    "ltr3_encode_same_strand_transcription_detected",
    "ltr3_encode_opposite_strand_transcription_detected",
    "ltr3_encode_bidirectional_transcription_detected",
    "ltr5_encode_ccre_overlap",
    "ltr3_encode_ccre_overlap",
    "ltr5_encode_ccre_n_overlaps",
    "ltr3_encode_ccre_n_overlaps",

    "has_ltr5_tfchip_overlap",
    "has_ltr3_tfchip_overlap",
    "ltr5_tfchip_n_peak_clusters",
    "ltr3_tfchip_n_peak_clusters",
    "ltr5_tfchip_n_tfs",
    "ltr3_tfchip_n_tfs",
    "ltr5_tfchip_tf_list",
    "ltr3_tfchip_tf_list"
  )

  if (include_comparison_group) {
    compact_cols <- c(compact_cols, "comparison_group")
  }

  compact_cols <- compact_cols[compact_cols %in% colnames(features)]
  features[, compact_cols, drop = FALSE]
}
