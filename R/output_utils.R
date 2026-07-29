# ------------------------------------------------------------
# HERVariumR output-simplification utilities
# ------------------------------------------------------------

.match_output_level <- function(output_level = c("standard", "extended")) {
  match.arg(output_level)
}

.prepare_output_layout <- function(output_dir,
                                   output_level = c("standard", "extended")) {
  output_level <- .match_output_level(output_level)

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  details_dir <- file.path(output_dir, "details")
  plots_dir <- file.path(output_dir, "plots")

  if (identical(output_level, "extended")) {
    if (!dir.exists(details_dir)) dir.create(details_dir, recursive = TRUE)
    if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
  }

  list(
    output_level = output_level,
    output_dir = output_dir,
    details_dir = details_dir,
    plots_dir = plots_dir,
    write_details = identical(output_level, "extended")
  )
}

.write_hervarium_tsv <- function(x, file) {
  if (is.null(x)) {
    x <- data.frame()
  }

  write.table(
    x,
    file = file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  invisible(file)
}

.write_detail_tsv <- function(x, filename, layout) {
  if (!isTRUE(layout$write_details)) {
    return(invisible(NULL))
  }

  if (!dir.exists(layout$details_dir)) {
    dir.create(layout$details_dir, recursive = TRUE)
  }

  .write_hervarium_tsv(x, file.path(layout$details_dir, filename))
}

.make_transient_layer_dir <- function(layout, label) {
  if (isTRUE(layout$write_details)) {
    path <- file.path(layout$details_dir, label)
    if (!dir.exists(path)) dir.create(path, recursive = TRUE)
    return(path)
  }

  path <- tempfile(pattern = paste0("hervarium_", label, "_"))
  dir.create(path, recursive = TRUE)
  path
}

.clean_result_feature_label <- function(x) {
  z <- as.character(x)

  replacements <- c(
    "has_internal_robust_multilayer_encode_evidence" = "Internal: robust multilayer evidence",
    "has_internal_multilayer_encode_evidence" = "Internal: multilayer evidence",
    "has_internal_encode_evidence" = "Internal: external evidence",
    "ltr5_encode_evidence_including_ccre" = "5\u2032 LTR: external evidence incl. cCRE",
    "ltr3_encode_evidence_including_ccre" = "3\u2032 LTR: external evidence incl. cCRE",
    "ltr5_encode_ccre_overlap" = "5\u2032 LTR: cCRE overlap",
    "ltr3_encode_ccre_overlap" = "3\u2032 LTR: cCRE overlap",
    "has_ltr5_tfchip_overlap" = "5\u2032 LTR: TF ChIP-seq overlap",
    "has_ltr3_tfchip_overlap" = "3\u2032 LTR: TF ChIP-seq overlap",
    "has_any_terminal_exon_domain" = "Terminal-exon domain",
    "has_terminal_exon_domain_protein_coding" = "Terminal-exon domain: protein-coding",
    "has_terminal_exon_domain_lncRNA" = "Terminal-exon domain: lncRNA",
    "has_both_ltrs" = "Both LTRs present",
    "has_ltr5" = "5\u2032 LTR present",
    "has_ltr3" = "3\u2032 LTR present",
    "domain_count" = "Internal domain count",
    "max_domain_cov" = "Maximum internal-domain coverage",
    "ltr5_tfbm_burden" = "5\u2032 LTR TFBM burden",
    "ltr3_tfbm_burden" = "3\u2032 LTR TFBM burden",
    "rbp_burden" = "RBP motif burden",
    "rbp_unique" = "Unique RBP motifs"
  )

  exact <- z %in% names(replacements)
  z[exact] <- unname(replacements[z[exact]])

  z <- gsub("^internal_encode_", "Internal: ", z)
  z <- gsub("^internal_", "Internal: ", z)
  z <- gsub("^ltr5_encode_", "5\u2032 LTR: ", z)
  z <- gsub("^ltr3_encode_", "3\u2032 LTR: ", z)
  z <- gsub("^ltr5_tfchip_", "5\u2032 LTR TF ChIP-seq: ", z)
  z <- gsub("^ltr3_tfchip_", "3\u2032 LTR TF ChIP-seq: ", z)
  z <- gsub("^ltr5_", "5\u2032 LTR: ", z)
  z <- gsub("^ltr3_", "3\u2032 LTR: ", z)
  z <- gsub("^has_", "", z)
  z <- gsub("_", " ", z, fixed = TRUE)

  z <- gsub("encode", "external evidence", z, ignore.case = TRUE)
  z <- gsub("ccre", "cCRE", z, ignore.case = TRUE)
  z <- gsub("tfchip", "TF ChIP-seq", z, ignore.case = TRUE)
  z <- gsub("tfbm", "TFBM", z, ignore.case = TRUE)
  z <- gsub("rbp", "RBP", z, ignore.case = TRUE)
  z <- gsub("stat1stat2", "STAT1::STAT2", z, ignore.case = TRUE)
  z <- gsub("stat1", "STAT1", z, ignore.case = TRUE)
  z <- gsub("h3k27ac", "H3K27ac", z, ignore.case = TRUE)
  z <- gsub("h3k4me3", "H3K4me3", z, ignore.case = TRUE)
  z <- gsub("dnase", "DNase", z, ignore.case = TRUE)
  z <- gsub("protein coding", "protein-coding", z, fixed = TRUE)
  z <- gsub(" +", " ", z)
  trimws(z)
}

.infer_result_layer <- function(feature) {
  f <- tolower(as.character(feature))

  if (grepl("subfamily", f)) return("subfamily")
  if (grepl("tfchip.*tf|tf.*tfchip", f)) return("tfchip_tf")
  if (grepl("tfchip", f)) return("tfchip")
  if (grepl("tfbm", f)) return("tfbm")
  if (grepl("rbp", f)) return("rbp")
  if (grepl("stat1|stat2|ifn|irf", f)) return("ifn_ltr")
  if (grepl("terminal", f)) return("terminal_exon")
  if (grepl("encode|ccre|dnase|h3k27ac|h3k4me3|transcription", f)) {
    return("external_evidence")
  }
  if (grepl("domain|gag|pol|env|accessory|protease|integrase|rnase|dutpase|chromodomain", f)) {
    return("internal_domains")
  }
  if (grepl("ltr5|ltr3|both_ltr|both ltr", f)) return("ltr_structure")
  if (grepl("intergenic|overlaps|exon|intron|utr|cds|lncrna|protein_coding|protein-coding|gene", f)) {
    return("transcript_context")
  }

  "core_identity"
}

.order_result_rows <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)

  padj_sort <- suppressWarnings(as.numeric(df$padj))
  p_sort <- suppressWarnings(as.numeric(df$p_value))
  effect_sort <- abs(suppressWarnings(as.numeric(df$effect_size)))

  padj_sort[is.na(padj_sort)] <- Inf
  p_sort[is.na(p_sort)] <- Inf
  effect_sort[is.na(effect_sort)] <- -Inf

  df[order(padj_sort, p_sort, -effect_sort), , drop = FALSE]
}

.binary_results_to_key <- function(df,
                                   layer = NULL,
                                   test_type = "binary_fisher",
                                   max_rows = NULL,
                                   label_prefix = NULL,
                                   pseudocount = 0.5) {
  if (is.null(df) || nrow(df) == 0) return(data.frame())

  required <- c(
    "feature", "foreground_yes", "foreground_no",
    "background_yes", "background_no", "p_value", "padj", "direction"
  )
  if (!all(required %in% colnames(df))) return(data.frame())

  df <- df[!is.na(df$p_value), , drop = FALSE]
  if (nrow(df) == 0) return(data.frame())

  a <- suppressWarnings(as.numeric(df$foreground_yes)) + pseudocount
  b <- suppressWarnings(as.numeric(df$foreground_no)) + pseudocount
  c <- suppressWarnings(as.numeric(df$background_yes)) + pseudocount
  d <- suppressWarnings(as.numeric(df$background_no)) + pseudocount
  effect <- log2((a / b) / (c / d))

  fg_total <- suppressWarnings(as.numeric(df$foreground_yes)) +
    suppressWarnings(as.numeric(df$foreground_no))
  bg_total <- suppressWarnings(as.numeric(df$background_yes)) +
    suppressWarnings(as.numeric(df$background_no))

  fg_value <- ifelse(
    fg_total > 0,
    100 * suppressWarnings(as.numeric(df$foreground_yes)) / fg_total,
    NA_real_
  )
  bg_value <- ifelse(
    bg_total > 0,
    100 * suppressWarnings(as.numeric(df$background_yes)) / bg_total,
    NA_real_
  )

  feature_id <- as.character(df$feature)
  feature_label <- .clean_result_feature_label(feature_id)

  if (!is.null(label_prefix)) {
    feature_id <- paste(label_prefix, feature_id, sep = ":")
    feature_label <- paste(label_prefix, feature_label, sep = ": ")
  }

  inferred_layer <- if (is.null(layer)) {
    vapply(feature_id, .infer_result_layer, character(1))
  } else {
    rep(layer, nrow(df))
  }

  out <- data.frame(
    layer = inferred_layer,
    feature = feature_id,
    feature_label = feature_label,
    test_type = test_type,
    effect_type = "log2_odds_ratio",
    foreground_value = fg_value,
    background_value = bg_value,
    effect_size = effect,
    p_value = suppressWarnings(as.numeric(df$p_value)),
    padj = suppressWarnings(as.numeric(df$padj)),
    direction = as.character(df$direction),
    significant = !is.na(df$padj) & df$padj < 0.05,
    stringsAsFactors = FALSE
  )

  out <- .order_result_rows(out)
  if (!is.null(max_rows) && is.finite(max_rows)) out <- head(out, max_rows)
  out
}

.numeric_results_to_key <- function(df,
                                    layer = NULL,
                                    test_type = "numeric_wilcox",
                                    max_rows = NULL,
                                    label_prefix = NULL) {
  if (is.null(df) || nrow(df) == 0) return(data.frame())

  required <- c(
    "feature", "foreground_median", "background_median",
    "delta_median", "p_value", "padj", "direction"
  )
  if (!all(required %in% colnames(df))) return(data.frame())

  df <- df[!is.na(df$p_value), , drop = FALSE]
  if (nrow(df) == 0) return(data.frame())

  feature_id <- as.character(df$feature)
  feature_label <- .clean_result_feature_label(feature_id)

  if (!is.null(label_prefix)) {
    feature_id <- paste(label_prefix, feature_id, sep = ":")
    feature_label <- paste(label_prefix, feature_label, sep = ": ")
  }

  inferred_layer <- if (is.null(layer)) {
    vapply(feature_id, .infer_result_layer, character(1))
  } else {
    rep(layer, nrow(df))
  }

  out <- data.frame(
    layer = inferred_layer,
    feature = feature_id,
    feature_label = feature_label,
    test_type = test_type,
    effect_type = "delta_median",
    foreground_value = suppressWarnings(as.numeric(df$foreground_median)),
    background_value = suppressWarnings(as.numeric(df$background_median)),
    effect_size = suppressWarnings(as.numeric(df$delta_median)),
    p_value = suppressWarnings(as.numeric(df$p_value)),
    padj = suppressWarnings(as.numeric(df$padj)),
    direction = as.character(df$direction),
    significant = !is.na(df$padj) & df$padj < 0.05,
    stringsAsFactors = FALSE
  )

  out <- .order_result_rows(out)
  if (!is.null(max_rows) && is.finite(max_rows)) out <- head(out, max_rows)
  out
}

build_comparison_key_results <- function(binary_results,
                                         numeric_results,
                                         subfamily_results = NULL,
                                         domain_type_results = NULL,
                                         tfbm_ltr5_results = NULL,
                                         tfbm_ltr3_results = NULL,
                                         rbp_motif_results = NULL,
                                         tfchip_ltr5_results = NULL,
                                         tfchip_ltr3_results = NULL,
                                         detailed_top_n = 15) {
  pieces <- list(
    .binary_results_to_key(binary_results),
    .numeric_results_to_key(numeric_results),
    .binary_results_to_key(
      subfamily_results,
      layer = "subfamily",
      test_type = "category_enrichment"
    ),
    .binary_results_to_key(
      domain_type_results,
      layer = "internal_domains",
      test_type = "category_enrichment"
    ),
    .binary_results_to_key(
      tfchip_ltr5_results,
      layer = "tfchip_tf",
      test_type = "motif_enrichment",
      max_rows = detailed_top_n,
      label_prefix = "5\u2032 LTR"
    ),
    .binary_results_to_key(
      tfchip_ltr3_results,
      layer = "tfchip_tf",
      test_type = "motif_enrichment",
      max_rows = detailed_top_n,
      label_prefix = "3\u2032 LTR"
    ),
    .binary_results_to_key(
      tfbm_ltr5_results,
      layer = "tfbm",
      test_type = "motif_enrichment",
      max_rows = detailed_top_n,
      label_prefix = "5\u2032 LTR"
    ),
    .binary_results_to_key(
      tfbm_ltr3_results,
      layer = "tfbm",
      test_type = "motif_enrichment",
      max_rows = detailed_top_n,
      label_prefix = "3\u2032 LTR"
    ),
    .binary_results_to_key(
      rbp_motif_results,
      layer = "rbp",
      test_type = "motif_enrichment",
      max_rows = detailed_top_n
    )
  )

  pieces <- pieces[vapply(pieces, function(z) is.data.frame(z) && nrow(z) > 0, logical(1))]

  if (length(pieces) == 0) {
    return(data.frame(
      rank = integer(0),
      layer = character(0),
      feature = character(0),
      feature_label = character(0),
      test_type = character(0),
      effect_type = character(0),
      foreground_value = numeric(0),
      background_value = numeric(0),
      effect_size = numeric(0),
      p_value = numeric(0),
      padj = numeric(0),
      direction = character(0),
      significant = logical(0),
      stringsAsFactors = FALSE
    ))
  }

  out <- do.call(rbind, pieces)
  out <- unique(out)
  out <- .order_result_rows(out)
  out$rank <- seq_len(nrow(out))
  out <- out[, c(
    "rank", "layer", "feature", "feature_label", "test_type",
    "effect_type", "foreground_value", "background_value",
    "effect_size", "p_value", "padj", "direction", "significant"
  ), drop = FALSE]
  rownames(out) <- NULL
  out
}
