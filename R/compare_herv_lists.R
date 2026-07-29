#' Compare two HERV lists
#'
#' Compares a foreground HERV list against a background HERV universe using
#' HERVariumR annotation layers. Foreground IDs are removed from the supplied
#' background internally, so the statistical comparison uses mutually exclusive
#' groups. The workflow includes core HERV, transcript-context, internal-domain,
#' 5′ and 3′ LTR, motif-burden, and bundled external-regulatory evidence.
#' Detailed TFBM and RBP motif enrichment can optionally be added from
#' user-supplied files.
#'
#' @param foreground_ids Character vector of foreground HERV IDs.
#' @param background_ids Character vector defining the background universe.
#'   Foreground IDs are removed internally.
#' @param foreground_name Character label for the foreground group.
#' @param background_name Character label for the effective non-foreground
#'   background group.
#' @param annotation_file Optional path to the main transcript-context
#'   annotation file. If `NULL`, the bundled HERVariumR resource is used.
#' @param ifn_stat1_file Optional path to the STAT1 LTR summary file. If `NULL`,
#'   the bundled HERVariumR resource is used.
#' @param ifn_stat1stat2_irf_file Optional path to the STAT1/STAT2/IRF LTR
#'   summary file. If `NULL`, the bundled HERVariumR resource is used.
#' @param last_exon_file Optional path to the terminal-exon domain annotation
#'   file. If `NULL`, the bundled HERVariumR resource is used.
#' @param add_external_evidence Logical. Whether to add compact ENCODE/UCSC
#'   44-organ evidence and LTR-specific cCRE context.
#' @param internal_encode_file Optional path to
#'   `internal_encode_44organs_compact.tsv.gz`. If `NULL`, the bundled resource
#'   is used.
#' @param ltr_encode_file Optional path to
#'   `ltr_encode_44organs_compact.tsv.gz`. If `NULL`, the bundled resource is
#'   used.
#' @param add_tfchip Logical. Whether to add compact LTR TF ChIP-seq
#'   peak-cluster overlap evidence.
#' @param tfchip_file Optional path to `ltr_tfchip_compact.tsv.gz`. If `NULL`,
#'   the bundled resource is used.
#' @param domain_coverage_cutoff Numeric. Minimum HMM profile coverage required
#'   for an internal retroviral domain to be retained.
#' @param keep_raw_domain_columns Logical. Whether raw, unfiltered domain
#'   columns are retained internally during annotation and comparison.
#' @param output_dir Directory where compact outputs, objects, dashboards, and
#'   optional detailed results are written.
#' @param add_tfbm_details Logical. Whether to add detailed LTR TFBM motif
#'   enrichment from a user-supplied FIMO file.
#' @param fimo_file Path to the detailed LTR TFBM FIMO file. Required when
#'   `add_tfbm_details = TRUE`; this large file is not bundled.
#' @param qvalue_cutoff Numeric. Maximum q-value retained for detailed TFBM hits.
#' @param add_rbp_details Logical. Whether to add detailed RBP motif enrichment
#'   from a user-supplied motif file.
#' @param rbp_file Path to the detailed RBP motif file. Required when
#'   `add_rbp_details = TRUE`; this large file is not bundled.
#' @param rbp_qvalue_cutoff Numeric. Maximum q-value retained for detailed RBP
#'   hits.
#' @param use_awk Logical. Whether to use `awk` for memory-efficient filtering
#'   of large motif files.
#' @param min_feature_count Integer. Minimum total foreground-plus-background
#'   count required to retain a categorical enrichment feature.
#' @param output_level Output verbosity: `"standard"` or `"extended"`.
#'   Standard mode writes compact primary outputs and the unified key-results
#'   table. Extended mode additionally writes detailed tables and supports
#'   static plots.
#' @param make_dashboard Logical. Whether to generate the concise interactive
#'   HTML dashboard.
#' @param make_plots Logical or `NULL`. Whether to generate static comparison
#'   plots. `NULL` defaults to `FALSE` in standard mode and `TRUE` in extended
#'   mode.
#' @param plot_top_n Integer. Number of top features displayed in plots and
#'   dashboard sections.
#' @param plot_all_binary_features Logical. Whether the static binary-feature
#'   plot may include all tested binary features.
#' @param plot_padj_cutoff Numeric. FDR cutoff used to classify plotted results
#'   as significant.
#' @param plot_show_only_significant Logical. Whether static plots show only
#'   FDR-significant results.
#' @param plot_pseudocount Numeric. Pseudocount used only when plotting log2
#'   odds ratios with zero contingency-table cells.
#' @param verbose Logical. Whether to print progress messages.
#' @param show_progress Logical. Whether to display progress bars for
#'   long-running motif-file operations.
#'
#' @return A `HERVarium_comparison` object containing mutually exclusive input
#'   groups, compact HERV-level features, the unified key-results table,
#'   underlying comparison tables, optional detailed regulatory enrichments,
#'   annotation metadata, and generated-file paths.
#' @export
compare_herv_lists <- function(foreground_ids,
                               background_ids,
                               foreground_name = "foreground",
                               background_name = "background",
                               annotation_file = NULL,
                               ifn_stat1_file = NULL,
                               ifn_stat1stat2_irf_file = NULL,
                               last_exon_file = NULL,

                               # External evidence layers
                               add_external_evidence = TRUE,
                               internal_encode_file = NULL,
                               ltr_encode_file = NULL,
                               add_tfchip = TRUE,
                               tfchip_file = NULL,

                               # Domain filtering
                               domain_coverage_cutoff = 0.40,
                               keep_raw_domain_columns = TRUE,

                               output_dir = "HERVariumR_comparison",

                               # Optional detailed LTR TFBM layer
                               add_tfbm_details = FALSE,
                               fimo_file = NULL,
                               qvalue_cutoff = 1,

                               # Optional detailed RBP layer
                               add_rbp_details = FALSE,
                               rbp_file = NULL,
                               rbp_qvalue_cutoff = 1,

                               use_awk = TRUE,

                               # Statistical options
                               min_feature_count = 1,

                               # Output/report options
                               output_level = c("standard", "extended"),
                               make_dashboard = TRUE,

                               # Plotting options
                               make_plots = NULL,
                               plot_top_n = 25,
                               plot_all_binary_features = TRUE,
                               plot_padj_cutoff = 0.05,
                               plot_show_only_significant = FALSE,
                               plot_pseudocount = 0.5,

                               # Verbosity
                               verbose = TRUE,
                               show_progress = TRUE) {

  # ------------------------------------------------------------
  # 0. Checks and output layout
  # ------------------------------------------------------------

  output_level <- .match_output_level(output_level)
  if (is.null(make_plots)) {
    make_plots <- identical(output_level, "extended")
  }

  layout <- .prepare_output_layout(output_dir, output_level)
  if (isTRUE(make_plots) && !dir.exists(layout$plots_dir)) {
    dir.create(layout$plots_dir, recursive = TRUE)
  }

  transient_dirs <- character(0)
  on.exit({
    if (length(transient_dirs) > 0) {
      unlink(transient_dirs, recursive = TRUE, force = TRUE)
    }
  }, add = TRUE)

  if (missing(foreground_ids) || length(foreground_ids) == 0) {
    stop("Please provide foreground_ids.")
  }

  if (missing(background_ids) || length(background_ids) == 0) {
    stop("Please provide background_ids.")
  }

  foreground_ids <- .clean_herv_id_vector(foreground_ids)
  n_invalid_foreground_ids <- attr(foreground_ids, "n_invalid_removed")
  n_duplicate_foreground_ids <- attr(foreground_ids, "n_duplicates_removed")

  background_ids <- .clean_herv_id_vector(background_ids)
  n_invalid_background_ids <- attr(background_ids, "n_invalid_removed")
  n_duplicate_background_ids <- attr(background_ids, "n_duplicates_removed")

  if (length(foreground_ids) == 0) {
    stop("No valid foreground IDs remain after removing NA, empty, '.', and 'NA' values.")
  }

  if (length(background_ids) == 0) {
    stop("No valid background IDs remain after removing NA, empty, '.', and 'NA' values.")
  }

  if (verbose && n_invalid_foreground_ids > 0) {
    message("Dropped invalid foreground IDs: ", n_invalid_foreground_ids, " (NA, empty, '.', or 'NA').")
  }

  if (verbose && n_invalid_background_ids > 0) {
    message("Dropped invalid background IDs: ", n_invalid_background_ids, " (NA, empty, '.', or 'NA').")
  }

  if (verbose && n_duplicate_foreground_ids > 0) {
    message("Dropped duplicated foreground IDs: ", n_duplicate_foreground_ids)
  }

  if (verbose && n_duplicate_background_ids > 0) {
    message("Dropped duplicated background IDs: ", n_duplicate_background_ids)
  }

  # Important: remove foreground from background to avoid duplicated group assignment
  background_ids <- setdiff(background_ids, foreground_ids)

  if (length(background_ids) == 0) {
    stop("After removing foreground IDs, background_ids is empty.")
  }

  all_ids <- unique(c(foreground_ids, background_ids))

  if (verbose) {
    message("Starting HERVariumR list comparison")
    message("Output level: ", output_level)
    message("Foreground: ", foreground_name, " (", length(foreground_ids), " IDs)")
    message("Background: ", background_name, " (", length(background_ids), " IDs)")
    message("Total unique IDs to annotate: ", length(all_ids))
    message("Output directory: ", output_dir)
  }

  # ------------------------------------------------------------
  # 1. Annotate combined universe
  # ------------------------------------------------------------

  anno <- annotate_hervs(
    herv_ids = all_ids,
    annotation_file = annotation_file,
    ifn_stat1_file = ifn_stat1_file,
    ifn_stat1stat2_irf_file = ifn_stat1stat2_irf_file,
    last_exon_file = last_exon_file,
    domain_coverage_cutoff = domain_coverage_cutoff,
    keep_raw_domain_columns = keep_raw_domain_columns,
    output_dir = NULL,
    verbose = verbose,
    .return_full_features = TRUE
  )

  if (add_external_evidence) {
    if (verbose) {
      message("Adding HERVarium external-evidence layers for comparison.")
    }

    anno <- add_encode_44organ_layers(
      anno,
      internal_encode_file = internal_encode_file,
      ltr_encode_file = ltr_encode_file,
      output_dir = NULL,
      verbose = verbose,
      .return_full_features = TRUE
    )
  }

  if (add_tfchip) {
    if (verbose) {
      message("Adding LTR TF ChIP-seq evidence layer for comparison.")
    }

    anno <- add_ltr_tfchip_layer(
      anno,
      tfchip_file = tfchip_file,
      output_dir = NULL,
      verbose = verbose,
      .return_full_features = TRUE
    )
  }

  features <- anno$features

  features$comparison_group <- ifelse(
    features$HERV_id %in% foreground_ids,
    foreground_name,
    background_name
  )

  # Keep only foreground/background after matching
  features <- features[
    features$comparison_group %in% c(foreground_name, background_name),
  ]

  # ------------------------------------------------------------
  # 2. Write input summary
  # ------------------------------------------------------------

  input_summary <- data.frame(
    group = c(foreground_name, background_name),
    n_input_ids = c(length(foreground_ids), length(background_ids)),
    n_matched_ids = c(
      sum(features$comparison_group == foreground_name),
      sum(features$comparison_group == background_name)
    ),
    stringsAsFactors = FALSE
  )

  .write_hervarium_tsv(
    input_summary,
    file.path(output_dir, "01_input_summary.tsv")
  )

  .write_hervarium_tsv(
    select_compact_herv_features(features, include_comparison_group = TRUE),
    file.path(output_dir, "02_features_with_comparison_group.tsv")
  )

  # ------------------------------------------------------------
  # Helper functions
  # ------------------------------------------------------------

  safe_p_adjust <- function(p) {
    out <- rep(NA_real_, length(p))
    ok <- !is.na(p)
    out[ok] <- p.adjust(p[ok], method = "BH")
    out
  }

  fisher_feature <- function(feature_name, yes_vector, group_vector) {

    yes_vector <- as.logical(yes_vector)
    yes_vector[is.na(yes_vector)] <- FALSE

    fg <- group_vector == foreground_name
    bg <- group_vector == background_name

    a <- sum(yes_vector[fg], na.rm = TRUE)
    b <- sum(!yes_vector[fg], na.rm = TRUE)
    c <- sum(yes_vector[bg], na.rm = TRUE)
    d <- sum(!yes_vector[bg], na.rm = TRUE)

    mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)

    test <- tryCatch(
      fisher.test(mat),
      error = function(e) NULL
    )

    if (is.null(test)) {
      odds_ratio <- NA_real_
      p_value <- NA_real_
    } else {
      odds_ratio <- unname(test$estimate)
      p_value <- test$p.value
    }

    fg_percent <- ifelse((a + b) > 0, 100 * a / (a + b), NA_real_)
    bg_percent <- ifelse((c + d) > 0, 100 * c / (c + d), NA_real_)

    direction <- ifelse(
      is.na(fg_percent) | is.na(bg_percent),
      NA_character_,
      ifelse(fg_percent > bg_percent, "enriched_in_foreground",
             ifelse(fg_percent < bg_percent, "depleted_in_foreground", "no_difference"))
    )

    data.frame(
      feature = feature_name,
      foreground_yes = a,
      foreground_no = b,
      foreground_total = a + b,
      foreground_percent = fg_percent,
      background_yes = c,
      background_no = d,
      background_total = c + d,
      background_percent = bg_percent,
      odds_ratio = odds_ratio,
      p_value = p_value,
      direction = direction,
      stringsAsFactors = FALSE
    )
  }

  wilcox_feature <- function(feature_name, value_vector, group_vector) {

    value_vector <- suppressWarnings(as.numeric(value_vector))

    fg_values <- value_vector[group_vector == foreground_name]
    bg_values <- value_vector[group_vector == background_name]

    fg_values <- fg_values[!is.na(fg_values)]
    bg_values <- bg_values[!is.na(bg_values)]

    if (length(fg_values) == 0 || length(bg_values) == 0) {
      p_value <- NA_real_
    } else {
      p_value <- tryCatch(
        wilcox.test(fg_values, bg_values, exact = FALSE)$p.value,
        error = function(e) NA_real_
      )
    }

    med_fg <- ifelse(length(fg_values) > 0, median(fg_values), NA_real_)
    med_bg <- ifelse(length(bg_values) > 0, median(bg_values), NA_real_)
    mean_fg <- ifelse(length(fg_values) > 0, mean(fg_values), NA_real_)
    mean_bg <- ifelse(length(bg_values) > 0, mean(bg_values), NA_real_)

    direction <- ifelse(
      is.na(med_fg) | is.na(med_bg),
      NA_character_,
      ifelse(
        med_fg > med_bg,
        "higher_in_foreground",
        ifelse(
          med_fg < med_bg,
          "lower_in_foreground",
          ifelse(
            is.na(mean_fg) | is.na(mean_bg),
            "no_difference",
            ifelse(
              mean_fg > mean_bg,
              "higher_in_foreground_mean_only",
              ifelse(
                mean_fg < mean_bg,
                "lower_in_foreground_mean_only",
                "no_difference"
              )
            )
          )
        )
      )
    )

    data.frame(
      feature = feature_name,
      foreground_n = length(fg_values),
      background_n = length(bg_values),
      foreground_median = med_fg,
      background_median = med_bg,
      delta_median = med_fg - med_bg,
      foreground_mean = mean_fg,
      background_mean = mean_bg,
      delta_mean = mean_fg - mean_bg,
      p_value = p_value,
      direction = direction,
      stringsAsFactors = FALSE
    )
  }

  parse_domains_type <- function(z, HERV_id, group) {

    if (is.na(z) || z == "" || z == ".") {
      return(NULL)
    }

    entries <- unlist(strsplit(as.character(z), ";", fixed = TRUE))
    entries <- entries[entries != "" & entries != "."]

    if (length(entries) == 0) {
      return(NULL)
    }

    out <- lapply(entries, function(e) {

      parts1 <- unlist(strsplit(e, "|", fixed = TRUE))
      if (length(parts1) != 2) return(NULL)

      gene_class <- parts1[1]
      domain_and_cov <- parts1[2]

      parts2 <- unlist(strsplit(domain_and_cov, ":", fixed = TRUE))
      if (length(parts2) != 2) return(NULL)

      domain_type <- parts2[1]
      coverage <- suppressWarnings(as.numeric(parts2[2]))

      if (is.na(coverage)) return(NULL)

      data.frame(
        HERV_id = HERV_id,
        comparison_group = group,
        gene_class = gene_class,
        domain_type = domain_type,
        coverage = coverage,
        raw_domain_entry = e,
        stringsAsFactors = FALSE
      )
    })

    do.call(rbind, out)
  }

  clean_domain_type <- function(z) {
    z <- as.character(z)
    z[z == "ENV"] <- "Env"
    z[z == "GAG"] <- "Gag"
    z[z == "AP"] <- "Protease"
    z[z == "DUT"] <- "dUTPase"
    z[z == "INT"] <- "Integrase"
    z[z == "RT"] <- "RT"
    z[z == "RH"] <- "RNaseH"
    z
  }

  # ------------------------------------------------------------
  # 3. Binary feature enrichment
  # ------------------------------------------------------------

  binary_features <- c(
    "has_domain",
    "has_gag",
    "has_pol",
    "has_env",
    "has_accessory",
    "has_complete_gag_pol_env",
    "has_ltr5",
    "has_ltr3",
    "has_both_ltrs",
    "has_ltr5_stat1_motif",
    "has_ltr3_stat1_motif",
    "has_ltr5_stat1stat2_motif",
    "has_ltr3_stat1stat2_motif",
    "has_any_terminal_exon_domain",
    "has_terminal_exon_domain_protein_coding",
    "has_terminal_exon_domain_lncRNA",
    "is_intergenic",
    "overlaps_gene",
    "overlaps_exon",
    "overlaps_intron",
    "overlaps_cds",
    "overlaps_utr",
    "overlaps_lncRNA",
    "overlaps_protein_coding",

    # Compact internal ENCODE/UCSC evidence.
    # Internal-region cCRE columns are intentionally not tested by default.
    "has_internal_encode_evidence",
    "has_internal_multilayer_encode_evidence",
    "has_internal_robust_multilayer_encode_evidence",

    # Compact LTR-specific ENCODE/UCSC and cCRE evidence.
    "ltr5_encode_evidence_including_ccre",
    "ltr3_encode_evidence_including_ccre",
    "ltr5_encode_ccre_overlap",
    "ltr3_encode_ccre_overlap",

    # Compact LTR-specific TF ChIP-seq peak-cluster evidence.
    "has_ltr5_tfchip_overlap",
    "has_ltr3_tfchip_overlap"
  )

  binary_features <- binary_features[binary_features %in% colnames(features)]

  binary_results <- do.call(
    rbind,
    lapply(binary_features, function(f) {
      fisher_feature(f, features[[f]], features$comparison_group)
    })
  )

  if (!is.null(binary_results) && nrow(binary_results) > 0) {
    binary_results$padj <- safe_p_adjust(binary_results$p_value)
    binary_results <- binary_results[order(binary_results$padj, binary_results$p_value), ]
  }

  .write_detail_tsv(binary_results, "03_comparison_binary_features.tsv", layout)

  external_evidence_feature_pattern <- paste(
    c(
      "encode",
      "ccre",
      "dnase",
      "h3k27ac",
      "h3k4me3",
      "transcription",
      "tfchip"
    ),
    collapse = "|"
  )

  external_binary_results <- binary_results[
    grepl(external_evidence_feature_pattern, binary_results$feature),
  ]

  .write_detail_tsv(external_binary_results, "16_comparison_external_evidence_binary_features.tsv", layout)




  # ------------------------------------------------------------
  # 4. Numeric feature comparison
  # ------------------------------------------------------------

  numeric_features <- c(
    "domain_count",
    "max_domain_cov",
    "ltr5_tfbm_burden",
    "ltr3_tfbm_burden",
    "rbp_burden",
    "rbp_unique",
    "max_terminal_domain_coverage_protein_coding",
    "max_terminal_domain_coverage_lncRNA",

    # Compact internal ENCODE/UCSC detected-layer count.
    "internal_encode_detected_layer_count",

    # Compact LTR-specific ENCODE/UCSC detected-layer and cCRE counts.
    "ltr5_encode_detected_layer_count",
    "ltr3_encode_detected_layer_count",
    "ltr5_encode_ccre_n_overlaps",
    "ltr3_encode_ccre_n_overlaps",

    # Compact LTR-specific TF ChIP-seq burden.
    "ltr5_tfchip_n_peak_clusters",
    "ltr3_tfchip_n_peak_clusters",
    "ltr5_tfchip_n_tfs",
    "ltr3_tfchip_n_tfs"
  )

  numeric_features <- numeric_features[numeric_features %in% colnames(features)]

  numeric_results <- do.call(
    rbind,
    lapply(numeric_features, function(f) {
      wilcox_feature(f, features[[f]], features$comparison_group)
    })
  )

  if (!is.null(numeric_results) && nrow(numeric_results) > 0) {
    numeric_results$padj <- safe_p_adjust(numeric_results$p_value)
    numeric_results <- numeric_results[order(numeric_results$padj, numeric_results$p_value), ]
  }

  .write_detail_tsv(numeric_results, "04_comparison_numeric_features.tsv", layout)

  external_numeric_results <- numeric_results[
    grepl(external_evidence_feature_pattern, numeric_results$feature),
  ]

  .write_detail_tsv(external_numeric_results, "17_comparison_external_evidence_numeric_features.tsv", layout)


  # ------------------------------------------------------------
  # 5. Subfamily enrichment
  # ------------------------------------------------------------

  subfamilies <- sort(unique(features$subfamily))
  subfamilies <- subfamilies[!is.na(subfamilies) & subfamilies != "" & subfamilies != "."]

  subfamily_results <- do.call(
    rbind,
    lapply(subfamilies, function(sf) {
      fisher_feature(
        feature_name = sf,
        yes_vector = features$subfamily == sf,
        group_vector = features$comparison_group
      )
    })
  )

  if (!is.null(subfamily_results) && nrow(subfamily_results) > 0) {
    subfamily_results <- subfamily_results[
      subfamily_results$foreground_yes + subfamily_results$background_yes >= min_feature_count,
    ]
    subfamily_results$padj <- safe_p_adjust(subfamily_results$p_value)
    subfamily_results <- subfamily_results[order(subfamily_results$padj, subfamily_results$p_value), ]
  }

  .write_detail_tsv(subfamily_results, "05_comparison_subfamily_enrichment.tsv", layout)

  # ------------------------------------------------------------
  # 6. Domain type enrichment
  # ------------------------------------------------------------

  domain_hits <- NULL

  if ("domains_type" %in% colnames(features)) {
    domain_hits <- do.call(
      rbind,
      lapply(seq_len(nrow(features)), function(i) {
        parse_domains_type(
          z = features$domains_type[i],
          HERV_id = features$HERV_id[i],
          group = features$comparison_group[i]
        )
      })
    )
  }

  domain_type_results <- data.frame()

  if (!is.null(domain_hits) && nrow(domain_hits) > 0) {

    domain_hits$domain_type <- clean_domain_type(domain_hits$domain_type)

    domain_presence <- unique(domain_hits[, c("HERV_id", "comparison_group", "domain_type")])

    domain_types <- sort(unique(domain_presence$domain_type))

    domain_type_results <- do.call(
      rbind,
      lapply(domain_types, function(dt) {
        yes_ids <- unique(domain_presence$HERV_id[domain_presence$domain_type == dt])

        fisher_feature(
          feature_name = dt,
          yes_vector = features$HERV_id %in% yes_ids,
          group_vector = features$comparison_group
        )
      })
    )

    domain_type_results <- domain_type_results[
      domain_type_results$foreground_yes + domain_type_results$background_yes >= min_feature_count,
    ]

    domain_type_results$padj <- safe_p_adjust(domain_type_results$p_value)
    domain_type_results <- domain_type_results[
      order(domain_type_results$padj, domain_type_results$p_value),
    ]

    .write_detail_tsv(domain_hits, "06_domain_hits_long.tsv", layout)
  }

  .write_detail_tsv(domain_type_results, "07_comparison_domain_type_enrichment.tsv", layout)

  # ------------------------------------------------------------
  # 7. IFN/LTR curated feature enrichment
  # ------------------------------------------------------------

  ifn_features <- c(
    "has_ltr5_stat1_motif",
    "has_ltr3_stat1_motif",
    "has_ltr5_stat1stat2_motif",
    "has_ltr3_stat1stat2_motif"
  )

  ifn_features <- ifn_features[ifn_features %in% colnames(features)]

  ifn_results <- do.call(
    rbind,
    lapply(ifn_features, function(f) {
      fisher_feature(f, features[[f]], features$comparison_group)
    })
  )

  if (!is.null(ifn_results) && nrow(ifn_results) > 0) {
    ifn_results$padj <- safe_p_adjust(ifn_results$p_value)
    ifn_results <- ifn_results[order(ifn_results$padj, ifn_results$p_value), ]
  }

  .write_detail_tsv(ifn_results, "08_comparison_ifn_ltr_features.tsv", layout)

  # ------------------------------------------------------------
  # 8. LTR TFBM burden dedicated table
  # ------------------------------------------------------------

  ltr_tfbm_burden_features <- c(
    "ltr5_tfbm_burden",
    "ltr3_tfbm_burden"
  )

  ltr_tfbm_burden_features <- ltr_tfbm_burden_features[
    ltr_tfbm_burden_features %in% colnames(features)
  ]

  ltr_tfbm_burden_results <- do.call(
    rbind,
    lapply(ltr_tfbm_burden_features, function(f) {
      wilcox_feature(f, features[[f]], features$comparison_group)
    })
  )

  if (!is.null(ltr_tfbm_burden_results) && nrow(ltr_tfbm_burden_results) > 0) {
    ltr_tfbm_burden_results$padj <- safe_p_adjust(ltr_tfbm_burden_results$p_value)
    ltr_tfbm_burden_results <- ltr_tfbm_burden_results[
      order(ltr_tfbm_burden_results$padj, ltr_tfbm_burden_results$p_value),
    ]
  }

  .write_detail_tsv(ltr_tfbm_burden_results, "09_comparison_ltr_tfbm_burden.tsv", layout)

  # ------------------------------------------------------------
  # 9. RBP burden dedicated table
  # ------------------------------------------------------------

  rbp_burden_features <- c(
    "rbp_burden",
    "rbp_unique"
  )

  rbp_burden_features <- rbp_burden_features[
    rbp_burden_features %in% colnames(features)
  ]

  rbp_burden_results <- do.call(
    rbind,
    lapply(rbp_burden_features, function(f) {
      wilcox_feature(f, features[[f]], features$comparison_group)
    })
  )

  if (!is.null(rbp_burden_results) && nrow(rbp_burden_results) > 0) {
    rbp_burden_results$padj <- safe_p_adjust(rbp_burden_results$p_value)
    rbp_burden_results <- rbp_burden_results[
      order(rbp_burden_results$padj, rbp_burden_results$p_value),
    ]
  }

  .write_detail_tsv(rbp_burden_results, "10_comparison_rbp_burden.tsv", layout)

  # ------------------------------------------------------------
  # 9b. Optional detailed RBP layer
  # ------------------------------------------------------------

  if (add_rbp_details) {

    if (is.null(rbp_file)) {
      stop("add_rbp_details = TRUE requires rbp_file.")
    }

    if (verbose) {
      message("Adding detailed RBP layer for comparison.")
    }

    anno$features <- features

    rbp_layer_dir <- .make_transient_layer_dir(layout, "detailed_rbp")
    if (!isTRUE(layout$write_details)) transient_dirs <- c(transient_dirs, rbp_layer_dir)

    anno <- add_rbp_details(
      anno,
      rbp_file = rbp_file,
      output_dir = rbp_layer_dir,
      qvalue_cutoff = rbp_qvalue_cutoff,
      use_awk = use_awk,
      show_progress = show_progress,
      .return_full_features = TRUE
    )

    features <- anno$features

    # Restore comparison group after feature table update
    features$comparison_group <- ifelse(
      features$HERV_id %in% foreground_ids,
      foreground_name,
      background_name
    )

  }


  # ------------------------------------------------------------
  # 10. Optional detailed LTR TFBM motif enrichment
  # ------------------------------------------------------------

  tfbm_ltr5_results <- data.frame()
  tfbm_ltr3_results <- data.frame()
  tfbm_any_results <- data.frame()

  if (add_tfbm_details) {

    if (is.null(fimo_file)) {
      stop("add_tfbm_details = TRUE requires fimo_file.")
    }

    if (verbose) {
      message("Adding detailed TFBM layer for comparison.")
    }

    anno$features <- features

    tfbm_layer_dir <- .make_transient_layer_dir(layout, "detailed_tfbm")
    if (!isTRUE(layout$write_details)) transient_dirs <- c(transient_dirs, tfbm_layer_dir)

    anno <- add_ltr_tfbm_details(
      anno,
      fimo_file = fimo_file,
      output_dir = tfbm_layer_dir,
      qvalue_cutoff = qvalue_cutoff,
      use_awk = use_awk,
      show_progress = show_progress,
      .return_full_features = TRUE
    )

    features <- anno$features
    features$comparison_group <- ifelse(
      features$HERV_id %in% foreground_ids,
      foreground_name,
      background_name
    )

    tfbm_hits <- anno$ltr_tfbm_hits

    if (!is.null(tfbm_hits) && nrow(tfbm_hits) > 0) {

      # Remove possible duplicated subfamily column from FIMO file
      if ("subfamily" %in% colnames(tfbm_hits)) {
        tfbm_hits$subfamily <- NULL
      }

      group_info <- unique(features[, c("HERV_id", "comparison_group")])

      tfbm_hits <- merge(
        tfbm_hits,
        group_info,
        by = "HERV_id",
        all.x = TRUE
      )

      motif_enrichment <- function(tfbm_hits_subset, label) {

        motif_presence <- unique(
          tfbm_hits_subset[, c("HERV_id", "motif_alt_id")]
        )

        motifs <- sort(unique(motif_presence$motif_alt_id))
        motifs <- motifs[!is.na(motifs) & motifs != "" & motifs != "."]

        if (length(motifs) == 0) {
          return(data.frame())
        }

        out <- do.call(
          rbind,
          lapply(motifs, function(motif) {
            yes_ids <- unique(motif_presence$HERV_id[
              motif_presence$motif_alt_id == motif
            ])

            fisher_feature(
              feature_name = motif,
              yes_vector = features$HERV_id %in% yes_ids,
              group_vector = features$comparison_group
            )
          })
        )

        out <- out[
          out$foreground_yes + out$background_yes >= min_feature_count,
        ]

        out$comparison_layer <- label
        out$padj <- safe_p_adjust(out$p_value)
        out <- out[order(out$padj, out$p_value), ]
        out
      }

      tfbm_ltr5_results <- motif_enrichment(
        tfbm_hits[tfbm_hits$ltr_position == "LTR5", ],
        "LTR5"
      )

      tfbm_ltr3_results <- motif_enrichment(
        tfbm_hits[tfbm_hits$ltr_position == "LTR3", ],
        "LTR3"
      )

      tfbm_any_results <- motif_enrichment(
        tfbm_hits,
        "Any_LTR"
      )
    }
  }

  .write_detail_tsv(tfbm_ltr5_results, "11_comparison_ltr5_tfbm_motif_enrichment.tsv", layout)

  .write_detail_tsv(tfbm_ltr3_results, "12_comparison_ltr3_tfbm_motif_enrichment.tsv", layout)

  .write_detail_tsv(tfbm_any_results, "13_comparison_any_ltr_tfbm_motif_enrichment.tsv", layout)

  # ------------------------------------------------------------
  # 11. Optional RBP motif enrichment if motif-name columns exist
  # ------------------------------------------------------------

  rbp_motif_results <- data.frame()

  rbp_candidate_cols <- c(
    "rbp_names_all",
    "top_rbp_names_preview",
    "rbp_names",
    "rbp_motifs",
    "rbp_motif_names",
    "rbp_proteins",
    "rbp_tf_names_all"
  )

  rbp_col <- rbp_candidate_cols[rbp_candidate_cols %in% colnames(features)]
  rbp_col <- ifelse(length(rbp_col) > 0, rbp_col[1], NA_character_)

  if (!is.na(rbp_col)) {

    rbp_long <- do.call(
      rbind,
      lapply(seq_len(nrow(features)), function(i) {

        z <- features[[rbp_col]][i]

        if (is.na(z) || z == "" || z == ".") {
          return(NULL)
        }

        rbps <- unlist(strsplit(as.character(z), ";", fixed = TRUE))
        rbps <- rbps[rbps != "" & rbps != "."]

        if (length(rbps) == 0) {
          return(NULL)
        }

        data.frame(
          HERV_id = features$HERV_id[i],
          comparison_group = features$comparison_group[i],
          rbp_motif = unique(rbps),
          stringsAsFactors = FALSE
        )
      })
    )

    if (!is.null(rbp_long) && nrow(rbp_long) > 0) {

      rbps <- sort(unique(rbp_long$rbp_motif))

      rbp_motif_results <- do.call(
        rbind,
        lapply(rbps, function(rbp) {
          yes_ids <- unique(rbp_long$HERV_id[rbp_long$rbp_motif == rbp])

          fisher_feature(
            feature_name = rbp,
            yes_vector = features$HERV_id %in% yes_ids,
            group_vector = features$comparison_group
          )
        })
      )

      rbp_motif_results <- rbp_motif_results[
        rbp_motif_results$foreground_yes + rbp_motif_results$background_yes >= min_feature_count,
      ]

      rbp_motif_results$padj <- safe_p_adjust(rbp_motif_results$p_value)
      rbp_motif_results <- rbp_motif_results[
        order(rbp_motif_results$padj, rbp_motif_results$p_value),
      ]

      .write_detail_tsv(rbp_long, "14_rbp_motifs_long.tsv", layout)
    }
  } else {
    if (verbose) {
      message("No per-HERV RBP motif-name column found. Only RBP burden was compared.")
    }
  }

  .write_detail_tsv(rbp_motif_results, "15_comparison_rbp_motif_enrichment.tsv", layout)

  # ------------------------------------------------------------
  # 12. LTR TF ChIP-seq TF-name enrichment
  # ------------------------------------------------------------

  tfchip_burden_features <- c(
    "ltr5_tfchip_n_peak_clusters",
    "ltr3_tfchip_n_peak_clusters",
    "ltr5_tfchip_n_tfs",
    "ltr3_tfchip_n_tfs"
  )

  tfchip_burden_features <- tfchip_burden_features[
    tfchip_burden_features %in% colnames(features)
  ]

  tfchip_burden_results <- data.frame()

  if (length(tfchip_burden_features) > 0) {
    tfchip_burden_results <- do.call(
      rbind,
      lapply(tfchip_burden_features, function(f) {
        wilcox_feature(f, features[[f]], features$comparison_group)
      })
    )

    if (!is.null(tfchip_burden_results) && nrow(tfchip_burden_results) > 0) {
      tfchip_burden_results$padj <- safe_p_adjust(tfchip_burden_results$p_value)
      tfchip_burden_results <- tfchip_burden_results[
        order(tfchip_burden_results$padj, tfchip_burden_results$p_value),
      ]
    }
  }

  .write_detail_tsv(tfchip_burden_results, "18_comparison_tfchip_burden.tsv", layout)

  split_tfchip_names <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & x != "" & x != "."]
    if (length(x) == 0) return(character(0))

    x <- paste(x, collapse = ",")
    z <- unlist(strsplit(x, "[,;]", perl = TRUE))
    z <- trimws(z)
    z <- z[!is.na(z) & z != "" & z != "."]
    unique(z)
  }

  tfchip_presence_from_cols <- function(cols) {
    cols <- cols[cols %in% colnames(features)]

    if (length(cols) == 0) {
      return(data.frame())
    }

    out <- do.call(
      rbind,
      lapply(seq_len(nrow(features)), function(i) {
        tf_names <- unique(unlist(lapply(cols, function(col) {
          split_tfchip_names(features[[col]][i])
        })))

        tf_names <- tf_names[!is.na(tf_names) & tf_names != "" & tf_names != "."]

        if (length(tf_names) == 0) {
          return(NULL)
        }

        data.frame(
          HERV_id = features$HERV_id[i],
          comparison_group = features$comparison_group[i],
          tf_name = tf_names,
          stringsAsFactors = FALSE
        )
      })
    )

    if (is.null(out)) {
      return(data.frame())
    }

    unique(out)
  }

  tfchip_tf_enrichment <- function(tfchip_presence, label) {
    if (is.null(tfchip_presence) || nrow(tfchip_presence) == 0) {
      return(data.frame())
    }

    tf_names <- sort(unique(tfchip_presence$tf_name))
    tf_names <- tf_names[!is.na(tf_names) & tf_names != "" & tf_names != "."]

    if (length(tf_names) == 0) {
      return(data.frame())
    }

    out <- do.call(
      rbind,
      lapply(tf_names, function(tf_name) {
        yes_ids <- unique(tfchip_presence$HERV_id[
          tfchip_presence$tf_name == tf_name
        ])

        fisher_feature(
          feature_name = tf_name,
          yes_vector = features$HERV_id %in% yes_ids,
          group_vector = features$comparison_group
        )
      })
    )

    out <- out[
      out$foreground_yes + out$background_yes >= min_feature_count,
    ]

    out$comparison_layer <- label
    out$padj <- safe_p_adjust(out$p_value)
    out <- out[order(out$padj, out$p_value), ]
    out
  }

  tfchip_ltr5_presence <- tfchip_presence_from_cols(c(
    "ltr5_tfchip_tf_list",
    "ltr5_tfchip_tf_list_full"
  ))

  tfchip_ltr3_presence <- tfchip_presence_from_cols(c(
    "ltr3_tfchip_tf_list",
    "ltr3_tfchip_tf_list_full"
  ))

  tfchip_any_presence <- unique(rbind(tfchip_ltr5_presence, tfchip_ltr3_presence))

  tfchip_ltr5_results <- tfchip_tf_enrichment(tfchip_ltr5_presence, "LTR5")
  tfchip_ltr3_results <- tfchip_tf_enrichment(tfchip_ltr3_presence, "LTR3")
  tfchip_any_results <- tfchip_tf_enrichment(tfchip_any_presence, "Any_LTR")

  .write_detail_tsv(tfchip_ltr5_results, "19_comparison_ltr5_tfchip_tf_enrichment.tsv", layout)

  .write_detail_tsv(tfchip_ltr3_results, "20_comparison_ltr3_tfchip_tf_enrichment.tsv", layout)

  .write_detail_tsv(tfchip_any_results, "21_comparison_any_ltr_tfchip_tf_enrichment.tsv", layout)

  # ------------------------------------------------------------
  # 13. Unified key-results table and compact outputs
  # ------------------------------------------------------------

  key_results <- build_comparison_key_results(
    binary_results = binary_results,
    numeric_results = numeric_results,
    subfamily_results = subfamily_results,
    domain_type_results = domain_type_results,
    tfbm_ltr5_results = tfbm_ltr5_results,
    tfbm_ltr3_results = tfbm_ltr3_results,
    rbp_motif_results = rbp_motif_results,
    tfchip_ltr5_results = tfchip_ltr5_results,
    tfchip_ltr3_results = tfchip_ltr3_results,
    detailed_top_n = min(plot_top_n, 15)
  )

  key_results_file <- file.path(output_dir, "03_comparison_key_results.tsv")
  .write_hervarium_tsv(key_results, key_results_file)

  compact_features <- select_compact_herv_features(
    features,
    include_comparison_group = TRUE
  )

  features_file <- file.path(output_dir, "02_features_with_comparison_group.tsv")
  .write_hervarium_tsv(compact_features, features_file)

  anno$features <- select_compact_herv_features(features)
  anno$output_dir <- output_dir
  anno$output_level <- output_level

  missing_file <- file.path(output_dir, "missing_herv_ids.txt")
  if (length(anno$missing_ids) > 0) {
    writeLines(anno$missing_ids, con = missing_file)
  } else if (file.exists(missing_file)) {
    unlink(missing_file)
  }

  if (isTRUE(layout$write_details)) {
    matched_report <- anno$summary
    if ("domains_gene" %in% colnames(matched_report)) matched_report$domains_gene <- NULL
    .write_detail_tsv(matched_report, "matched_annotation_summary.tsv", layout)
    .write_detail_tsv(anno$stats, "annotation_stats.tsv", layout)
  }

  # ------------------------------------------------------------
  # 14. Return object
  # ------------------------------------------------------------

  object_file <- file.path(output_dir, "hervarium_comparison_object.rds")
  dashboard_file <- file.path(output_dir, "hervarium_comparison_dashboard.html")

  out <- list(
    foreground_name = foreground_name,
    background_name = background_name,
    foreground_ids = foreground_ids,
    background_ids = background_ids,
    annotation = anno,
    features = compact_features,
    input_summary = input_summary,
    key_results = key_results,
    binary_features = binary_results,
    numeric_features = numeric_results,
    external_evidence_binary_features = external_binary_results,
    external_evidence_numeric_features = external_numeric_results,
    subfamily_enrichment = subfamily_results,
    domain_hits = domain_hits,
    domain_type_enrichment = domain_type_results,
    domain_coverage_cutoff = domain_coverage_cutoff,
    ifn_ltr_features = ifn_results,
    ltr_tfbm_burden = ltr_tfbm_burden_results,
    rbp_burden = rbp_burden_results,
    tfchip_burden = tfchip_burden_results,
    ltr5_tfbm_motif_enrichment = tfbm_ltr5_results,
    ltr3_tfbm_motif_enrichment = tfbm_ltr3_results,
    any_ltr_tfbm_motif_enrichment = tfbm_any_results,
    rbp_motif_enrichment = rbp_motif_results,
    ltr5_tfchip_tf_enrichment = tfchip_ltr5_results,
    ltr3_tfchip_tf_enrichment = tfchip_ltr3_results,
    any_ltr_tfchip_tf_enrichment = tfchip_any_results,
    output_level = output_level,
    output_dir = output_dir,
    output_files = list(
      input_summary = file.path(output_dir, "01_input_summary.tsv"),
      compact_features = features_file,
      key_results = key_results_file,
      object = object_file,
      dashboard = if (isTRUE(make_dashboard)) dashboard_file else NULL,
      details_dir = if (isTRUE(layout$write_details)) layout$details_dir else NULL,
      plots_dir = if (isTRUE(make_plots)) layout$plots_dir else NULL
    )
  )

  class(out) <- "HERVarium_comparison"
  saveRDS(out, object_file)

  # ------------------------------------------------------------
  # 15. Optional extended static plots and concise dashboard
  # ------------------------------------------------------------

  if (isTRUE(make_plots)) {
    if (exists("plot_herv_comparison", mode = "function")) {
      if (verbose) message("Generating extended HERVariumR comparison plots.")

      plot_herv_comparison(
        out,
        output_dir = output_dir,
        top_n = plot_top_n,
        plot_all_binary_features = plot_all_binary_features,
        padj_cutoff = plot_padj_cutoff,
        show_only_significant = plot_show_only_significant,
        pseudocount = plot_pseudocount
      )
    } else {
      warning(
        "make_plots = TRUE, but plot_herv_comparison() is not available. ",
        "Source R/plot_herv_comparison.R or include it in the package."
      )
    }
  }

  if (isTRUE(make_dashboard)) {
    generate_hervarium_dashboard(
      out,
      output_dir = output_dir,
      report_file = "hervarium_comparison_dashboard.html",
      top_n = plot_top_n
    )
  }

  saveRDS(out, object_file)

  if (verbose) {
    message("HERVariumR comparison finished.")
    message("Key results: ", key_results_file)
    if (isTRUE(make_dashboard)) message("Dashboard: ", dashboard_file)
    message("Results written to: ", output_dir)
  }

  out
}


print.HERVarium_comparison <- function(x, ...) {
  cat("HERVarium comparison result\n")
  cat("---------------------------\n")
  cat("Foreground:", x$foreground_name, "\n")
  cat("Background:", x$background_name, "\n")
  cat("\nInput summary:\n")
  print(x$input_summary)

  if (!is.null(x$output_dir)) {
    cat("\nOutput directory:\n")
    cat(x$output_dir, "\n")
  }

  invisible(x)
}