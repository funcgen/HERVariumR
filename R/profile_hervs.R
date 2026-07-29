#' Profile a single list of HERV IDs
#'
#' Annotates and summarizes a single list of HERV internal-region IDs using
#' HERVariumR annotation layers. Standard output is intentionally compact:
#' a HERV-level feature table, an input summary, an RDS object, and a concise
#' interactive dashboard. Extended output additionally writes detailed tables
#' and optional static plots under `details/` and `plots/`.
#'
#' @param herv_ids Character vector of HERV IDs to annotate.
#' @param annotation_file Optional path to the main transcript-context annotation
#'   file. If `NULL`, the bundled HERVariumR resource is used.
#' @param output_dir Directory where output tables, objects, plots, and reports
#'   are written.
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
#'   columns are retained internally during annotation.
#' @param add_tfbm_details Logical. Whether to add detailed LTR TFBM annotations
#'   from a user-supplied FIMO file.
#' @param fimo_file Path to the detailed LTR TFBM FIMO file. Required when
#'   `add_tfbm_details = TRUE`; this large file is not bundled.
#' @param qvalue_cutoff Numeric. Maximum q-value retained for detailed TFBM hits.
#' @param add_rbp_details Logical. Whether to add detailed RBP motif annotations
#'   from a user-supplied FIMO-like file.
#' @param rbp_file Path to the detailed RBP motif file. Required when
#'   `add_rbp_details = TRUE`; this large file is not bundled.
#' @param rbp_qvalue_cutoff Numeric. Maximum q-value retained for detailed RBP
#'   hits.
#' @param use_awk Logical. Whether to use `awk` for memory-efficient filtering
#'   of large motif files.
#' @param output_level Output verbosity: `"standard"` or `"extended"`.
#'   Standard mode writes the compact primary outputs. Extended mode also writes
#'   detailed tables and supports static plots.
#' @param make_dashboard Logical. Whether to generate the concise interactive
#'   HTML dashboard.
#' @param make_plots Logical or `NULL`. Whether to generate static summary
#'   plots. `NULL` defaults to `FALSE` in standard mode and `TRUE` in extended
#'   mode.
#' @param make_ltr_tfbm_plots Logical or `NULL`. Whether to generate detailed
#'   LTR TFBM static plots when detailed TFBM results are available. `NULL`
#'   follows the selected output level.
#' @param make_internal_domain_plot Logical or `NULL`. Whether to generate the
#'   internal-domain conservation plot. `NULL` follows the selected output level.
#' @param top_n Integer. Number of top categories displayed in general plots
#'   and the dashboard.
#' @param top_n_motifs Integer. Number of top motifs displayed in detailed
#'   motif outputs.
#' @param top_n_subfamilies Integer. Number of top HERV subfamilies displayed
#'   in subfamily-level plots.
#' @param internal_domain_size_by Character. Bubble-size mode for the internal
#'   domain plot: `"fixed"` or `"domain_count"`.
#' @param verbose Logical. Whether to print progress messages.
#' @param show_progress Logical. Whether to display progress bars for
#'   long-running motif-file operations.
#'
#' @return A `HERVarium_annotation` object containing the compact HERV-level
#'   feature table, matched annotation, missing IDs, summary statistics,
#'   optional detailed regulatory results, output metadata, and generated-file
#'   paths.
#' @export
profile_hervs <- function(herv_ids,
                          annotation_file = NULL,
                          output_dir = "HERVariumR_results",

                          # Optional annotation layers
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

                          # Optional heavy TFBM layer
                          add_tfbm_details = FALSE,
                          fimo_file = NULL,
                          qvalue_cutoff = 1,

                          # Optional detailed RBP layer
                          add_rbp_details = FALSE,
                          rbp_file = NULL,
                          rbp_qvalue_cutoff = 1,

                          use_awk = TRUE,

                          # Output/report options
                          output_level = c("standard", "extended"),
                          make_dashboard = TRUE,
                          make_plots = NULL,
                          make_ltr_tfbm_plots = NULL,
                          make_internal_domain_plot = NULL,
                          top_n = 25,
                          top_n_motifs = 25,
                          top_n_subfamilies = 20,
                          internal_domain_size_by = c("fixed", "domain_count"),

                          # Verbosity
                          verbose = TRUE,
                          show_progress = TRUE) {

  output_level <- .match_output_level(output_level)
  internal_domain_size_by <- match.arg(internal_domain_size_by)

  if (is.null(make_plots)) {
    make_plots <- identical(output_level, "extended")
  }
  if (is.null(make_ltr_tfbm_plots)) {
    make_ltr_tfbm_plots <- identical(output_level, "extended")
  }
  if (is.null(make_internal_domain_plot)) {
    make_internal_domain_plot <- identical(output_level, "extended")
  }

  if (missing(herv_ids) || length(herv_ids) == 0) {
    stop("Please provide a vector of HERV IDs.")
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

  if (verbose) {
    message("Starting HERVariumR single-list profiling")
    message("Output level: ", output_level)
    message("Output directory: ", output_dir)
  }

  # ------------------------------------------------------------
  # 1. Core annotation (kept in memory; no intermediate files)
  # ------------------------------------------------------------

  if (verbose) message("\n[1/7] Annotating HERV list")

  anno <- annotate_hervs(
    herv_ids = herv_ids,
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

  # ------------------------------------------------------------
  # 2. External evidence
  # ------------------------------------------------------------

  if (add_external_evidence) {
    if (verbose) message("\n[2/7] Adding HERVarium external-evidence layers")

    anno <- add_encode_44organ_layers(
      anno,
      internal_encode_file = internal_encode_file,
      ltr_encode_file = ltr_encode_file,
      output_dir = NULL,
      verbose = verbose,
      .return_full_features = TRUE
    )
  } else if (verbose) {
    message("\n[2/7] Skipping HERVarium external-evidence layers")
  }

  # ------------------------------------------------------------
  # 3. LTR TF ChIP-seq evidence
  # ------------------------------------------------------------

  if (add_tfchip) {
    if (verbose) message("\n[3/7] Adding LTR TF ChIP-seq evidence layer")

    anno <- add_ltr_tfchip_layer(
      anno,
      tfchip_file = tfchip_file,
      output_dir = NULL,
      verbose = verbose,
      .return_full_features = TRUE
    )
  } else if (verbose) {
    message("\n[3/7] Skipping LTR TF ChIP-seq evidence layer")
  }

  # ------------------------------------------------------------
  # 4. Optional detailed LTR TFBM layer
  # ------------------------------------------------------------

  if (add_tfbm_details) {
    if (is.null(fimo_file)) {
      stop("add_tfbm_details = TRUE requires fimo_file.")
    }

    if (verbose) message("\n[4/7] Adding detailed LTR TFBM layer")

    tfbm_dir <- .make_transient_layer_dir(layout, "detailed_tfbm")
    if (!isTRUE(layout$write_details)) transient_dirs <- c(transient_dirs, tfbm_dir)

    anno <- add_ltr_tfbm_details(
      anno,
      fimo_file = fimo_file,
      output_dir = tfbm_dir,
      qvalue_cutoff = qvalue_cutoff,
      top_n_motifs = top_n_motifs,
      use_awk = use_awk,
      show_progress = show_progress,
      .return_full_features = TRUE
    )
  } else if (verbose) {
    message("\n[4/7] Skipping detailed LTR TFBM layer")
  }

  # ------------------------------------------------------------
  # 5. Optional detailed RBP layer
  # ------------------------------------------------------------

  if (add_rbp_details) {
    if (is.null(rbp_file)) {
      stop("add_rbp_details = TRUE requires rbp_file.")
    }

    if (verbose) message("\n[5/7] Adding detailed RBP layer")

    rbp_dir <- .make_transient_layer_dir(layout, "detailed_rbp")
    if (!isTRUE(layout$write_details)) transient_dirs <- c(transient_dirs, rbp_dir)

    anno <- add_rbp_details(
      anno,
      rbp_file = rbp_file,
      output_dir = rbp_dir,
      qvalue_cutoff = rbp_qvalue_cutoff,
      use_awk = use_awk,
      show_progress = show_progress,
      .return_full_features = TRUE
    )
  } else if (verbose) {
    message("\n[5/7] Skipping detailed RBP layer")
  }

  # ------------------------------------------------------------
  # 6. Final compact outputs
  # ------------------------------------------------------------

  if (verbose) message("\n[6/7] Writing compact standard outputs")

  full_features <- anno$features
  plot_anno <- anno
  plot_anno$features <- full_features
  anno$features <- select_compact_herv_features(full_features)
  anno$output_dir <- output_dir
  anno$output_level <- output_level

  compact_file <- file.path(output_dir, "herv_features_compact.tsv")
  object_file <- file.path(output_dir, "hervarium_annotation_object.rds")
  summary_file <- file.path(output_dir, "input_summary.tsv")
  missing_file <- file.path(output_dir, "missing_herv_ids.txt")

  .write_hervarium_tsv(anno$features, compact_file)

  count_bool <- function(col) {
    if (!col %in% colnames(anno$features)) return(NA_integer_)
    sum(.as_herv_bool(anno$features[[col]]), na.rm = TRUE)
  }

  input_summary <- data.frame(
    metric = c(
      "input_ids", "matched_hervs", "missing_ids",
      "with_internal_domains", "with_ltr5", "with_ltr3", "with_both_ltrs",
      "with_internal_external_evidence",
      "with_ltr5_external_evidence", "with_ltr3_external_evidence",
      "with_ltr5_tfchip", "with_ltr3_tfchip"
    ),
    value = c(
      anno$stats$n_input[1], nrow(anno$features), length(anno$missing_ids),
      count_bool("has_domain"), count_bool("has_ltr5"), count_bool("has_ltr3"),
      count_bool("has_both_ltrs"), count_bool("has_internal_encode_evidence"),
      count_bool("ltr5_encode_evidence_including_ccre"),
      count_bool("ltr3_encode_evidence_including_ccre"),
      count_bool("has_ltr5_tfchip_overlap"), count_bool("has_ltr3_tfchip_overlap")
    ),
    stringsAsFactors = FALSE
  )
  anno$input_summary <- input_summary
  .write_hervarium_tsv(input_summary, summary_file)

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

  anno$output_files <- list(
    compact_features = compact_file,
    input_summary = summary_file,
    object = object_file,
    dashboard = if (isTRUE(make_dashboard)) {
      file.path(output_dir, "hervarium_annotation_dashboard.html")
    } else {
      NULL
    },
    details_dir = if (isTRUE(layout$write_details)) layout$details_dir else NULL,
    plots_dir = if (isTRUE(make_plots)) layout$plots_dir else NULL
  )

  saveRDS(anno, object_file)

  # ------------------------------------------------------------
  # 7. Optional extended plots and concise dashboard
  # ------------------------------------------------------------

  if (isTRUE(make_plots)) {
    if (verbose) message("\n[7/7] Generating static extended plots")

    plot_herv_annotation(
      plot_anno,
      output_dir = layout$plots_dir,
      top_n = top_n
    )

    if (isTRUE(make_internal_domain_plot) && any(anno$features$has_domain, na.rm = TRUE)) {
      plot_internal_domain_bubbles(
        plot_anno,
        output_dir = layout$plots_dir,
        top_n_subfamilies = top_n_subfamilies,
        top_n_domains = top_n,
        size_by = internal_domain_size_by
      )
    }

    if (add_tfbm_details && isTRUE(make_ltr_tfbm_plots) &&
        !is.null(anno$ltr_tfbm_hits) && nrow(anno$ltr_tfbm_hits) > 0) {
      plot_ltr_tfbm(
        plot_anno,
        output_dir = layout$plots_dir,
        top_n_motifs = top_n_motifs,
        top_n_subfamilies = top_n_subfamilies
      )
    }
  } else if (verbose) {
    message("\n[7/7] Static plots disabled for standard output")
  }

  if (isTRUE(make_dashboard)) {
    generate_hervarium_dashboard(
      anno,
      output_dir = output_dir,
      report_file = "hervarium_annotation_dashboard.html",
      top_n = top_n
    )
  }

  # Save once more so the persisted object contains final output metadata.
  saveRDS(anno, object_file)

  if (verbose) {
    message("\nHERVariumR profiling finished.")
    message("Main table: ", compact_file)
    if (isTRUE(make_dashboard)) {
      message("Dashboard: ", file.path(output_dir, "hervarium_annotation_dashboard.html"))
    }
    message("Results written to: ", output_dir)
  }

  anno
}
