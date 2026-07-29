#' Add HERVarium ENCODE/UCSC 44-organ external-evidence layers
#'
#' Adds compact internal-region and LTR-level ENCODE/UCSC 44-organ evidence to
#' a `HERVarium_annotation` object or HERV-level feature table.
#'
#' Internal evidence is joined by `HERV_id` through the compact precomputed
#' bridge. Associated LTR evidence is joined separately through `ltr5_name` and
#' `ltr3_name`. In the signal hierarchy, high evidence is a subset of robust
#' evidence, which is a subset of detected evidence.
#'
#' @param x A `HERVarium_annotation` object or HERV-level `data.frame`.
#' @param internal_encode_file Optional path to
#'   `internal_encode_44organs_compact.tsv.gz`. If `NULL`, the bundled resource
#'   is used.
#' @param ltr_encode_file Optional path to
#'   `ltr_encode_44organs_compact.tsv.gz`. If `NULL`, the bundled resource is
#'   used.
#' @param output_dir Optional output directory. For a `HERVarium_annotation`
#'   object, `x$output_dir` is used when available.
#' @param verbose Logical. Whether to print progress messages.
#' @param .return_full_features Internal use only. Whether to retain and return
#'   the full intermediate feature matrix instead of applying compact
#'   user-facing pruning.
#'
#' @return The updated object or feature `data.frame`. By default, its
#'   user-facing feature table is compact; `.return_full_features = TRUE` is
#'   reserved for internal workflow composition.
#' @export
add_encode_44organ_layers <- function(x,
                                      internal_encode_file = NULL,
                                      ltr_encode_file = NULL,
                                      output_dir = NULL,
                                      verbose = TRUE,
                                      .return_full_features = FALSE) {

  if (inherits(x, "HERVarium_annotation")) {
    features <- x$features

    if (is.null(output_dir)) {
      output_dir <- x$output_dir
    }
  } else {
    features <- x
  }

  if (!is.data.frame(features)) {
    stop("x must be a HERVarium_annotation object or a feature data.frame.")
  }

  required_feature_cols <- c("HERV_id", "ltr5_name", "ltr3_name")
  missing_feature_cols <- setdiff(required_feature_cols, colnames(features))

  if (length(missing_feature_cols) > 0) {
    stop(
      "The following required feature columns are missing: ",
      paste(missing_feature_cols, collapse = ", ")
    )
  }

  if (is.null(internal_encode_file)) {
    internal_encode_file <- system.file(
      "extdata",
      "internal_encode_44organs_compact.tsv.gz",
      package = "HERVariumR"
    )
  }

  if (is.null(ltr_encode_file)) {
    ltr_encode_file <- system.file(
      "extdata",
      "ltr_encode_44organs_compact.tsv.gz",
      package = "HERVariumR"
    )
  }

  if (internal_encode_file == "" || !file.exists(internal_encode_file)) {
    stop(
      "Internal ENCODE/UCSC compact evidence file not found. Expected: ",
      "internal_encode_44organs_compact.tsv.gz"
    )
  }

  if (ltr_encode_file == "" || !file.exists(ltr_encode_file)) {
    stop(
      "LTR ENCODE/UCSC compact evidence file not found. Expected: ",
      "ltr_encode_44organs_compact.tsv.gz"
    )
  }

  if (verbose) {
    message("Adding HERVarium ENCODE/UCSC 44-organ external-evidence layers")
  }

  # ------------------------------------------------------------
  # Helpers
  # ------------------------------------------------------------

  as_bool <- function(x) {
    if (is.logical(x)) {
      x[is.na(x)] <- FALSE
      return(x)
    }

    x <- as.character(x)
    out <- tolower(x) %in% c("true", "t", "1", "yes")
    out[is.na(out)] <- FALSE
    out
  }

  row_any_bool <- function(...) {
    mats <- list(...)
    if (length(mats) == 0) return(logical(0))
    Reduce(`|`, lapply(mats, as_bool))
  }

  row_all_bool <- function(...) {
    mats <- list(...)
    if (length(mats) == 0) return(logical(0))
    Reduce(`&`, lapply(mats, as_bool))
  }

  row_max_numeric <- function(...) {
    z <- data.frame(..., check.names = FALSE)
    z[] <- lapply(z, function(v) suppressWarnings(as.numeric(v)))
    apply(z, 1, function(v) {
      if (all(is.na(v))) return(NA_real_)
      max(v, na.rm = TRUE)
    })
  }

  row_sum_numeric <- function(...) {
    z <- data.frame(..., check.names = FALSE)
    z[] <- lapply(z, function(v) suppressWarnings(as.numeric(v)))
    apply(z, 1, function(v) {
      if (all(is.na(v))) return(NA_real_)
      sum(v, na.rm = TRUE)
    })
  }

  read_compact <- function(path) {
    read.delim(
      path,
      sep = "\t",
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE,
      na.strings = c("NA", "")
    )
  }

  add_ltr_table <- function(features, ltr_tbl, ltr_col, prefix) {
    idx <- match(features[[ltr_col]], ltr_tbl$sequence_name)

    skip_cols <- c("sequence_name", "chrom", "start", "end", "strand", "region_length")
    add_cols <- setdiff(colnames(ltr_tbl), skip_cols)

    for (col in add_cols) {
      new_col <- paste0(prefix, "_", col)
      features[[new_col]] <- ltr_tbl[[col]][idx]
    }

    features
  }

  # ------------------------------------------------------------
  # 1. Internal evidence
  # ------------------------------------------------------------

  internal <- read_compact(internal_encode_file)

  if (!"HERV_id" %in% colnames(internal)) {
    stop("internal_encode_file must contain a HERV_id column.")
  }

  internal_add_cols <- setdiff(
    colnames(internal),
    c("HERV_id", "locid", "chrom", "start", "end", "strand")
  )

  idx_internal <- match(features$HERV_id, internal$HERV_id)

  for (col in internal_add_cols) {
    features[[col]] <- internal[[col]][idx_internal]
  }

  features$has_internal_encode_evidence <-
    as_bool(features$internal_encode_any_evidence)

  features$has_internal_encode_evidence_including_ccre <-
    as_bool(features$internal_encode_any_evidence_including_ccre)

  features$has_internal_multilayer_encode_evidence <-
    as_bool(features$internal_encode_multilayer_detected_evidence)

  features$has_internal_robust_multilayer_encode_evidence <-
    as_bool(features$internal_encode_multilayer_robust_evidence)

  features$has_internal_ccre_overlap <-
    as_bool(features$internal_ccre_overlap)

  # ------------------------------------------------------------
  # 2. LTR evidence
  # ------------------------------------------------------------

  ltr <- read_compact(ltr_encode_file)

  if (!"sequence_name" %in% colnames(ltr)) {
    stop("ltr_encode_file must contain a sequence_name column.")
  }

  features <- add_ltr_table(
    features = features,
    ltr_tbl = ltr,
    ltr_col = "ltr5_name",
    prefix = "ltr5_encode"
  )

  features <- add_ltr_table(
    features = features,
    ltr_tbl = ltr,
    ltr_col = "ltr3_name",
    prefix = "ltr3_encode"
  )

  # ------------------------------------------------------------
  # 3. HERV-level combined LTR summaries
  # ------------------------------------------------------------

  features$has_ltr5_encode_evidence <-
    suppressWarnings(as.numeric(features$ltr5_encode_detected_layer_count)) > 0

  features$has_ltr3_encode_evidence <-
    suppressWarnings(as.numeric(features$ltr3_encode_detected_layer_count)) > 0

  features$has_any_ltr_encode_evidence <- row_any_bool(
    features$has_ltr5_encode_evidence,
    features$has_ltr3_encode_evidence
  )

  features$has_both_ltr_encode_evidence <- row_all_bool(
    features$has_ltr5_encode_evidence,
    features$has_ltr3_encode_evidence
  )

  features$has_any_ltr_encode_evidence_including_ccre <- row_any_bool(
    features$has_ltr5_encode_evidence,
    features$has_ltr3_encode_evidence,
    features$ltr5_encode_ccre_overlap,
    features$ltr3_encode_ccre_overlap
  )

  features$has_any_ltr_ccre_overlap <- row_any_bool(
    features$ltr5_encode_ccre_overlap,
    features$ltr3_encode_ccre_overlap
  )

  features$has_both_ltr_ccre_overlap <- row_all_bool(
    features$ltr5_encode_ccre_overlap,
    features$ltr3_encode_ccre_overlap
  )

  features$any_ltr_encode_detected_layer_count_max <- row_max_numeric(
    features$ltr5_encode_detected_layer_count,
    features$ltr3_encode_detected_layer_count
  )

  features$any_ltr_encode_robust_layer_count_max <- row_max_numeric(
    features$ltr5_encode_robust_layer_count,
    features$ltr3_encode_robust_layer_count
  )

  features$any_ltr_encode_high_layer_count_max <- row_max_numeric(
    features$ltr5_encode_high_layer_count,
    features$ltr3_encode_high_layer_count
  )

  features$total_ltr_encode_detected_layer_count <- row_sum_numeric(
    features$ltr5_encode_detected_layer_count,
    features$ltr3_encode_detected_layer_count
  )

  features$total_ltr_ccre_n_overlaps <- row_sum_numeric(
    features$ltr5_encode_ccre_n_overlaps,
    features$ltr3_encode_ccre_n_overlaps
  )

  # Add compact LTR-specific evidence columns used by default outputs
  # and compact comparison tests.
  features <- .add_compact_derived_columns(features)

  # ------------------------------------------------------------
  # 4. Write updated feature table if possible
  # ------------------------------------------------------------

  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }

    write.table(
      select_compact_herv_features(features),
      file = file.path(output_dir, "herv_features_compact.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }

  out_features <- if (.return_full_features) {
    features
  } else {
    select_compact_herv_features(features)
  }

  if (inherits(x, "HERVarium_annotation")) {
    x$features <- out_features
    x$external_evidence_files <- list(
      internal_encode_file = internal_encode_file,
      ltr_encode_file = ltr_encode_file
    )

    if (verbose) {
      message("External-evidence layers added to HERVarium_annotation object.")
    }

    return(x)
  }

  out_features
}
