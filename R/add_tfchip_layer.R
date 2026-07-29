#' Add LTR TF ChIP-seq overlap evidence
#'
#' Adds compact TF ChIP-seq peak-cluster overlap annotations for the associated
#' 5′ and 3′ LTRs of each HERV. The layer is joined through the LTR names stored
#' in the HERVariumR feature table.
#'
#' TF ChIP-seq overlap is experimental positional context. It is distinct from
#' sequence-based TFBM prediction and does not by itself prove autonomous TF
#' recruitment by a HERV sequence.
#'
#' @param x A `HERVarium_annotation` object or HERV-level `data.frame`.
#' @param tfchip_file Optional path to `ltr_tfchip_compact.tsv.gz`. If `NULL`,
#'   the bundled resource is used.
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
add_ltr_tfchip_layer <- function(x,
                                 tfchip_file = NULL,
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

  if (is.null(tfchip_file)) {
    tfchip_file <- system.file(
      "extdata",
      "ltr_tfchip_compact.tsv.gz",
      package = "HERVariumR"
    )
  }

  if (tfchip_file == "" || !file.exists(tfchip_file)) {
    stop(
      "LTR TF ChIP-seq compact evidence file not found. Expected: ",
      "ltr_tfchip_compact.tsv.gz"
    )
  }

  if (verbose) {
    message("Adding HERVarium LTR TF ChIP-seq evidence layer")
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

  row_sum_numeric <- function(...) {
    z <- data.frame(..., check.names = FALSE)
    z[] <- lapply(z, function(v) suppressWarnings(as.numeric(v)))
    apply(z, 1, function(v) {
      if (all(is.na(v))) return(NA_real_)
      sum(v, na.rm = TRUE)
    })
  }

  row_max_numeric <- function(...) {
    z <- data.frame(..., check.names = FALSE)
    z[] <- lapply(z, function(v) suppressWarnings(as.numeric(v)))
    apply(z, 1, function(v) {
      if (all(is.na(v))) return(NA_real_)
      max(v, na.rm = TRUE)
    })
  }

  split_tf_names <- function(x) {
    x <- as.character(x)
    x <- x[!is.na(x) & x != "" & x != "."]
    if (length(x) == 0) return(character(0))

    x <- paste(x, collapse = ",")
    z <- unlist(strsplit(x, "[,;]", perl = TRUE))
    z <- trimws(z)
    z <- z[!is.na(z) & z != "" & z != "."]
    unique(z)
  }

  collapse_tf_names <- function(..., max_items = 30) {
    vals <- list(...)
    z <- unique(unlist(lapply(vals, split_tf_names)))
    z <- z[!is.na(z) & z != "" & z != "."]
    if (length(z) == 0) return(NA_character_)

    n_total <- length(z)
    if (n_total > max_items) {
      z <- c(z[seq_len(max_items)], paste0("...+", n_total - max_items, " more"))
    }
    paste(z, collapse = ",")
  }

  count_tf_names <- function(...) {
    vals <- list(...)
    z <- unique(unlist(lapply(vals, split_tf_names)))
    z <- z[!is.na(z) & z != "" & z != "."]
    length(z)
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

  make_ltr_col_name <- function(prefix, col) {
    if (col == "tfchip_overlap_detected") {
      return(paste0(prefix, "_overlap"))
    }

    if (col == "tfchip_overlap") {
      return(paste0(prefix, "_overlap"))
    }

    if (grepl("^tfchip_", col)) {
      return(paste0(prefix, "_", sub("^tfchip_", "", col)))
    }

    paste0(prefix, "_", col)
  }

  add_ltr_table <- function(features, tfchip_tbl, ltr_col, prefix) {
    idx <- match(features[[ltr_col]], tfchip_tbl$sequence_name)

    skip_cols <- c(
      "sequence_name",
      "tfchip_source_sequence_name",
      "chrom",
      "ltr_start",
      "ltr_end",
      "ltr_strand",
      "start",
      "end",
      "strand"
    )

    add_cols <- setdiff(colnames(tfchip_tbl), skip_cols)

    for (col in add_cols) {
      new_col <- make_ltr_col_name(prefix, col)
      features[[new_col]] <- tfchip_tbl[[col]][idx]
    }

    features
  }

  choose_top_tf <- function(tf5, bp5, tf3, bp3) {
    bp5 <- suppressWarnings(as.numeric(bp5))
    bp3 <- suppressWarnings(as.numeric(bp3))

    out <- rep(NA_character_, length(tf5))

    use5 <- !is.na(bp5) & (is.na(bp3) | bp5 >= bp3) & !is.na(tf5) & tf5 != ""
    use3 <- !is.na(bp3) & (is.na(bp5) | bp3 > bp5) & !is.na(tf3) & tf3 != ""

    out[use5] <- as.character(tf5[use5])
    out[use3] <- as.character(tf3[use3])
    out
  }

  # ------------------------------------------------------------
  # 1. Read and standardize compact TF ChIP-seq asset
  # ------------------------------------------------------------

  tfchip <- read_compact(tfchip_file)

  if (!"sequence_name" %in% colnames(tfchip)) {
    stop("tfchip_file must contain a sequence_name column.")
  }

  if ("tfchip_overlap_detected" %in% colnames(tfchip) &&
      !"tfchip_overlap" %in% colnames(tfchip)) {
    tfchip$tfchip_overlap <- tfchip$tfchip_overlap_detected
  }

  # Normalize the canonical TF-list field at the source-layer boundary.
  # Some compact assets contain both a truly complete list
  # (`tfchip_tf_list_full`) and shortened display/list variants containing
  # strings such as "...(+41)". Downstream code should always use the full
  # list whenever it is available.
  looks_truncated_tf_list <- function(x) {
    x <- as.character(x)
    !is.na(x) & grepl("\\.\\.\\.\\(\\+?[0-9]+(?: more)?\\)", x, perl = TRUE)
  }

  tf_list_source <- if ("tfchip_tf_list_full" %in% colnames(tfchip)) {
    "tfchip_tf_list_full"
  } else if ("tfchip_tf_list" %in% colnames(tfchip)) {
    "tfchip_tf_list"
  } else if ("tfchip_tf_list_display" %in% colnames(tfchip)) {
    "tfchip_tf_list_display"
  } else {
    NA_character_
  }

  if (!is.na(tf_list_source)) {
    tfchip$tfchip_tf_list <- tfchip[[tf_list_source]]
  }

  if ("tfchip_tf_list_full" %in% colnames(tfchip) &&
      any(looks_truncated_tf_list(tfchip$tfchip_tf_list_full), na.rm = TRUE)) {
    warning(
      "The bundled/custom tfchip_tf_list_full column still contains truncated values. ",
      "The compact TF ChIP-seq asset must be rebuilt from the original source table ",
      "to recover complete TF lists.",
      call. = FALSE
    )
  } else if (!is.na(tf_list_source) && tf_list_source != "tfchip_tf_list_full" &&
             any(looks_truncated_tf_list(tfchip$tfchip_tf_list), na.rm = TRUE)) {
    warning(
      "No complete tfchip_tf_list_full column was available; TF ChIP-seq lists may remain truncated.",
      call. = FALSE
    )
  }

  # ------------------------------------------------------------
  # 2. Add 5' and 3' LTR TF ChIP-seq evidence
  # ------------------------------------------------------------

  features <- add_ltr_table(
    features = features,
    tfchip_tbl = tfchip,
    ltr_col = "ltr5_name",
    prefix = "ltr5_tfchip"
  )

  features <- add_ltr_table(
    features = features,
    tfchip_tbl = tfchip,
    ltr_col = "ltr3_name",
    prefix = "ltr3_tfchip"
  )

  # Ensure boolean LTR overlap columns exist even if the source naming changes.
  if (!"ltr5_tfchip_overlap" %in% colnames(features)) {
    features$ltr5_tfchip_overlap <- FALSE
  }

  if (!"ltr3_tfchip_overlap" %in% colnames(features)) {
    features$ltr3_tfchip_overlap <- FALSE
  }

  features$ltr5_tfchip_overlap <- as_bool(features$ltr5_tfchip_overlap)
  features$ltr3_tfchip_overlap <- as_bool(features$ltr3_tfchip_overlap)

  # ------------------------------------------------------------
  # 3. HERV-level combined LTR TF ChIP-seq summaries
  # ------------------------------------------------------------

  features$has_ltr5_tfchip_overlap <- features$ltr5_tfchip_overlap
  features$has_ltr3_tfchip_overlap <- features$ltr3_tfchip_overlap

  features$has_any_ltr_tfchip_overlap <- row_any_bool(
    features$ltr5_tfchip_overlap,
    features$ltr3_tfchip_overlap
  )

  features$has_both_ltr_tfchip_overlap <- row_all_bool(
    features$ltr5_tfchip_overlap,
    features$ltr3_tfchip_overlap
  )

  if (all(c("ltr5_tfchip_n_peak_clusters", "ltr3_tfchip_n_peak_clusters") %in% colnames(features))) {
    features$total_ltr_tfchip_n_peak_clusters <- row_sum_numeric(
      features$ltr5_tfchip_n_peak_clusters,
      features$ltr3_tfchip_n_peak_clusters
    )

    features$any_ltr_tfchip_n_peak_clusters_max <- row_max_numeric(
      features$ltr5_tfchip_n_peak_clusters,
      features$ltr3_tfchip_n_peak_clusters
    )
  }

  if (all(c("ltr5_tfchip_n_tfs", "ltr3_tfchip_n_tfs") %in% colnames(features))) {
    features$total_ltr_tfchip_n_tfs <- row_sum_numeric(
      features$ltr5_tfchip_n_tfs,
      features$ltr3_tfchip_n_tfs
    )

    features$any_ltr_tfchip_n_tfs_max <- row_max_numeric(
      features$ltr5_tfchip_n_tfs,
      features$ltr3_tfchip_n_tfs
    )
  }

  if (all(c("ltr5_tfchip_total_overlap_bp", "ltr3_tfchip_total_overlap_bp") %in% colnames(features))) {
    features$total_ltr_tfchip_overlap_bp <- row_sum_numeric(
      features$ltr5_tfchip_total_overlap_bp,
      features$ltr3_tfchip_total_overlap_bp
    )
  }

  if (all(c("ltr5_tfchip_max_fraction_of_region", "ltr3_tfchip_max_fraction_of_region") %in% colnames(features))) {
    features$any_ltr_tfchip_max_fraction_of_region <- row_max_numeric(
      features$ltr5_tfchip_max_fraction_of_region,
      features$ltr3_tfchip_max_fraction_of_region
    )
  }

  if (all(c("ltr5_tfchip_max_fraction_of_peak", "ltr3_tfchip_max_fraction_of_peak") %in% colnames(features))) {
    features$any_ltr_tfchip_max_fraction_of_peak <- row_max_numeric(
      features$ltr5_tfchip_max_fraction_of_peak,
      features$ltr3_tfchip_max_fraction_of_peak
    )
  }

  if (all(c("ltr5_tfchip_max_cluster_support", "ltr3_tfchip_max_cluster_support") %in% colnames(features))) {
    features$any_ltr_tfchip_max_cluster_support <- row_max_numeric(
      features$ltr5_tfchip_max_cluster_support,
      features$ltr3_tfchip_max_cluster_support
    )
  }

  if (all(c("ltr5_tfchip_max_biosamples", "ltr3_tfchip_max_biosamples") %in% colnames(features))) {
    features$any_ltr_tfchip_max_biosamples <- row_max_numeric(
      features$ltr5_tfchip_max_biosamples,
      features$ltr3_tfchip_max_biosamples
    )
  }

  # Union of TF names across both LTRs. Prefer the explicit complete list,
  # then the normalized canonical list.
  tf5_col <- if ("ltr5_tfchip_tf_list_full" %in% colnames(features)) {
    "ltr5_tfchip_tf_list_full"
  } else if ("ltr5_tfchip_tf_list" %in% colnames(features)) {
    "ltr5_tfchip_tf_list"
  } else {
    NA_character_
  }

  tf3_col <- if ("ltr3_tfchip_tf_list_full" %in% colnames(features)) {
    "ltr3_tfchip_tf_list_full"
  } else if ("ltr3_tfchip_tf_list" %in% colnames(features)) {
    "ltr3_tfchip_tf_list"
  } else {
    NA_character_
  }

  if (!is.na(tf5_col) && !is.na(tf3_col)) {
    features$any_ltr_tfchip_n_tfs <- vapply(
      seq_len(nrow(features)),
      function(i) count_tf_names(features[[tf5_col]][i], features[[tf3_col]][i]),
      integer(1)
    )

    features$any_ltr_tfchip_tf_list_display <- vapply(
      seq_len(nrow(features)),
      function(i) collapse_tf_names(features[[tf5_col]][i], features[[tf3_col]][i]),
      character(1)
    )
  }

  if (all(c(
    "ltr5_tfchip_top_tf_by_overlap",
    "ltr5_tfchip_top_tf_overlap_bp",
    "ltr3_tfchip_top_tf_by_overlap",
    "ltr3_tfchip_top_tf_overlap_bp"
  ) %in% colnames(features))) {
    features$any_ltr_tfchip_top_tf_by_overlap <- choose_top_tf(
      tf5 = features$ltr5_tfchip_top_tf_by_overlap,
      bp5 = features$ltr5_tfchip_top_tf_overlap_bp,
      tf3 = features$ltr3_tfchip_top_tf_by_overlap,
      bp3 = features$ltr3_tfchip_top_tf_overlap_bp
    )

    features$any_ltr_tfchip_top_tf_overlap_bp <- row_max_numeric(
      features$ltr5_tfchip_top_tf_overlap_bp,
      features$ltr3_tfchip_top_tf_overlap_bp
    )
  }

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
    x$tfchip_file <- tfchip_file

    if (verbose) {
      message("LTR TF ChIP-seq layer added to HERVarium_annotation object.")
    }

    return(x)
  }

  out_features
}
