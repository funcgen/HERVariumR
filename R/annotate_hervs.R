#' Annotate a list of HERV IDs
#'
#' Performs the core HERVariumR annotation step for a list of HERV IDs. This
#' function loads bundled HERVariumR annotation resources by default, filters
#' internal retroviral domains according to a coverage threshold, adds bundled
#' annotation layers such as IFN-related LTR features and terminal-exon domain
#' context, and returns a `HERVarium_annotation` object.
#'
#' @param herv_ids Character vector of HERV IDs to annotate.
#' @param annotation_file Optional path to the main transcript-context annotation file. If `NULL`, the bundled HERVariumR annotation file is used.
#' @param id_columns Character vector of length two giving the columns used to match input IDs. Defaults to `c("HERV_id", "locid")`.
#' @param ifn_stat1_file Optional path to the STAT1 LTR summary file. If `NULL`, the bundled HERVariumR file is used.
#' @param ifn_stat1stat2_irf_file Optional path to the STAT1/STAT2/IRF LTR summary file. If `NULL`, the bundled HERVariumR file is used.
#' @param last_exon_file Optional path to the terminal-exon domain annotation file. If `NULL`, the bundled HERVariumR file is used.
#' @param domain_coverage_cutoff Numeric. Minimum HMM profile coverage required for an internal retroviral domain to be considered valid.
#' @param keep_raw_domain_columns Logical. Whether to keep unfiltered domain annotation columns with `_raw` suffixes.
#' @param output_dir Optional output directory where annotation tables will be written.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @return A `HERVarium_annotation` object containing matched annotation rows, compact feature summaries, missing IDs, summary statistics, and output metadata.
#' @export
annotate_hervs <- function(herv_ids,
                           annotation_file = NULL,
                           id_columns = c("HERV_id", "locid"),
                           ifn_stat1_file = NULL,
                           ifn_stat1stat2_irf_file = NULL,
                           last_exon_file = NULL,
                           domain_coverage_cutoff = 0.40,
                           keep_raw_domain_columns = TRUE,
                           output_dir = NULL,
                           verbose = TRUE) {
  
  hervarium_file <- function(filename) {
    path <- system.file("extdata", filename, package = "HERVariumR")
    
    if (path == "") {
      stop(
        "Could not find bundled HERVariumR file: ", filename,
        call. = FALSE
      )
    }
    
    path
  }
  
  
  .get_default_annotation_file <- function(annotation_file = NULL) {
    if (!is.null(annotation_file)) {
      return(annotation_file)
    }
    
    hervarium_file("transcript_context.with_herv_id.tsv.gz")
  }
  
  
  .get_default_ifn_stat1_file <- function(ifn_stat1_file = NULL) {
    if (!is.null(ifn_stat1_file)) {
      return(ifn_stat1_file)
    }
    
    hervarium_file("LTR_IFN_STAT1_summary.tsv")
  }
  
  
  .get_default_ifn_stat1stat2_irf_file <- function(ifn_stat1stat2_irf_file = NULL) {
    if (!is.null(ifn_stat1stat2_irf_file)) {
      return(ifn_stat1stat2_irf_file)
    }
    
    hervarium_file("LTR_IFN_STAT1STAT2_IRF_summary.tsv")
  }
  
  
  .get_default_last_exon_file <- function(last_exon_file = NULL) {
    if (!is.null(last_exon_file)) {
      return(last_exon_file)
    }
    
    hervarium_file("HERV_domains_transcript_context_last_exon.xlsx")
  }
  
  if (missing(herv_ids) || length(herv_ids) == 0) {
    stop("Please provide a vector of HERV IDs.")
  }
  
  herv_ids <- unique(as.character(herv_ids))
  
  annotation_file <- .get_default_annotation_file(annotation_file)
  ifn_stat1_file <- .get_default_ifn_stat1_file(ifn_stat1_file)
  ifn_stat1stat2_irf_file <- .get_default_ifn_stat1stat2_irf_file(ifn_stat1stat2_irf_file)
  last_exon_file <- .get_default_last_exon_file(last_exon_file)
  
  annot <- load_transcript_context(annotation_file)
  
  missing_cols <- setdiff(id_columns, colnames(annot))
  if (length(missing_cols) > 0) {
    stop(
      "The following ID columns are missing from the annotation table: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  matched <- annot[
    annot[[id_columns[1]]] %in% herv_ids |
      annot[[id_columns[2]]] %in% herv_ids,
  ]
  
  matched <- filter_herv_domains_by_coverage(
    matched,
    coverage_cutoff = domain_coverage_cutoff,
    keep_raw_domain_columns = keep_raw_domain_columns
  )
  
  features <- summarize_herv_features(matched)
  
  features <- add_ltr_ifn_annotations(
    features = features,
    stat1_file = ifn_stat1_file,
    stat1stat2_irf_file = ifn_stat1stat2_irf_file
  )
  
  features <- add_terminal_exon_domain_annotations(
    features = features,
    last_exon_file = last_exon_file
  )
  
  features <- apply_terminal_domain_coverage_cutoff(
    features,
    coverage_cutoff = domain_coverage_cutoff
  )
  
  matched_ids <- unique(c(
    matched[[id_columns[1]]],
    matched[[id_columns[2]]]
  ))
  
  missing_ids <- setdiff(herv_ids, matched_ids)
  
  stats <- data.frame(
    n_input = length(herv_ids),
    n_matched_rows = nrow(matched),
    n_unique_hervs = length(unique(matched$HERV_id)),
    n_missing = length(missing_ids),
    domain_coverage_cutoff = domain_coverage_cutoff,
    n_with_domains = sum(matched$domain_count > 0, na.rm = TRUE),
    n_with_ltr5 = sum(matched$has_ltr5 == 1, na.rm = TRUE),
    n_with_ltr3 = sum(matched$has_ltr3 == 1, na.rm = TRUE),
    n_with_any_ltr = sum(matched$has_ltr5 == 1 | matched$has_ltr3 == 1, na.rm = TRUE),
    mean_rbp_burden = mean(matched$rbp_burden, na.rm = TRUE),
    mean_ltr5_tfbm_burden = mean(matched$ltr5_tfbm_burden, na.rm = TRUE),
    mean_ltr3_tfbm_burden = mean(matched$ltr3_tfbm_burden, na.rm = TRUE)
  )
  
  if (verbose) {
    message("Input HERV IDs: ", length(herv_ids))
    message("Matched rows: ", nrow(matched))
    message("Unique matched HERVs: ", length(unique(matched$HERV_id)))
    message("Missing IDs: ", length(missing_ids))
  }
  
  if (!is.null(output_dir)) {
    if (!dir.exists(output_dir)) {
      dir.create(output_dir, recursive = TRUE)
    }
    
    matched_report <- matched
    
    if ("domains_gene" %in% colnames(matched_report)) {
      matched_report$domains_gene <- NULL
    }
    
    write.table(
      matched_report,
      file = file.path(output_dir, "herv_annotation_summary.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    write.table(
      features,
      file = file.path(output_dir, "herv_features_compact.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    write.table(
      stats,
      file = file.path(output_dir, "annotation_stats.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    writeLines(
      missing_ids,
      con = file.path(output_dir, "missing_herv_ids.txt")
    )
    
    if (verbose) {
      message("Results written to: ", output_dir)
    }
  }
  
  out <- list(
    summary = matched,
    features = features,
    missing_ids = missing_ids,
    stats = stats,
    output_dir = output_dir
  )
  out$domain_coverage_cutoff <- domain_coverage_cutoff
  
  class(out) <- "HERVarium_annotation"
  
  return(out)
}


print.HERVarium_annotation <- function(x, ...) {
  cat("HERVarium annotation result\n")
  cat("---------------------------\n")
  print(x$stats)
  
  if (!is.null(x$output_dir)) {
    cat("\nOutput directory:\n")
    cat(x$output_dir, "\n")
  }
  
  if (length(x$missing_ids) > 0) {
    cat("\nMissing IDs:\n")
    print(x$missing_ids)
  }
  
  invisible(x)
}