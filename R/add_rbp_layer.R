add_rbp_details <- function(x,
                            rbp_file,
                            output_dir = NULL,
                            qvalue_cutoff = 1,
                            use_awk = TRUE,
                            show_progress = TRUE,
                            .return_full_features = FALSE) {
  
  if (!file.exists(rbp_file)) {
    stop("RBP FIMO file not found: ", rbp_file)
  }
  
  if (inherits(x, "HERVarium_annotation")) {
    features <- x$features
  } else {
    features <- x
  }
  
  if (is.null(output_dir) && inherits(x, "HERVarium_annotation")) {
    output_dir <- x$output_dir
  }
  
  if (is.null(output_dir)) {
    stop("Please provide output_dir or use a HERVarium_annotation object with output_dir.")
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  required_cols <- c("HERV_id", "locid")
  missing_cols <- setdiff(required_cols, colnames(features))
  
  if (length(missing_cols) > 0) {
    stop(
      "The following required columns are missing from features: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  progress_total <- 6
  progress_step <- 0
  
  if (show_progress) {
    pb <- txtProgressBar(min = 0, max = progress_total, style = 3)
    on.exit(close(pb), add = TRUE)
  }
  
  update_progress <- function(label = NULL) {
    progress_step <<- progress_step + 1
    if (show_progress) {
      setTxtProgressBar(pb, progress_step)
      if (!is.null(label)) {
        message("\n", label)
      }
    }
  }
  
  update_progress("Preparing internal HERV list")
  
  rbp_map <- unique(
    features[
      !is.na(features$locid) &
        features$locid != "" &
        features$locid != ".",
      c("HERV_id", "locid")
    ]
  )
  
  colnames(rbp_map) <- c("HERV_id", "sequence_name")
  
  if (nrow(rbp_map) == 0) {
    warning("No valid locid values found in the annotation object.")
    return(x)
  }
  
  wanted_hervs_file <- tempfile(fileext = ".txt")
  writeLines(unique(rbp_map$sequence_name), wanted_hervs_file)
  
  filtered_rbp_file <- file.path(output_dir, "rbp_fimo_hits.tsv")
  
  update_progress("Filtering RBP FIMO file")
  
  if (use_awk) {
    
    # RBP_fimo.tsv columns:
    # motif_id motif_alt_id sequence_name start stop strand score p-value q-value matched_sequence
    # sequence_name is column 3.
    cmd <- paste(
      "awk 'BEGIN{FS=OFS=\"\\t\"}",
      "NR==FNR{keep[$1]=1; next}",
      "FNR==1{print; next}",
      "($3 in keep){print}'",
      shQuote(wanted_hervs_file),
      shQuote(rbp_file),
      ">",
      shQuote(filtered_rbp_file)
    )
    
    status <- system(cmd)
    
    if (status != 0) {
      stop("awk filtering failed. Try use_awk = FALSE.")
    }
    
    update_progress("Reading filtered RBP hits")
    
    rbp_hits <- read.delim(
      filtered_rbp_file,
      sep = "\t",
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
  } else {
    
    update_progress("Reading full RBP file with data.table")
    
    if (!requireNamespace("data.table", quietly = TRUE)) {
      stop("Package 'data.table' is required when use_awk = FALSE.")
    }
    
    rbp_all <- data.table::fread(
      rbp_file,
      sep = "\t",
      header = TRUE,
      data.table = FALSE
    )
    
    rbp_hits <- rbp_all[rbp_all$sequence_name %in% rbp_map$sequence_name, ]
    
    write.table(
      rbp_hits,
      file = filtered_rbp_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  
  if (nrow(rbp_hits) == 0) {
    warning("No RBP motif hits found for the selected HERVs.")
    return(x)
  }
  
  update_progress("Applying q-value filter and linking hits to HERVs")
  
  if ("q-value" %in% colnames(rbp_hits)) {
    rbp_hits <- rbp_hits[
      is.na(rbp_hits[["q-value"]]) |
        rbp_hits[["q-value"]] <= qvalue_cutoff,
    ]
  }
  
  rbp_hits <- merge(
    rbp_hits,
    rbp_map,
    by = "sequence_name",
    all.x = TRUE
  )
  
  write.table(
    rbp_hits,
    file = filtered_rbp_file,
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  collapse_unique_full <- function(z) {
    z <- unique(na.omit(as.character(z)))
    z <- z[z != "" & z != "."]
    if (length(z) == 0) return(NA_character_)
    paste(z, collapse = ";")
  }
  
  collapse_unique_preview <- function(z, max_items = 30) {
    z <- unique(na.omit(as.character(z)))
    z <- z[z != "" & z != "."]
    if (length(z) == 0) return(NA_character_)
    
    n_total <- length(z)
    
    if (n_total > max_items) {
      z <- c(z[seq_len(max_items)], paste0("...+", n_total - max_items, " more"))
    }
    
    paste(z, collapse = ";")
  }
  
  update_progress("Summarising RBP hits per HERV")
  
  rbp_split <- split(rbp_hits, rbp_hits$HERV_id)
  
  herv_rbp_summary <- do.call(
    rbind,
    lapply(names(rbp_split), function(herv) {
      z <- rbp_split[[herv]]
      
      rbp_ordered <- names(sort(table(z$motif_alt_id), decreasing = TRUE))
      
      data.frame(
        HERV_id = herv,
        n_rbp_hits_detailed = nrow(z),
        n_unique_rbp_motifs_detailed = length(unique(z$motif_id)),
        n_unique_rbp_names_detailed = length(unique(z$motif_alt_id)),
        rbp_names_all = collapse_unique_full(rbp_ordered),
        top_rbp_names_preview = collapse_unique_preview(rbp_ordered),
        best_rbp_p_value = suppressWarnings(min(z[["p-value"]], na.rm = TRUE)),
        best_rbp_q_value = if ("q-value" %in% colnames(z)) {
          suppressWarnings(min(z[["q-value"]], na.rm = TRUE))
        } else {
          NA_real_
        },
        best_rbp_score = suppressWarnings(max(z$score, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )
    })
  )
  
  rownames(herv_rbp_summary) <- NULL
  
  write.table(
    herv_rbp_summary,
    file = file.path(output_dir, "herv_rbp_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  update_progress("Adding RBP summaries to feature table")
  
  idx <- match(features$HERV_id, herv_rbp_summary$HERV_id)
  has_match <- !is.na(idx)
  
  features$n_rbp_hits_detailed <- 0
  features$n_unique_rbp_motifs_detailed <- 0
  features$n_unique_rbp_names_detailed <- 0
  features$rbp_names_all <- NA_character_
  features$top_rbp_names_preview <- NA_character_
  features$best_rbp_p_value <- NA_real_
  features$best_rbp_q_value <- NA_real_
  features$best_rbp_score <- NA_real_
  
  features$n_rbp_hits_detailed[has_match] <-
    herv_rbp_summary$n_rbp_hits_detailed[idx[has_match]]
  
  features$n_unique_rbp_motifs_detailed[has_match] <-
    herv_rbp_summary$n_unique_rbp_motifs_detailed[idx[has_match]]
  
  features$n_unique_rbp_names_detailed[has_match] <-
    herv_rbp_summary$n_unique_rbp_names_detailed[idx[has_match]]
  
  features$rbp_names_all[has_match] <-
    herv_rbp_summary$rbp_names_all[idx[has_match]]
  
  features$top_rbp_names_preview[has_match] <-
    herv_rbp_summary$top_rbp_names_preview[idx[has_match]]
  
  features$best_rbp_p_value[has_match] <-
    herv_rbp_summary$best_rbp_p_value[idx[has_match]]
  
  features$best_rbp_q_value[has_match] <-
    herv_rbp_summary$best_rbp_q_value[idx[has_match]]
  
  features$best_rbp_score[has_match] <-
    herv_rbp_summary$best_rbp_score[idx[has_match]]
  
  update_progress("Writing final RBP outputs")
  
  write.table(
    select_compact_herv_features(features),
    file = file.path(output_dir, "herv_features_compact.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  out_features <- if (.return_full_features) {
    features
  } else {
    select_compact_herv_features(features)
  }

  if (inherits(x, "HERVarium_annotation")) {
    x$features <- out_features
    x$rbp_hits <- rbp_hits
    x$herv_rbp_summary <- herv_rbp_summary
    
    message("Detailed RBP annotations written to: ", output_dir)
    return(x)
  }
  
  out_features
}