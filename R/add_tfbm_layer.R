add_ltr_tfbm_details <- function(x,
                                 fimo_file,
                                 output_dir = NULL,
                                 qvalue_cutoff = 1,
                                 top_n_motifs = 25,
                                 use_awk = TRUE,
                                 show_progress = TRUE,
                                 .return_full_features = FALSE) {
  
  if (!file.exists(fimo_file)) {
    stop("FIMO file not found: ", fimo_file)
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
  
  progress_total <- 8
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
  
  update_progress("Preparing LTR list")
  
  ltr_map <- data.frame(
    HERV_id = c(features$HERV_id, features$HERV_id),
    ltr_position = c(rep("LTR5", nrow(features)), rep("LTR3", nrow(features))),
    sequence_name = c(features$ltr5_name, features$ltr3_name),
    stringsAsFactors = FALSE
  )
  
  ltr_map <- ltr_map[
    !is.na(ltr_map$sequence_name) &
      ltr_map$sequence_name != "." &
      ltr_map$sequence_name != "",
  ]
  
  ltr_map <- unique(ltr_map)
  
  if (nrow(ltr_map) == 0) {
    warning("No valid LTR names found in the annotation object.")
    return(x)
  }
  
  wanted_ltrs_file <- tempfile(fileext = ".txt")
  writeLines(unique(ltr_map$sequence_name), wanted_ltrs_file)
  
  filtered_fimo_file <- file.path(output_dir, "ltr_tfbm_hits.tsv")
  
  update_progress("Filtering large FIMO file. This can take some time.")
  
  if (use_awk) {
    cmd <- paste(
      "awk 'BEGIN{FS=OFS=\"\\t\"}",
      "NR==FNR{keep[$1]=1; next}",
      "FNR==1{print; next}",
      "($4 in keep){print}'",
      shQuote(wanted_ltrs_file),
      shQuote(fimo_file),
      ">",
      shQuote(filtered_fimo_file)
    )
    
    status <- system(cmd)
    
    if (status != 0) {
      stop("awk filtering failed. Try use_awk = FALSE.")
    }
    
    update_progress("Reading filtered FIMO hits")
    
    tfbm_hits <- read.delim(
      filtered_fimo_file,
      sep = "\t",
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
  } else {
    update_progress("Reading full FIMO file with data.table")
    
    if (!requireNamespace("data.table", quietly = TRUE)) {
      stop("Package 'data.table' is required when use_awk = FALSE.")
    }
    
    fimo <- data.table::fread(
      fimo_file,
      sep = "\t",
      header = TRUE,
      data.table = FALSE
    )
    
    tfbm_hits <- fimo[fimo$sequence_name %in% ltr_map$sequence_name, ]
    
    write.table(
      tfbm_hits,
      file = filtered_fimo_file,
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
  }
  
  if (nrow(tfbm_hits) == 0) {
    warning("No TFBM hits found for the selected LTRs.")
    return(x)
  }
  
  update_progress("Applying q-value filter and linking hits to HERVs")
  
  if ("q-value" %in% colnames(tfbm_hits)) {
    tfbm_hits <- tfbm_hits[
      is.na(tfbm_hits[["q-value"]]) |
        tfbm_hits[["q-value"]] <= qvalue_cutoff,
    ]
  }
  
  tfbm_hits <- merge(
    tfbm_hits,
    ltr_map,
    by = "sequence_name",
    all.x = TRUE
  )
  
  write.table(
    tfbm_hits,
    file = filtered_fimo_file,
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
  
  update_progress("Summarising TFBM hits per LTR")
  
  ltr_split <- split(tfbm_hits, tfbm_hits$sequence_name)
  
  ltr_tfbm_summary <- do.call(
    rbind,
    lapply(names(ltr_split), function(ltr) {
      z <- ltr_split[[ltr]]
      
      data.frame(
        sequence_name = ltr,
        n_tfbm_hits = nrow(z),
        n_unique_motifs = length(unique(z$motif_id)),
        n_unique_tf_names = length(unique(z$motif_alt_id)),
        tf_names_all = collapse_unique_full(names(sort(table(z$motif_alt_id), decreasing = TRUE))),
        top_tf_names_preview = collapse_unique_preview(names(sort(table(z$motif_alt_id), decreasing = TRUE))),
        best_p_value = suppressWarnings(min(z[["p-value"]], na.rm = TRUE)),
        best_q_value = if ("q-value" %in% colnames(z)) {
          suppressWarnings(min(z[["q-value"]], na.rm = TRUE))
        } else {
          NA_real_
        },
        best_score = suppressWarnings(max(z$score, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )
    })
  )
  
  rownames(ltr_tfbm_summary) <- NULL
  
  write.table(
    ltr_tfbm_summary,
    file = file.path(output_dir, "ltr_tfbm_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  update_progress("Summarising TFBM hits per HERV")
  
  summarize_tfbm_subset <- function(z, prefix) {
    if (nrow(z) == 0) {
      out <- data.frame(
        n_tfbm_hits = 0,
        n_unique_motifs = 0,
        n_unique_tf_names = 0,
        tf_names_all = NA_character_,
        top_tf_names_preview = NA_character_,
        best_p_value = NA_real_,
        best_q_value = NA_real_,
        best_score = NA_real_,
        stringsAsFactors = FALSE
      )
    } else {
      tf_ordered <- names(sort(table(z$motif_alt_id), decreasing = TRUE))
      
      out <- data.frame(
        n_tfbm_hits = nrow(z),
        n_unique_motifs = length(unique(z$motif_id)),
        n_unique_tf_names = length(unique(z$motif_alt_id)),
        tf_names_all = collapse_unique_full(tf_ordered),
        top_tf_names_preview = collapse_unique_preview(tf_ordered),
        best_p_value = suppressWarnings(min(z[["p-value"]], na.rm = TRUE)),
        best_q_value = if ("q-value" %in% colnames(z)) {
          suppressWarnings(min(z[["q-value"]], na.rm = TRUE))
        } else {
          NA_real_
        },
        best_score = suppressWarnings(max(z$score, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )
    }
    
    colnames(out) <- paste0(prefix, "_", colnames(out))
    out
  }
  
  herv_split <- split(tfbm_hits, tfbm_hits$HERV_id)
  
  herv_tfbm_summary <- do.call(
    rbind,
    lapply(names(herv_split), function(herv) {
      z <- herv_split[[herv]]
      
      z_ltr5 <- z[z$ltr_position == "LTR5", , drop = FALSE]
      z_ltr3 <- z[z$ltr_position == "LTR3", , drop = FALSE]
      
      tf_ordered_total <- names(sort(table(z$motif_alt_id), decreasing = TRUE))
      
      total_summary <- data.frame(
        HERV_id = herv,
        n_ltrs_with_tfbm_hits = length(unique(z$sequence_name)),
        total_n_tfbm_hits = nrow(z),
        total_n_unique_motifs = length(unique(z$motif_id)),
        total_n_unique_tf_names = length(unique(z$motif_alt_id)),
        total_tf_names_all = collapse_unique_full(tf_ordered_total),
        total_top_tf_names_preview = collapse_unique_preview(tf_ordered_total),
        total_best_p_value = suppressWarnings(min(z[["p-value"]], na.rm = TRUE)),
        total_best_q_value = if ("q-value" %in% colnames(z)) {
          suppressWarnings(min(z[["q-value"]], na.rm = TRUE))
        } else {
          NA_real_
        },
        total_best_score = suppressWarnings(max(z$score, na.rm = TRUE)),
        stringsAsFactors = FALSE
      )
      
      cbind(
        total_summary,
        summarize_tfbm_subset(z_ltr5, "ltr5"),
        summarize_tfbm_subset(z_ltr3, "ltr3")
      )
    })
  )
  
  rownames(herv_tfbm_summary) <- NULL
  
  rownames(herv_tfbm_summary) <- NULL
  
  write.table(
    herv_tfbm_summary,
    file = file.path(output_dir, "herv_tfbm_summary.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  update_progress("Generating top motif table")
  
  motif_counts <- as.data.frame(sort(table(tfbm_hits$motif_alt_id), decreasing = TRUE))
  colnames(motif_counts) <- c("motif_alt_id", "n_hits")
  
  motif_herv_counts <- aggregate(
    HERV_id ~ motif_alt_id,
    data = unique(tfbm_hits[, c("motif_alt_id", "HERV_id")]),
    FUN = length
  )
  
  colnames(motif_herv_counts)[2] <- "n_hervs"
  
  top_tfbm_motifs <- merge(
    motif_counts,
    motif_herv_counts,
    by = "motif_alt_id",
    all.x = TRUE
  )
  
  top_tfbm_motifs <- top_tfbm_motifs[
    order(top_tfbm_motifs$n_hervs, top_tfbm_motifs$n_hits, decreasing = TRUE),
  ]
  
  write.table(
    top_tfbm_motifs,
    file = file.path(output_dir, "top_tfbm_motifs.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  update_progress("Adding detailed TFBM summaries to feature table")
  
  idx <- match(features$HERV_id, herv_tfbm_summary$HERV_id)
  
  features$total_n_tfbm_hits_detailed <- 0
  features$total_n_unique_tfbm_motifs_detailed <- 0
  features$total_n_unique_tfbm_tf_names_detailed <- 0
  features$total_tfbm_tf_names_all <- NA_character_
  features$total_top_tfbm_tf_names_preview <- NA_character_
  features$total_best_tfbm_p_value <- NA_real_
  features$total_best_tfbm_q_value <- NA_real_
  features$total_best_tfbm_score <- NA_real_
  
  features$ltr5_n_tfbm_hits_detailed <- 0
  features$ltr5_n_unique_tfbm_motifs_detailed <- 0
  features$ltr5_n_unique_tfbm_tf_names_detailed <- 0
  features$ltr5_tfbm_tf_names_all <- NA_character_
  features$ltr5_top_tfbm_tf_names_preview <- NA_character_
  features$ltr5_best_tfbm_p_value <- NA_real_
  features$ltr5_best_tfbm_q_value <- NA_real_
  features$ltr5_best_tfbm_score <- NA_real_
  
  features$ltr3_n_tfbm_hits_detailed <- 0
  features$ltr3_n_unique_tfbm_motifs_detailed <- 0
  features$ltr3_n_unique_tfbm_tf_names_detailed <- 0
  features$ltr3_tfbm_tf_names_all <- NA_character_
  features$ltr3_top_tfbm_tf_names_preview <- NA_character_
  features$ltr3_best_tfbm_p_value <- NA_real_
  features$ltr3_best_tfbm_q_value <- NA_real_
  features$ltr3_best_tfbm_score <- NA_real_
  
  has_match <- !is.na(idx)
  
  features$total_n_tfbm_hits_detailed[has_match] <-
    herv_tfbm_summary$total_n_tfbm_hits[idx[has_match]]
  
  features$total_n_unique_tfbm_motifs_detailed[has_match] <-
    herv_tfbm_summary$total_n_unique_motifs[idx[has_match]]
  
  features$total_n_unique_tfbm_tf_names_detailed[has_match] <-
    herv_tfbm_summary$total_n_unique_tf_names[idx[has_match]]
  
  features$total_tfbm_tf_names_all[has_match] <-
    herv_tfbm_summary$total_tf_names_all[idx[has_match]]
  
  features$total_top_tfbm_tf_names_preview[has_match] <-
    herv_tfbm_summary$total_top_tf_names_preview[idx[has_match]]
  
  features$total_best_tfbm_p_value[has_match] <-
    herv_tfbm_summary$total_best_p_value[idx[has_match]]
  
  features$total_best_tfbm_q_value[has_match] <-
    herv_tfbm_summary$total_best_q_value[idx[has_match]]
  
  features$total_best_tfbm_score[has_match] <-
    herv_tfbm_summary$total_best_score[idx[has_match]]
  
  
  features$ltr5_n_tfbm_hits_detailed[has_match] <-
    herv_tfbm_summary$ltr5_n_tfbm_hits[idx[has_match]]
  
  features$ltr5_n_unique_tfbm_motifs_detailed[has_match] <-
    herv_tfbm_summary$ltr5_n_unique_motifs[idx[has_match]]
  
  features$ltr5_n_unique_tfbm_tf_names_detailed[has_match] <-
    herv_tfbm_summary$ltr5_n_unique_tf_names[idx[has_match]]
  
  features$ltr5_tfbm_tf_names_all[has_match] <-
    herv_tfbm_summary$ltr5_tf_names_all[idx[has_match]]
  
  features$ltr5_top_tfbm_tf_names_preview[has_match] <-
    herv_tfbm_summary$ltr5_top_tf_names_preview[idx[has_match]]
  
  features$ltr5_best_tfbm_p_value[has_match] <-
    herv_tfbm_summary$ltr5_best_p_value[idx[has_match]]
  
  features$ltr5_best_tfbm_q_value[has_match] <-
    herv_tfbm_summary$ltr5_best_q_value[idx[has_match]]
  
  features$ltr5_best_tfbm_score[has_match] <-
    herv_tfbm_summary$ltr5_best_score[idx[has_match]]
  
  
  features$ltr3_n_tfbm_hits_detailed[has_match] <-
    herv_tfbm_summary$ltr3_n_tfbm_hits[idx[has_match]]
  
  features$ltr3_n_unique_tfbm_motifs_detailed[has_match] <-
    herv_tfbm_summary$ltr3_n_unique_motifs[idx[has_match]]
  
  features$ltr3_n_unique_tfbm_tf_names_detailed[has_match] <-
    herv_tfbm_summary$ltr3_n_unique_tf_names[idx[has_match]]
  
  features$ltr3_tfbm_tf_names_all[has_match] <-
    herv_tfbm_summary$ltr3_tf_names_all[idx[has_match]]
  
  features$ltr3_top_tfbm_tf_names_preview[has_match] <-
    herv_tfbm_summary$ltr3_top_tf_names_preview[idx[has_match]]
  
  features$ltr3_best_tfbm_p_value[has_match] <-
    herv_tfbm_summary$ltr3_best_p_value[idx[has_match]]
  
  features$ltr3_best_tfbm_q_value[has_match] <-
    herv_tfbm_summary$ltr3_best_q_value[idx[has_match]]
  
  features$ltr3_best_tfbm_score[has_match] <-
    herv_tfbm_summary$ltr3_best_score[idx[has_match]]
  
  update_progress("Writing final outputs")
  
  out_features <- if (.return_full_features) {
    features
  } else {
    select_compact_herv_features(features)
  }

  if (inherits(x, "HERVarium_annotation")) {
    x$features <- out_features
    x$ltr_tfbm_hits <- tfbm_hits
    x$ltr_tfbm_summary <- ltr_tfbm_summary
    x$herv_tfbm_summary <- herv_tfbm_summary
    x$top_tfbm_motifs <- top_tfbm_motifs
    
    write.table(
      select_compact_herv_features(features),
      file = file.path(output_dir, "herv_features_compact.tsv"),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    message("Detailed LTR TFBM annotations written to: ", output_dir)
    
    return(x)
  }
  
  out_features
}