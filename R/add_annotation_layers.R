add_ltr_ifn_annotations <- function(features,
                                    stat1_file = NULL,
                                    stat1stat2_irf_file = NULL) {
  
  if (is.null(stat1_file)) {
    stat1_file <- system.file(
      "extdata",
      "LTR_IFN_STAT1_summary.tsv",
      package = "HERVariumR"
    )
  }
  
  if (is.null(stat1stat2_irf_file)) {
    stat1stat2_irf_file <- system.file(
      "extdata",
      "LTR_IFN_STAT1STAT2_IRF_summary.tsv",
      package = "HERVariumR"
    )
  }
  
  add_ltr_table <- function(features, tbl, ltr_col, prefix, selected_cols) {
    idx <- match(features[[ltr_col]], tbl$sequence_name)
    
    for (col in selected_cols) {
      new_col <- paste0(prefix, "_", col)
      
      if (col %in% colnames(tbl)) {
        features[[new_col]] <- tbl[[col]][idx]
      } else {
        features[[new_col]] <- NA
      }
    }
    
    features
  }
  
  if (!is.null(stat1_file) && file.exists(stat1_file)) {
    stat1 <- read.delim(
      stat1_file,
      sep = "\t",
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    stat1_cols <- c(
      "n_stat1",
      "best_p_stat1",
      "best_score_stat1",
      "n_irf",
      "best_p_irf",
      "best_score_irf",
      "max_stat1_irf_in_75bp",
      "stat1_in_best_ifn_window",
      "irf_in_best_ifn_window",
      "ifn_stimulation_potential"
    )
    
    features <- add_ltr_table(
      features = features,
      tbl = stat1,
      ltr_col = "ltr5_name",
      prefix = "ltr5_stat1",
      selected_cols = stat1_cols
    )
    
    features <- add_ltr_table(
      features = features,
      tbl = stat1,
      ltr_col = "ltr3_name",
      prefix = "ltr3_stat1",
      selected_cols = stat1_cols
    )
  }
  
  if (!is.null(stat1stat2_irf_file) && file.exists(stat1stat2_irf_file)) {
    stat1stat2 <- read.delim(
      stat1stat2_irf_file,
      sep = "\t",
      header = TRUE,
      stringsAsFactors = FALSE,
      check.names = FALSE
    )
    
    stat1stat2_cols <- c(
      "n_stat1stat2",
      "best_p_stat1stat2",
      "best_score_stat1stat2",
      "n_irf",
      "best_p_irf",
      "best_score_irf",
      "max_stat1stat2_irf_in_75bp",
      "stat1stat2_in_best_ifn_window",
      "irf_in_best_ifn_window",
      "ifn_stimulation_potential"
    )
    
    features <- add_ltr_table(
      features = features,
      tbl = stat1stat2,
      ltr_col = "ltr5_name",
      prefix = "ltr5_stat1stat2_irf",
      selected_cols = stat1stat2_cols
    )
    
    features <- add_ltr_table(
      features = features,
      tbl = stat1stat2,
      ltr_col = "ltr3_name",
      prefix = "ltr3_stat1stat2_irf",
      selected_cols = stat1stat2_cols
    )
  }
  
  features$has_ltr5_stat1_motif <-
    !is.na(features$ltr5_stat1_n_stat1) &
    features$ltr5_stat1_n_stat1 > 0
  
  features$has_ltr3_stat1_motif <-
    !is.na(features$ltr3_stat1_n_stat1) &
    features$ltr3_stat1_n_stat1 > 0
  
  features$has_ltr5_stat1stat2_motif <-
    !is.na(features$ltr5_stat1stat2_irf_n_stat1stat2) &
    features$ltr5_stat1stat2_irf_n_stat1stat2 > 0
  
  features$has_ltr3_stat1stat2_motif <-
    !is.na(features$ltr3_stat1stat2_irf_n_stat1stat2) &
    features$ltr3_stat1stat2_irf_n_stat1stat2 > 0
  
  features$has_any_stat1_motif <-
    features$has_ltr5_stat1_motif |
    features$has_ltr3_stat1_motif
  
  features$has_any_stat1stat2_motif <-
    features$has_ltr5_stat1stat2_motif |
    features$has_ltr3_stat1stat2_motif
  
  features$has_any_ifn_related_ltr_motif <-
    features$has_any_stat1_motif |
    features$has_any_stat1stat2_motif
  
  features
}


add_terminal_exon_domain_annotations <- function(features,
                                                 last_exon_file = NULL) {
  
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is required. Install it with install.packages('readxl').")
  }
  
  if (is.null(last_exon_file)) {
    last_exon_file <- system.file(
      "extdata",
      "HERV_domains_transcript_context_last_exon.xlsx",
      package = "HERVariumR"
    )
  }
  
  if (is.null(last_exon_file) || !file.exists(last_exon_file)) {
    warning("Last exon annotation file not found. Skipping terminal exon layer.")
    return(features)
  }
  
  pc <- readxl::read_excel(last_exon_file, sheet = "protein_coding")
  lnc <- readxl::read_excel(last_exon_file, sheet = "lncRNA")
  
  pc$gene_type_layer <- "protein_coding"
  lnc$gene_type_layer <- "lncRNA"
  
  all_terminal <- rbind(
    as.data.frame(pc, stringsAsFactors = FALSE),
    as.data.frame(lnc, stringsAsFactors = FALSE)
  )
  
  is_terminal <- all_terminal$any_terminal_exon_overlap %in% c(
    TRUE, "TRUE", "True", "true", 1, "1"
  )
  
  all_terminal <- all_terminal[is_terminal, ]
  
  collapse_unique <- function(x) {
    x <- unique(na.omit(as.character(x)))
    x <- x[x != "" & x != "."]
    if (length(x) == 0) return(NA_character_)
    paste(x, collapse = ";")
  }
  
  summarize_layer <- function(df, layer_name) {
    df <- df[df$gene_type_layer == layer_name, ]
    
    if (nrow(df) == 0) {
      return(data.frame())
    }
    
    split_df <- split(df, df$Locus)
    
    out <- do.call(
      rbind,
      lapply(names(split_df), function(locus) {
        z <- split_df[[locus]]
        
        data.frame(
          locid = locus,
          n_terminal_domain_hits = nrow(z),
          terminal_domain_types = collapse_unique(z$DomainType),
          terminal_domain_genes = collapse_unique(z$Gene),
          terminal_domain_names = collapse_unique(z$Domain),
          terminal_domain_overlapping_genes = collapse_unique(z$gene_name),
          max_terminal_domain_coverage = suppressWarnings(max(z$Coverage_num, na.rm = TRUE)),
          stringsAsFactors = FALSE
        )
      })
    )
    
    rownames(out) <- NULL
    out
  }
  
  pc_summary <- summarize_layer(all_terminal, "protein_coding")
  lnc_summary <- summarize_layer(all_terminal, "lncRNA")
  
  add_summary <- function(features, summary_df, prefix) {
    features[[paste0("has_terminal_exon_domain_", prefix)]] <- FALSE
    features[[paste0("n_terminal_domain_hits_", prefix)]] <- 0
    features[[paste0("terminal_domain_types_", prefix)]] <- NA_character_
    features[[paste0("terminal_domain_genes_", prefix)]] <- NA_character_
    features[[paste0("terminal_domain_names_", prefix)]] <- NA_character_
    features[[paste0("terminal_domain_overlapping_genes_", prefix)]] <- NA_character_
    features[[paste0("max_terminal_domain_coverage_", prefix)]] <- NA_real_
    
    if (nrow(summary_df) == 0) {
      return(features)
    }
    
    idx <- match(features$locid, summary_df$locid)
    has_match <- !is.na(idx)
    
    features[[paste0("has_terminal_exon_domain_", prefix)]][has_match] <- TRUE
    features[[paste0("n_terminal_domain_hits_", prefix)]][has_match] <-
      summary_df$n_terminal_domain_hits[idx[has_match]]
    features[[paste0("terminal_domain_types_", prefix)]][has_match] <-
      summary_df$terminal_domain_types[idx[has_match]]
    features[[paste0("terminal_domain_genes_", prefix)]][has_match] <-
      summary_df$terminal_domain_genes[idx[has_match]]
    features[[paste0("terminal_domain_names_", prefix)]][has_match] <-
      summary_df$terminal_domain_names[idx[has_match]]
    features[[paste0("terminal_domain_overlapping_genes_", prefix)]][has_match] <-
      summary_df$terminal_domain_overlapping_genes[idx[has_match]]
    features[[paste0("max_terminal_domain_coverage_", prefix)]][has_match] <-
      summary_df$max_terminal_domain_coverage[idx[has_match]]
    
    features
  }
  
  features <- add_summary(features, pc_summary, "protein_coding")
  features <- add_summary(features, lnc_summary, "lncRNA")
  
  features$has_any_terminal_exon_domain <-
    features$has_terminal_exon_domain_protein_coding |
    features$has_terminal_exon_domain_lncRNA
  
  features
}