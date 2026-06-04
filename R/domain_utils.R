# ------------------------------------------------------------
# HERVariumR domain utilities
# ------------------------------------------------------------

.extract_domain_coverage <- function(x) {
  # Accepts entries like:
  # "pol|RNaseH:0.750"
  # "pol|RT|RT_gammaretroviridae:0.290"
  # "pol:0.437"
  
  x <- as.character(x)
  
  if (is.na(x) || x == "" || x == ".") {
    return(NA_real_)
  }
  
  cov <- sub("^.*:", "", x)
  suppressWarnings(as.numeric(cov))
}


.filter_domain_string_by_coverage <- function(x,
                                              coverage_cutoff = 0.40) {
  
  if (is.na(x) || x == "" || x == ".") {
    return(".")
  }
  
  entries <- unlist(strsplit(as.character(x), ";", fixed = TRUE))
  entries <- entries[entries != "" & entries != "."]
  
  if (length(entries) == 0) {
    return(".")
  }
  
  cov <- vapply(entries, .extract_domain_coverage, numeric(1))
  
  keep <- !is.na(cov) & cov >= coverage_cutoff
  
  if (!any(keep)) {
    return(".")
  }
  
  paste(entries[keep], collapse = ";")
}


.domain_coverages_from_string <- function(x) {
  
  if (is.na(x) || x == "" || x == ".") {
    return(numeric(0))
  }
  
  entries <- unlist(strsplit(as.character(x), ";", fixed = TRUE))
  entries <- entries[entries != "" & entries != "."]
  
  if (length(entries) == 0) {
    return(numeric(0))
  }
  
  cov <- vapply(entries, .extract_domain_coverage, numeric(1))
  cov[!is.na(cov)]
}


filter_herv_domains_by_coverage <- function(df,
                                            coverage_cutoff = 0.40,
                                            keep_raw_domain_columns = TRUE) {
  
  if (is.null(coverage_cutoff)) {
    return(df)
  }
  
  if (!is.numeric(coverage_cutoff) || length(coverage_cutoff) != 1) {
    stop("coverage_cutoff must be a single numeric value.")
  }
  
  if (is.na(coverage_cutoff) || coverage_cutoff < 0 || coverage_cutoff > 1) {
    stop("coverage_cutoff must be between 0 and 1.")
  }
  
  domain_cols <- c(
    "domains_gene",
    "domains_type",
    "domains_profile"
  )
  
  domain_cols <- domain_cols[domain_cols %in% colnames(df)]
  
  if (length(domain_cols) == 0) {
    warning("No domain annotation columns found to filter.")
    return(df)
  }
  
  if (keep_raw_domain_columns) {
    for (col in domain_cols) {
      raw_col <- paste0(col, "_raw")
      if (!raw_col %in% colnames(df)) {
        df[[raw_col]] <- df[[col]]
      }
    }
    
    if ("domain_count" %in% colnames(df) &&
        !"domain_count_raw" %in% colnames(df)) {
      df$domain_count_raw <- df$domain_count
    }
    
    if ("max_domain_cov" %in% colnames(df) &&
        !"max_domain_cov_raw" %in% colnames(df)) {
      df$max_domain_cov_raw <- df$max_domain_cov
    }
  }
  
  for (col in domain_cols) {
    df[[col]] <- vapply(
      df[[col]],
      .filter_domain_string_by_coverage,
      character(1),
      coverage_cutoff = coverage_cutoff
    )
  }
  
  # Recompute domain_count and max_domain_cov from the filtered domains_type
  # when possible. This is the most specific domain representation.
  count_source_col <- if ("domains_type" %in% colnames(df)) {
    "domains_type"
  } else {
    domain_cols[1]
  }
  
  cov_list <- lapply(df[[count_source_col]], .domain_coverages_from_string)
  
  df$domain_count <- vapply(cov_list, length, integer(1))
  
  df$max_domain_cov <- vapply(
    cov_list,
    function(z) {
      if (length(z) == 0) {
        return(NA_real_)
      }
      max(z, na.rm = TRUE)
    },
    numeric(1)
  )
  
  df
}


apply_terminal_domain_coverage_cutoff <- function(features,
                                                  coverage_cutoff = 0.40) {
  
  if (is.null(coverage_cutoff)) {
    return(features)
  }
  
  contexts <- c("protein_coding", "lncRNA")
  
  for (ctx in contexts) {
    
    has_col <- paste0("has_terminal_exon_domain_", ctx)
    n_col <- paste0("n_terminal_domain_hits_", ctx)
    max_col <- paste0("max_terminal_domain_coverage_", ctx)
    
    if (has_col %in% colnames(features) &&
        max_col %in% colnames(features)) {
      
      valid <- !is.na(features[[max_col]]) &
        features[[max_col]] >= coverage_cutoff
      
      features[[has_col]] <- features[[has_col]] & valid
      
      if (n_col %in% colnames(features)) {
        features[[n_col]][!valid] <- 0
      }
    }
  }
  
  terminal_cols <- c(
    "has_terminal_exon_domain_protein_coding",
    "has_terminal_exon_domain_lncRNA"
  )
  
  terminal_cols <- terminal_cols[terminal_cols %in% colnames(features)]
  
  if (length(terminal_cols) > 0) {
    features$has_any_terminal_exon_domain <- Reduce(
      `|`,
      lapply(terminal_cols, function(z) features[[z]])
    )
  }
  
  features
}