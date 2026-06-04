plot_internal_domain_bubbles <- function(x,
                                         output_dir = NULL,
                                         top_n_subfamilies = 10,
                                         top_n_domains = 25,
                                         jitter_height = 0.15,
                                         point_size = 3,
                                         alpha = 0.75,
                                         size_by = c("fixed", "domain_count")) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  
  # ------------------------------------------------------------
  # 1. Extract features table
  # ------------------------------------------------------------
  
  if (inherits(x, "HERVarium_annotation")) {
    features <- x$features
    
    if (is.null(output_dir)) {
      output_dir <- x$output_dir
    }
  } else {
    features <- x
  }
  
  if (is.null(output_dir)) {
    stop("Please provide output_dir or use a HERVarium_annotation object with output_dir.")
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  required_cols <- c("HERV_id", "subfamily", "domains_type")
  missing_cols <- setdiff(required_cols, colnames(features))
  
  if (length(missing_cols) > 0) {
    stop(
      "The following required columns are missing from features: ",
      paste(missing_cols, collapse = ", ")
    )
  }
  
  # ------------------------------------------------------------
  # 2. Parse domains_type
  # Expected format:
  # gag|GAG:0.144;pol|RNaseH:1.000;env|ENV:0.212
  # ------------------------------------------------------------
  
  parse_domains_type <- function(z, HERV_id, subfamily, domain_count) {
    
    if (is.na(z) || z == "" || z == ".") {
      return(NULL)
    }
    
    entries <- unlist(strsplit(as.character(z), ";", fixed = TRUE))
    entries <- entries[entries != "" & entries != "."]
    
    if (length(entries) == 0) {
      return(NULL)
    }
    
    out <- lapply(entries, function(e) {
      
      # Split gene class from the rest
      parts1 <- unlist(strsplit(e, "|", fixed = TRUE))
      
      if (length(parts1) != 2) {
        return(NULL)
      }
      
      gene_class <- parts1[1]
      domain_and_cov <- parts1[2]
      
      # Split domain type from coverage
      parts2 <- unlist(strsplit(domain_and_cov, ":", fixed = TRUE))
      
      if (length(parts2) != 2) {
        return(NULL)
      }
      
      domain_type <- parts2[1]
      coverage <- suppressWarnings(as.numeric(parts2[2]))
      
      if (is.na(coverage)) {
        return(NULL)
      }
      
      data.frame(
        HERV_id = HERV_id,
        subfamily = subfamily,
        domain_count = domain_count,
        gene_class = gene_class,
        domain_type = domain_type,
        coverage = coverage,
        raw_domain_entry = e,
        stringsAsFactors = FALSE
      )
    })
    
    do.call(rbind, out)
  }
  
  domain_hits <- do.call(
    rbind,
    lapply(seq_len(nrow(features)), function(i) {
      parse_domains_type(
        z = features$domains_type[i],
        HERV_id = features$HERV_id[i],
        subfamily = features$subfamily[i],
        domain_count = features$domain_count[i]
      )
    })
  )
  
  if (is.null(domain_hits) || nrow(domain_hits) == 0) {
    stop("No valid internal domain hits could be parsed from features$domains_type.")
  }
  
  # ------------------------------------------------------------
  # 3. Clean labels
  # ------------------------------------------------------------
  
  domain_hits$domain_type <- as.character(domain_hits$domain_type)
  
  domain_hits$domain_type[domain_hits$domain_type == "ENV"] <- "Env"
  domain_hits$domain_type[domain_hits$domain_type == "GAG"] <- "Gag"
  domain_hits$domain_type[domain_hits$domain_type == "AP"] <- "Protease"
  domain_hits$domain_type[domain_hits$domain_type == "DUT"] <- "dUTPase"
  domain_hits$domain_type[domain_hits$domain_type == "INT"] <- "Integrase"
  domain_hits$domain_type[domain_hits$domain_type == "RT"] <- "RT"
  domain_hits$domain_type[domain_hits$domain_type == "RH"] <- "RNaseH"
  
  domain_hits$subfamily[
    is.na(domain_hits$subfamily) |
      domain_hits$subfamily == "" |
      domain_hits$subfamily == "."
  ] <- "Unknown"
  
  # ------------------------------------------------------------
  # 4. Collapse subfamilies to top N + Other
  # ------------------------------------------------------------
  
  subfamily_info <- unique(domain_hits[, c("HERV_id", "subfamily")])
  subfamily_counts <- as.data.frame(table(subfamily_info$subfamily))
  colnames(subfamily_counts) <- c("subfamily", "n_hervs")
  subfamily_counts <- subfamily_counts[order(subfamily_counts$n_hervs, decreasing = TRUE), ]
  
  main_subfamilies <- head(subfamily_counts$subfamily, top_n_subfamilies)
  
  domain_hits$subfamily_group <- ifelse(
    domain_hits$subfamily %in% main_subfamilies,
    domain_hits$subfamily,
    "Other"
  )
  
  # ------------------------------------------------------------
  # 5. Restrict to top domain types
  # ------------------------------------------------------------
  
  domain_counts <- as.data.frame(table(domain_hits$domain_type))
  colnames(domain_counts) <- c("domain_type", "n_hits")
  domain_counts <- domain_counts[order(domain_counts$n_hits, decreasing = TRUE), ]
  
  selected_domains <- head(domain_counts$domain_type, top_n_domains)
  domain_hits <- domain_hits[domain_hits$domain_type %in% selected_domains, ]
  
  # Preferred biological order
  preferred_order <- c(
    "Gag",
    "Protease",
    "RT",
    "RNaseH",
    "Integrase",
    "dUTPase",
    "Chromodomain",
    "Env",
    "Accessory"
  )
  
  observed_order <- unique(domain_hits$domain_type)
  final_order <- c(
    preferred_order[preferred_order %in% observed_order],
    setdiff(observed_order, preferred_order)
  )
  
  domain_hits$domain_type <- factor(
    domain_hits$domain_type,
    levels = rev(final_order)
  )
  
  # ------------------------------------------------------------
  # 6. Plot: x = HMM coverage, y = domain type
  # ------------------------------------------------------------
  
  if (size_by == "fixed") {
    
    p <- ggplot2::ggplot(
      domain_hits,
      ggplot2::aes(
        x = coverage,
        y = domain_type,
        color = subfamily_group
      )
    ) +
      ggplot2::geom_jitter(
        height = jitter_height,
        width = 0,
        size = point_size,
        alpha = alpha
      ) +
      ggplot2::labs(
        size = NULL
      )
    
  } else if (size_by == "domain_count") {
    
    p <- ggplot2::ggplot(
      domain_hits,
      ggplot2::aes(
        x = coverage,
        y = domain_type,
        color = subfamily_group,
        size = domain_count
      )
    ) +
      ggplot2::geom_jitter(
        height = jitter_height,
        width = 0,
        alpha = alpha
      ) +
      ggplot2::scale_size_continuous(
        range = c(1.5, 6)
      ) +
      ggplot2::labs(
        size = "Domains per HERV"
      )
  }
  
  p <- p +
    ggplot2::scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.2)
    ) +
    ggplot2::labs(
      title = "Internal domain conservation landscape",
      subtitle = "Each bubble is one internal-domain hit; x-axis shows HMM profile coverage",
      x = "HMM profile coverage",
      y = "Domain type",
      color = "Subfamily"
    ) +
    ggplot2::theme_bw()
  
  ggplot2::ggsave(
    filename = file.path(output_dir, "internal_domain_conservation_bubbleplot.png"),
    plot = p,
    width = 10,
    height = 6,
    dpi = 300
  )
  
  # ------------------------------------------------------------
  # 7. Write plot-supporting table
  # ------------------------------------------------------------
  
  write.table(
    domain_hits,
    file = file.path(output_dir, "plot_table_internal_domain_conservation_bubbleplot.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  message("Internal domain conservation bubble plot written to: ", output_dir)
  
  invisible(p)
}