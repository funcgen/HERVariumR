plot_herv_annotation <- function(x, output_dir = NULL, top_n = 25) {
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  
  if (inherits(x, "HERVarium_annotation")) {
    features <- x$features
  } else {
    features <- x
  }
  
  if (is.null(output_dir)) output_dir <- x$output_dir
  if (is.null(output_dir)) stop("Please provide output_dir.")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  
  # 1. Domain class presence
  domain_long <- data.frame(
    domain = rep(c("gag", "pol", "env", "accessory"), each = nrow(features)),
    present = c(features$has_gag, features$has_pol, features$has_env, features$has_accessory)
  )
  
  domain_counts <- aggregate(present ~ domain, domain_long, sum)
  
  p <- ggplot2::ggplot(domain_counts, ggplot2::aes(x = reorder(domain, present), y = present)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Retroviral domain class presence",
      x = NULL,
      y = "Number of HERVs"
    ) +
    ggplot2::theme_bw()
  
  ggplot2::ggsave(file.path(output_dir, "domain_class_presence.png"), p, width = 7, height = 4, dpi = 300)
  
  # 1b. Specific internal domain type presence
  # Examples: Gag, Protease, RT, RNaseH, Integrase, Env, dUTPase
  
  parse_domain_types <- function(domains_type) {
    if (is.na(domains_type) || domains_type == "" || domains_type == ".") {
      return(character(0))
    }
    
    entries <- unlist(strsplit(as.character(domains_type), ";", fixed = TRUE))
    entries <- entries[entries != "" & entries != "."]
    
    domain_types <- vapply(entries, function(e) {
      parts1 <- unlist(strsplit(e, "|", fixed = TRUE))
      if (length(parts1) != 2) return(NA_character_)
      
      domain_and_cov <- parts1[2]
      parts2 <- unlist(strsplit(domain_and_cov, ":", fixed = TRUE))
      if (length(parts2) != 2) return(NA_character_)
      
      parts2[1]
    }, character(1))
    
    domain_types <- domain_types[!is.na(domain_types)]
    unique(domain_types)
  }
  
  clean_domain_type <- function(z) {
    z[z == "ENV"] <- "Env"
    z[z == "GAG"] <- "Gag"
    z[z == "AP"] <- "Protease"
    z[z == "DUT"] <- "dUTPase"
    z[z == "INT"] <- "Integrase"
    z[z == "RT"] <- "RT"
    z[z == "RNaseH"] <- "RNaseH"
    z
  }
  
  if ("domains_type" %in% colnames(features)) {
    
    domain_type_long <- do.call(
      rbind,
      lapply(seq_len(nrow(features)), function(i) {
        domain_types <- parse_domain_types(features$domains_type[i])
        
        if (length(domain_types) == 0) {
          return(NULL)
        }
        
        data.frame(
          HERV_id = features$HERV_id[i],
          domain_type = clean_domain_type(domain_types),
          stringsAsFactors = FALSE
        )
      })
    )
    
    if (!is.null(domain_type_long) && nrow(domain_type_long) > 0) {
      
      # Presence per HERV, not raw hit count
      domain_type_long <- unique(domain_type_long)
      
      domain_type_counts <- as.data.frame(table(domain_type_long$domain_type))
      colnames(domain_type_counts) <- c("domain_type", "n_hervs")
      domain_type_counts <- domain_type_counts[
        order(domain_type_counts$n_hervs, decreasing = TRUE),
      ]
      
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
      
      observed <- domain_type_counts$domain_type
      final_order <- c(
        preferred_order[preferred_order %in% observed],
        setdiff(observed, preferred_order)
      )
      
      domain_type_counts$domain_type <- factor(
        domain_type_counts$domain_type,
        levels = rev(final_order)
      )
      
      p <- ggplot2::ggplot(
        domain_type_counts,
        ggplot2::aes(x = domain_type, y = n_hervs)
      ) +
        ggplot2::geom_col() +
        ggplot2::coord_flip() +
        ggplot2::labs(
          title = "Specific internal domain type presence",
          x = NULL,
          y = "Number of HERVs"
        ) +
        ggplot2::theme_bw()
      
      ggplot2::ggsave(
        file.path(output_dir, "domain_type_presence.png"),
        p,
        width = 7,
        height = 4.5,
        dpi = 300
      )
      
      write.table(
        domain_type_counts,
        file = file.path(output_dir, "plot_table_domain_type_presence.tsv"),
        sep = "\t",
        quote = FALSE,
        row.names = FALSE
      )
    }
  }
  
  
  
  
  
  # 2. Domain count distribution
  p <- ggplot2::ggplot(features, ggplot2::aes(x = domain_count)) +
    ggplot2::geom_histogram(binwidth = 1, boundary = -0.5) +
    ggplot2::labs(
      title = "Domain count distribution",
      x = "Number of detected domains",
      y = "Number of HERVs"
    ) +
    ggplot2::theme_bw()
  
  ggplot2::ggsave(file.path(output_dir, "domain_count_distribution.png"), p, width = 7, height = 4, dpi = 300)
  
  # 3. Max domain coverage distribution
  p <- ggplot2::ggplot(features, ggplot2::aes(x = max_domain_cov)) +
    ggplot2::geom_histogram(bins = 30) +
    ggplot2::labs(
      title = "Maximum domain coverage distribution",
      x = "Maximum HMM profile coverage",
      y = "Number of HERVs"
    ) +
    ggplot2::theme_bw()
  
  ggplot2::ggsave(file.path(output_dir, "max_domain_coverage_distribution.png"), p, width = 7, height = 4, dpi = 300)
  
  # 4. LTR TFBM burden distribution
  # The compact feature contract retains the separate 5′ and 3′ LTR burdens.
  # Reconstruct the historical combined value locally for this plot only.
  if (!"total_ltr_tfbm_burden" %in% colnames(features)) {
    required_ltr_burden_cols <- c(
      "ltr5_tfbm_burden",
      "ltr3_tfbm_burden"
    )
    missing_ltr_burden_cols <- setdiff(
      required_ltr_burden_cols,
      colnames(features)
    )

    if (length(missing_ltr_burden_cols) > 0) {
      stop(
        "plot_herv_annotation() requires: ",
        paste(required_ltr_burden_cols, collapse = ", ")
      )
    }

    features$total_ltr_tfbm_burden <-
      suppressWarnings(as.numeric(features$ltr5_tfbm_burden)) +
      suppressWarnings(as.numeric(features$ltr3_tfbm_burden))
  }

  p <- ggplot2::ggplot(features, ggplot2::aes(x = total_ltr_tfbm_burden)) +
    ggplot2::geom_histogram(bins = 30) +
    ggplot2::labs(
      title = "LTR TFBM burden distribution",
      x = "Total LTR TFBM burden",
      y = "Number of HERVs"
    ) +
    ggplot2::theme_bw()
  
  ggplot2::ggsave(file.path(output_dir, "ltr_tfbm_burden_distribution.png"), p, width = 7, height = 4, dpi = 300)
  
  # 5. RBP burden distribution
  p <- ggplot2::ggplot(features, ggplot2::aes(x = rbp_burden)) +
    ggplot2::geom_histogram(bins = 30) +
    ggplot2::labs(
      title = "RBP burden distribution",
      x = "RBP motif burden",
      y = "Number of HERVs"
    ) +
    ggplot2::theme_bw()
  
  ggplot2::ggsave(file.path(output_dir, "rbp_burden_distribution.png"), p, width = 7, height = 4, dpi = 300)
  
  # 6. Top subfamilies
  subfamily_counts <- as.data.frame(table(features$subfamily))
  colnames(subfamily_counts) <- c("subfamily", "n")
  subfamily_counts <- subfamily_counts[order(subfamily_counts$n, decreasing = TRUE), ]
  subfamily_counts_top <- head(subfamily_counts, top_n)
  
  p <- ggplot2::ggplot(subfamily_counts_top, ggplot2::aes(x = reorder(subfamily, n), y = n)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste0("Top ", top_n, " HERV subfamilies"),
      x = NULL,
      y = "Number of HERVs"
    ) +
    ggplot2::theme_bw()
  
  ggplot2::ggsave(file.path(output_dir, "top_subfamilies.png"), p, width = 7, height = 6, dpi = 300)
  
  # 7. Genomic/transcriptomic context
  context_counts <- data.frame(
    context = c(
      "intergenic",
      "exon",
      "intron",
      "CDS",
      "UTR",
      "lncRNA",
      "protein_coding"
    ),
    n = c(
      sum(features$is_intergenic, na.rm = TRUE),
      sum(features$overlaps_exon, na.rm = TRUE),
      sum(features$overlaps_intron, na.rm = TRUE),
      sum(features$overlaps_cds, na.rm = TRUE),
      sum(features$overlaps_utr, na.rm = TRUE),
      sum(features$overlaps_lncRNA, na.rm = TRUE),
      sum(features$overlaps_protein_coding, na.rm = TRUE)
    )
  )
  
  p <- ggplot2::ggplot(context_counts, ggplot2::aes(x = reorder(context, n), y = n)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Genomic/transcriptomic context",
      x = NULL,
      y = "Number of HERVs"
    ) +
    ggplot2::theme_bw()
  
  ggplot2::ggsave(file.path(output_dir, "genomic_context_composition.png"), p, width = 7, height = 4, dpi = 300)
  
  # 8. IFN/LTR motif summary
  ifn_counts <- data.frame(
    feature = c(
      "Any IFN-related LTR motif",
      "Any STAT1 motif",
      "Any STAT1::STAT2 motif",
      "5' LTR STAT1 motif",
      "3' LTR STAT1 motif",
      "5' LTR STAT1::STAT2 motif",
      "3' LTR STAT1::STAT2 motif"
    ),
    n = c(
      sum(features$has_any_ifn_related_ltr_motif, na.rm = TRUE),
      sum(features$has_any_stat1_motif, na.rm = TRUE),
      sum(features$has_any_stat1stat2_motif, na.rm = TRUE),
      sum(features$has_ltr5_stat1_motif, na.rm = TRUE),
      sum(features$has_ltr3_stat1_motif, na.rm = TRUE),
      sum(features$has_ltr5_stat1stat2_motif, na.rm = TRUE),
      sum(features$has_ltr3_stat1stat2_motif, na.rm = TRUE)
    )
  )
  
  p <- ggplot2::ggplot(ifn_counts, ggplot2::aes(x = reorder(feature, n), y = n)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "IFN-related LTR motif summary",
      x = NULL,
      y = "Number of HERVs"
    ) +
    ggplot2::theme_bw()
  
  ggplot2::ggsave(file.path(output_dir, "ifn_ltr_motif_summary.png"), p, width = 8, height = 5, dpi = 300)
  
  # 9. Terminal exon conserved-domain summary
  terminal_counts <- data.frame(
    feature = c(
      "Any terminal exon conserved domain",
      "Protein-coding terminal exon",
      "lncRNA terminal exon"
    ),
    n = c(
      sum(features$has_any_terminal_exon_domain, na.rm = TRUE),
      sum(features$has_terminal_exon_domain_protein_coding, na.rm = TRUE),
      sum(features$has_terminal_exon_domain_lncRNA, na.rm = TRUE)
    )
  )
  
  p <- ggplot2::ggplot(terminal_counts, ggplot2::aes(x = reorder(feature, n), y = n)) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Terminal exon conserved-domain summary",
      x = NULL,
      y = "Number of HERVs"
    ) +
    ggplot2::theme_bw()
  
  ggplot2::ggsave(file.path(output_dir, "terminal_exon_domain_summary.png"), p, width = 8, height = 4, dpi = 300)
  
  # 10. IFN-related motifs by subfamily
  if ("has_any_ifn_related_ltr_motif" %in% colnames(features)) {
    sub_ifn <- aggregate(
      has_any_ifn_related_ltr_motif ~ subfamily,
      data = features,
      FUN = sum
    )
    
    sub_n <- as.data.frame(table(features$subfamily))
    colnames(sub_n) <- c("subfamily", "n_total")
    
    sub_ifn <- merge(sub_ifn, sub_n, by = "subfamily")
    sub_ifn$percent_ifn <- 100 * sub_ifn$has_any_ifn_related_ltr_motif / sub_ifn$n_total
    sub_ifn <- sub_ifn[order(sub_ifn$has_any_ifn_related_ltr_motif, decreasing = TRUE), ]
    sub_ifn_top <- head(sub_ifn, top_n)
    
    p <- ggplot2::ggplot(
      sub_ifn_top,
      ggplot2::aes(x = reorder(subfamily, has_any_ifn_related_ltr_motif),
                   y = has_any_ifn_related_ltr_motif)
    ) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = paste0("Top ", top_n, " subfamilies with IFN-related LTR motifs"),
        x = NULL,
        y = "Number of HERVs"
      ) +
      ggplot2::theme_bw()
    
    ggplot2::ggsave(file.path(output_dir, "subfamily_ifn_ltr_motif_summary.png"), p, width = 8, height = 6, dpi = 300)
  }
  
  # 11. Terminal exon domains by subfamily
  if ("has_any_terminal_exon_domain" %in% colnames(features)) {
    sub_terminal <- aggregate(
      has_any_terminal_exon_domain ~ subfamily,
      data = features,
      FUN = sum
    )
    
    sub_n <- as.data.frame(table(features$subfamily))
    colnames(sub_n) <- c("subfamily", "n_total")
    
    sub_terminal <- merge(sub_terminal, sub_n, by = "subfamily")
    sub_terminal$percent_terminal <- 100 * sub_terminal$has_any_terminal_exon_domain / sub_terminal$n_total
    sub_terminal <- sub_terminal[order(sub_terminal$has_any_terminal_exon_domain, decreasing = TRUE), ]
    sub_terminal_top <- head(sub_terminal, top_n)
    
    p <- ggplot2::ggplot(
      sub_terminal_top,
      ggplot2::aes(x = reorder(subfamily, has_any_terminal_exon_domain),
                   y = has_any_terminal_exon_domain)
    ) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::labs(
        title = paste0("Top ", top_n, " subfamilies with terminal-exon domains"),
        x = NULL,
        y = "Number of HERVs"
      ) +
      ggplot2::theme_bw()
    
    ggplot2::ggsave(file.path(output_dir, "subfamily_terminal_exon_domain_summary.png"), p, width = 8, height = 6, dpi = 300)
  }
  
  message("Summary plots written to: ", output_dir)
  invisible(output_dir)
}