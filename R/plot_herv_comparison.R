plot_herv_comparison <- function(x,
                                 output_dir = NULL,
                                 top_n = 20,
                                 plot_all_binary_features = TRUE,
                                 padj_cutoff = 0.05,
                                 show_only_significant = FALSE,
                                 pseudocount = 0.5) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  
  if (!inherits(x, "HERVarium_comparison")) {
    stop("x must be a HERVarium_comparison object returned by compare_herv_lists().")
  }
  
  if (is.null(output_dir)) {
    output_dir <- x$output_dir
  }
  
  if (is.null(output_dir)) {
    stop("Please provide output_dir or use an object with output_dir.")
  }
  
  plot_dir <- file.path(output_dir, "plots")
  
  if (!dir.exists(plot_dir)) {
    dir.create(plot_dir, recursive = TRUE)
  }
  
  safe_neglog10 <- function(p) {
    p <- suppressWarnings(as.numeric(p))
    p[p == 0] <- .Machine$double.xmin
    -log10(p)
  }
  
  add_significance_label <- function(df) {
    df$significant <- !is.na(df$padj) & df$padj < padj_cutoff
    df$significance_label <- ifelse(df$significant, paste0("FDR < ", padj_cutoff), "NS")
    df
  }
  
  # ------------------------------------------------------------
  # Helper for enrichment tables
  # ------------------------------------------------------------
  
  plot_enrichment_table <- function(df,
                                    filename,
                                    title,
                                    feature_label = "feature",
                                    max_items = top_n,
                                    pseudocount = 0.5) {
    
    if (is.null(df) || nrow(df) == 0) {
      return(invisible(NULL))
    }
    
    df <- df[!is.na(df$p_value), ]
    
    if (!("padj" %in% colnames(df))) {
      df$padj <- p.adjust(df$p_value, method = "BH")
    }
    
    df <- add_significance_label(df)
    
    if (show_only_significant) {
      df <- df[df$significant, ]
    }
    
    if (nrow(df) == 0) {
      return(invisible(NULL))
    }
    
    # Raw effect size from Fisher odds ratio
    df$log2_or_raw <- suppressWarnings(log2(df$odds_ratio))
    
    # Plotting effect size with Haldane-Anscombe correction.
    # This avoids Inf/-Inf when one cell is zero.
    required_count_cols <- c(
      "foreground_yes",
      "foreground_no",
      "background_yes",
      "background_no"
    )
    
    if (all(required_count_cols %in% colnames(df))) {
      a <- df$foreground_yes + pseudocount
      b <- df$foreground_no + pseudocount
      c <- df$background_yes + pseudocount
      d <- df$background_no + pseudocount
      
      df$odds_ratio_plot <- (a / b) / (c / d)
      df$log2_or_plot <- log2(df$odds_ratio_plot)
    } else {
      df$log2_or_plot <- df$log2_or_raw
      df$log2_or_plot[is.infinite(df$log2_or_plot)] <- NA_real_
    }
    
    df <- df[order(df$p_value), ]
    
    # If max_items is Inf, NULL, or larger than the table, keep all rows.
    if (!is.null(max_items) && is.finite(max_items)) {
      df <- head(df, max_items)
    }
    
    df$feature <- factor(df$feature, levels = rev(df$feature))
    
    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = feature,
        y = log2_or_plot,
        fill = significance_label
      )
    ) +
      ggplot2::geom_col() +
      ggplot2::coord_flip() +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
      ggplot2::labs(
        title = title,
        subtitle = paste0(
          x$foreground_name, " vs ", x$background_name,
          " | plotted log2 OR uses pseudocount = ", pseudocount
        ),
        x = NULL,
        y = "log2 odds ratio",
        fill = NULL
      ) +
      ggplot2::theme_bw()
    
    ggplot2::ggsave(
      file.path(plot_dir, filename),
      p,
      width = 8,
      height = max(5, 0.28 * nrow(df) + 2),
      dpi = 300
    )
    
    # Write plot table so the user can inspect exactly what was plotted
    table_filename <- sub("\\.png$", ".tsv", filename)
    
    write.table(
      df,
      file = file.path(plot_dir, table_filename),
      sep = "\t",
      quote = FALSE,
      row.names = FALSE
    )
    
    invisible(p)
  }
  
  
  # ------------------------------------------------------------
  # 1. Binary feature enrichment
  # ------------------------------------------------------------
  
  plot_enrichment_table(
    x$binary_features,
    filename = "binary_feature_enrichment.png",
    title = "Binary feature enrichment",
    max_items = if (plot_all_binary_features) Inf else top_n,
    pseudocount = pseudocount
  )
  
  # ------------------------------------------------------------
  # 2. Numeric feature shifts
  # ------------------------------------------------------------
  
  numeric_results <- x$numeric_features
  
  if (!is.null(numeric_results) && nrow(numeric_results) > 0) {
    
    df <- numeric_results
    df <- df[!is.na(df$p_value), ]
    df$neglog10_p <- safe_neglog10(df$p_value)
    df <- add_significance_label(df)
    
    if (show_only_significant) {
      df <- df[df$significant, ]
    }
    
    if (nrow(df) > 0) {
      df <- df[order(df$p_value), ]
      df <- head(df, top_n)
      df$feature <- factor(df$feature, levels = rev(df$feature))
      
      p <- ggplot2::ggplot(
        df,
        ggplot2::aes(
          x = feature,
          y = delta_median,
          fill = significance_label
        )
      ) +
        ggplot2::geom_col() +
        ggplot2::coord_flip() +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2::labs(
          title = "Numeric feature shifts",
          subtitle = paste0("Positive values indicate higher median in ", x$foreground_name),
          x = NULL,
          y = "Delta median",
          fill = NULL
        ) +
        ggplot2::theme_bw()
      
      ggplot2::ggsave(
        file.path(plot_dir, "numeric_feature_shifts.png"),
        p,
        width = 8,
        height = 5,
        dpi = 300
      )
    }
  }
  

  # ------------------------------------------------------------
  # 3. Subfamily enrichment
  # ------------------------------------------------------------
  
  plot_enrichment_table(
    x$subfamily_enrichment,
    filename = "subfamily_enrichment.png",
    title = "HERV subfamily enrichment"
  )
  
  # ------------------------------------------------------------
  # 4. Domain type enrichment
  # ------------------------------------------------------------
  
  plot_enrichment_table(
    x$domain_type_enrichment,
    filename = "domain_type_enrichment.png",
    title = "Internal domain type enrichment"
  )
  
  # ------------------------------------------------------------
  # 5. IFN/LTR curated features
  # ------------------------------------------------------------
  
  plot_enrichment_table(
    x$ifn_ltr_features,
    filename = "ifn_ltr_feature_enrichment.png",
    title = "IFN-related LTR feature enrichment"
  )
  
  # ------------------------------------------------------------
  # 6. LTR TFBM burden
  # ------------------------------------------------------------
  
  ltr_burden <- x$ltr_tfbm_burden
  
  if (!is.null(ltr_burden) && nrow(ltr_burden) > 0) {
    
    df <- ltr_burden
    df <- df[!is.na(df$p_value), ]
    df <- add_significance_label(df)
    
    if (nrow(df) > 0) {
      df$feature <- factor(df$feature, levels = rev(df$feature))
      
      p <- ggplot2::ggplot(
        df,
        ggplot2::aes(
          x = feature,
          y = delta_median,
          fill = significance_label
        )
      ) +
        ggplot2::geom_col() +
        ggplot2::coord_flip() +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2::labs(
          title = "LTR TFBM burden comparison",
          subtitle = paste0("Positive values indicate higher median in ", x$foreground_name),
          x = NULL,
          y = "Delta median",
          fill = NULL
        ) +
        ggplot2::theme_bw()
      
      ggplot2::ggsave(
        file.path(plot_dir, "ltr_tfbm_burden_comparison.png"),
        p,
        width = 7,
        height = 4,
        dpi = 300
      )
    }
  }
  
  # ------------------------------------------------------------
  # 7. RBP burden
  # ------------------------------------------------------------
  
  rbp_burden <- x$rbp_burden
  
  if (!is.null(rbp_burden) && nrow(rbp_burden) > 0) {
    
    df <- rbp_burden
    df <- df[!is.na(df$p_value), ]
    df <- add_significance_label(df)
    
    if (nrow(df) > 0) {
      df$feature <- factor(df$feature, levels = rev(df$feature))
      
      p <- ggplot2::ggplot(
        df,
        ggplot2::aes(
          x = feature,
          y = delta_median,
          fill = significance_label
        )
      ) +
        ggplot2::geom_col() +
        ggplot2::coord_flip() +
        ggplot2::geom_hline(yintercept = 0, linetype = "dashed") +
        ggplot2::labs(
          title = "RBP burden comparison",
          subtitle = paste0("Positive values indicate higher median in ", x$foreground_name),
          x = NULL,
          y = "Delta median",
          fill = NULL
        ) +
        ggplot2::theme_bw()
      
      ggplot2::ggsave(
        file.path(plot_dir, "rbp_burden_comparison.png"),
        p,
        width = 7,
        height = 4,
        dpi = 300
      )
    }
  }
  
  # ------------------------------------------------------------
  # 8. Detailed LTR TFBM motif enrichment
  # ------------------------------------------------------------
  
  plot_enrichment_table(
    x$ltr5_tfbm_motif_enrichment,
    filename = "ltr5_tfbm_motif_enrichment.png",
    title = "5' LTR TFBM motif enrichment"
  )
  
  plot_enrichment_table(
    x$ltr3_tfbm_motif_enrichment,
    filename = "ltr3_tfbm_motif_enrichment.png",
    title = "3' LTR TFBM motif enrichment"
  )
  
  plot_enrichment_table(
    x$any_ltr_tfbm_motif_enrichment,
    filename = "any_ltr_tfbm_motif_enrichment.png",
    title = "Any-LTR TFBM motif enrichment"
  )
  
  # ------------------------------------------------------------
  # 9. Detailed RBP motif enrichment
  # ------------------------------------------------------------
  
  plot_enrichment_table(
    x$rbp_motif_enrichment,
    filename = "rbp_motif_enrichment.png",
    title = "RBP motif enrichment"
  )
  
  message("HERVarium comparison plots written to: ", plot_dir)
  
  invisible(plot_dir)
}