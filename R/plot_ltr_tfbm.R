plot_ltr_tfbm <- function(x,
                          output_dir = NULL,
                          top_n_motifs = 25,
                          top_n_subfamilies = 20) {
  
  if (!requireNamespace("ggplot2", quietly = TRUE)) {
    stop("Package 'ggplot2' is required.")
  }
  
  if (inherits(x, "HERVarium_annotation")) {
    features <- x$features
    tfbm_hits <- x$ltr_tfbm_hits
    herv_tfbm_summary <- x$herv_tfbm_summary
    top_tfbm_motifs <- x$top_tfbm_motifs
  } else {
    stop("x must be a HERVarium_annotation object after add_ltr_tfbm_details().")
  }
  
  if (is.null(output_dir)) {
    output_dir <- x$output_dir
  }
  
  if (is.null(output_dir)) {
    stop("Please provide output_dir or use an object with output_dir.")
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  if (is.null(tfbm_hits) || nrow(tfbm_hits) == 0) {
    stop("No detailed TFBM hits found. Run add_ltr_tfbm_details() first.")
  }
  
  if (is.null(herv_tfbm_summary) || nrow(herv_tfbm_summary) == 0) {
    stop("No HERV-level TFBM summary found. Run add_ltr_tfbm_details() first.")
  }
  
  # The FIMO file already has a column called "subfamily".
  # We remove it here because we want the subfamily derived from HERV_id/features,
  # otherwise merge() creates subfamily.x/subfamily.y and downstream code breaks.
  if ("subfamily" %in% colnames(tfbm_hits)) {
    tfbm_hits$subfamily <- NULL
  }
  
  # ------------------------------------------------------------
  # 1. LTR5 vs LTR3 TFBM burden per HERV
  # ------------------------------------------------------------
  
  p <- ggplot2::ggplot(
    herv_tfbm_summary,
    ggplot2::aes(
      x = ltr5_n_tfbm_hits,
      y = ltr3_n_tfbm_hits
    )
  ) +
    ggplot2::geom_point(alpha = 0.75) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
    ggplot2::labs(
      title = "LTR5 versus LTR3 TFBM burden",
      x = "5' LTR TFBM hits",
      y = "3' LTR TFBM hits"
    ) +
    ggplot2::theme_bw()
  
  ggplot2::ggsave(
    file.path(output_dir, "ltr5_vs_ltr3_tfbm_burden.png"),
    p,
    width = 6,
    height = 5,
    dpi = 300
  )
  
  # ------------------------------------------------------------
  # 2. LTR asymmetry distribution
  # ------------------------------------------------------------
  
  herv_tfbm_summary$ltr_tfbm_asymmetry <-
    herv_tfbm_summary$ltr5_n_tfbm_hits -
    herv_tfbm_summary$ltr3_n_tfbm_hits
  
  p <- ggplot2::ggplot(
    herv_tfbm_summary,
    ggplot2::aes(x = ltr_tfbm_asymmetry)
  ) +
    ggplot2::geom_histogram(bins = 30) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
    ggplot2::labs(
      title = "LTR TFBM burden asymmetry",
      subtitle = "Positive values indicate higher 5' LTR burden; negative values indicate higher 3' LTR burden",
      x = "5' LTR hits - 3' LTR hits",
      y = "Number of HERVs"
    ) +
    ggplot2::theme_bw()
  
  ggplot2::ggsave(
    file.path(output_dir, "ltr_tfbm_asymmetry_distribution.png"),
    p,
    width = 7,
    height = 4,
    dpi = 300
  )
  
  # ------------------------------------------------------------
  # 3. Top motifs by number of HERVs
  # ------------------------------------------------------------
  
  motif_herv_counts <- unique(tfbm_hits[, c("motif_alt_id", "HERV_id")])
  motif_herv_counts <- as.data.frame(table(motif_herv_counts$motif_alt_id))
  colnames(motif_herv_counts) <- c("motif_alt_id", "n_hervs")
  motif_herv_counts <- motif_herv_counts[
    order(motif_herv_counts$n_hervs, decreasing = TRUE),
  ]
  
  top_motifs <- head(motif_herv_counts, top_n_motifs)
  
  p <- ggplot2::ggplot(
    top_motifs,
    ggplot2::aes(x = reorder(motif_alt_id, n_hervs), y = n_hervs)
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste0("Top ", top_n_motifs, " TF motifs by number of HERVs"),
      x = NULL,
      y = "Number of HERVs carrying motif"
    ) +
    ggplot2::theme_bw()
  
  ggplot2::ggsave(
    file.path(output_dir, "top_tfbm_motifs_by_hervs.png"),
    p,
    width = 7,
    height = 7,
    dpi = 300
  )
  
  # ------------------------------------------------------------
  # 4. Top motifs split by LTR5/LTR3
  # ------------------------------------------------------------
  
  motif_ltr_counts <- unique(tfbm_hits[, c("motif_alt_id", "HERV_id", "ltr_position")])
  motif_ltr_counts <- as.data.frame(table(
    motif_ltr_counts$motif_alt_id,
    motif_ltr_counts$ltr_position
  ))
  
  colnames(motif_ltr_counts) <- c("motif_alt_id", "ltr_position", "n_hervs")
  
  selected_motifs <- top_motifs$motif_alt_id
  motif_ltr_counts <- motif_ltr_counts[
    motif_ltr_counts$motif_alt_id %in% selected_motifs,
  ]
  
  p <- ggplot2::ggplot(
    motif_ltr_counts,
    ggplot2::aes(
      x = reorder(motif_alt_id, n_hervs),
      y = n_hervs,
      fill = ltr_position
    )
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = "Top TF motifs split by LTR position",
      x = NULL,
      y = "Number of HERVs carrying motif",
      fill = "LTR"
    ) +
    ggplot2::theme_bw()
  
  ggplot2::ggsave(
    file.path(output_dir, "top_tfbm_motifs_ltr5_ltr3.png"),
    p,
    width = 8,
    height = 7,
    dpi = 300
  )
  
  # ------------------------------------------------------------
  # 5. Subfamily-level LTR5/LTR3 mean burden
  # ------------------------------------------------------------
  
  subfamily_info <- unique(features[, c("HERV_id", "subfamily")])
  
  hsum <- merge(
    herv_tfbm_summary,
    subfamily_info,
    by = "HERV_id",
    all.x = TRUE
  )
  
  subfamily_counts <- as.data.frame(table(hsum$subfamily))
  colnames(subfamily_counts) <- c("subfamily", "n_hervs")
  subfamily_counts <- subfamily_counts[
    order(subfamily_counts$n_hervs, decreasing = TRUE),
  ]
  
  selected_subfamilies <- head(subfamily_counts$subfamily, top_n_subfamilies)
  
  hsum_top <- hsum[hsum$subfamily %in% selected_subfamilies, ]
  
  sub_ltr <- do.call(
    rbind,
    lapply(split(hsum_top, hsum_top$subfamily), function(z) {
      data.frame(
        subfamily = unique(z$subfamily),
        n_hervs = nrow(z),
        mean_ltr5_hits = mean(z$ltr5_n_tfbm_hits, na.rm = TRUE),
        mean_ltr3_hits = mean(z$ltr3_n_tfbm_hits, na.rm = TRUE),
        mean_total_hits = mean(z$total_n_tfbm_hits, na.rm = TRUE),
        stringsAsFactors = FALSE
      )
    })
  )
  
  sub_ltr_long <- data.frame(
    subfamily = rep(sub_ltr$subfamily, times = 2),
    ltr_position = rep(c("LTR5", "LTR3"), each = nrow(sub_ltr)),
    mean_hits = c(sub_ltr$mean_ltr5_hits, sub_ltr$mean_ltr3_hits),
    n_hervs = rep(sub_ltr$n_hervs, times = 2),
    stringsAsFactors = FALSE
  )
  
  p <- ggplot2::ggplot(
    sub_ltr_long,
    ggplot2::aes(
      x = reorder(subfamily, mean_hits),
      y = mean_hits,
      fill = ltr_position
    )
  ) +
    ggplot2::geom_col(position = "dodge") +
    ggplot2::coord_flip() +
    ggplot2::labs(
      title = paste0("Mean LTR TFBM burden by subfamily"),
      subtitle = paste0("Top ", top_n_subfamilies, " subfamilies by representation"),
      x = NULL,
      y = "Mean number of TFBM hits",
      fill = "LTR"
    ) +
    ggplot2::theme_bw()
  
  ggplot2::ggsave(
    file.path(output_dir, "subfamily_ltr5_ltr3_mean_tfbm_burden.png"),
    p,
    width = 8,
    height = 7,
    dpi = 300
  )
  
  # ------------------------------------------------------------
  # 6. Dot plot: TF motifs × subfamilies
  # ------------------------------------------------------------
  
  top_tf_names <- top_motifs$motif_alt_id
  top_subfamilies <- selected_subfamilies
  
  tf_sub <- merge(
    tfbm_hits,
    subfamily_info,
    by = "HERV_id",
    all.x = TRUE
  )
  
  tf_sub <- tf_sub[
    tf_sub$motif_alt_id %in% top_tf_names &
      tf_sub$subfamily %in% top_subfamilies,
  ]
  
  tf_sub_unique <- unique(tf_sub[, c("motif_alt_id", "subfamily", "HERV_id")])
  
  tf_sub_counts <- as.data.frame(table(
    tf_sub_unique$motif_alt_id,
    tf_sub_unique$subfamily
  ))
  
  colnames(tf_sub_counts) <- c("motif_alt_id", "subfamily", "n_hervs")
  
  subfamily_denominator <- as.data.frame(table(subfamily_info$subfamily))
  colnames(subfamily_denominator) <- c("subfamily", "n_total_subfamily")
  
  tf_sub_counts <- merge(
    tf_sub_counts,
    subfamily_denominator,
    by = "subfamily",
    all.x = TRUE
  )
  
  tf_sub_counts$percent_hervs <-
    100 * tf_sub_counts$n_hervs / tf_sub_counts$n_total_subfamily
  
  tf_sub_counts <- tf_sub_counts[tf_sub_counts$n_hervs > 0, ]
  
  p <- ggplot2::ggplot(
    tf_sub_counts,
    ggplot2::aes(
      x = subfamily,
      y = motif_alt_id,
      size = percent_hervs,
      alpha = n_hervs
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::labs(
      title = "TF motif landscape by HERV subfamily",
      subtitle = "Dot size indicates percent of HERVs in the subfamily carrying the motif",
      x = "HERV subfamily",
      y = "TF motif",
      size = "% HERVs",
      alpha = "n HERVs"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  
  ggplot2::ggsave(
    file.path(output_dir, "tfbm_motif_by_subfamily_dotplot.png"),
    p,
    width = 10,
    height = 8,
    dpi = 300
  )
  
  # ------------------------------------------------------------
  # 7. LTR-position-specific dot plot: motif × LTR position
  # ------------------------------------------------------------
  
  tf_ltr_unique <- unique(tfbm_hits[, c("motif_alt_id", "ltr_position", "HERV_id")])
  tf_ltr_unique <- tf_ltr_unique[tf_ltr_unique$motif_alt_id %in% top_tf_names, ]
  
  tf_ltr_counts <- as.data.frame(table(
    tf_ltr_unique$motif_alt_id,
    tf_ltr_unique$ltr_position
  ))
  
  colnames(tf_ltr_counts) <- c("motif_alt_id", "ltr_position", "n_hervs")
  tf_ltr_counts <- tf_ltr_counts[tf_ltr_counts$n_hervs > 0, ]
  
  p <- ggplot2::ggplot(
    tf_ltr_counts,
    ggplot2::aes(
      x = ltr_position,
      y = motif_alt_id,
      size = n_hervs
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::labs(
      title = "LTR-position-specific TF motif landscape",
      x = "LTR position",
      y = "TF motif",
      size = "Number of HERVs"
    ) +
    ggplot2::theme_bw()
  
  ggplot2::ggsave(
    file.path(output_dir, "tfbm_motif_by_ltr_position_dotplot.png"),
    p,
    width = 6,
    height = 8,
    dpi = 300
  )
  
  
  # ------------------------------------------------------------
  # 8. Top TF motifs by subfamily and LTR position
  # ------------------------------------------------------------
  
  tf_sub_ltr <- merge(
    tfbm_hits,
    subfamily_info,
    by = "HERV_id",
    all.x = TRUE
  )
  
  tf_sub_ltr <- tf_sub_ltr[
    !is.na(tf_sub_ltr$subfamily) &
      !is.na(tf_sub_ltr$motif_alt_id) &
      !is.na(tf_sub_ltr$ltr_position),
  ]
  
  # Count motif presence by number of unique HERVs,
  # not raw number of motif hits
  tf_sub_ltr_unique <- unique(
    tf_sub_ltr[, c("subfamily", "ltr_position", "motif_alt_id", "HERV_id")]
  )
  
  tf_sub_ltr_counts <- as.data.frame(table(
    tf_sub_ltr_unique$subfamily,
    tf_sub_ltr_unique$ltr_position,
    tf_sub_ltr_unique$motif_alt_id
  ))
  
  colnames(tf_sub_ltr_counts) <- c(
    "subfamily",
    "ltr_position",
    "motif_alt_id",
    "n_hervs"
  )
  
  tf_sub_ltr_counts <- tf_sub_ltr_counts[tf_sub_ltr_counts$n_hervs > 0, ]
  
  # Denominator: number of HERVs per subfamily in the input list
  subfamily_denominator <- as.data.frame(table(subfamily_info$subfamily))
  colnames(subfamily_denominator) <- c("subfamily", "n_total_subfamily")
  
  tf_sub_ltr_counts <- merge(
    tf_sub_ltr_counts,
    subfamily_denominator,
    by = "subfamily",
    all.x = TRUE
  )
  
  tf_sub_ltr_counts$percent_hervs <-
    100 * tf_sub_ltr_counts$n_hervs / tf_sub_ltr_counts$n_total_subfamily
  
  # Restrict to top represented subfamilies
  top_subfamilies_for_tfbm <- head(
    subfamily_counts$subfamily,
    top_n_subfamilies
  )
  
  tf_sub_ltr_counts_top <- tf_sub_ltr_counts[
    tf_sub_ltr_counts$subfamily %in% top_subfamilies_for_tfbm,
  ]
  
  # Select top motifs separately for LTR5 and LTR3
  top_motifs_by_ltr <- do.call(
    rbind,
    lapply(split(tf_sub_ltr_counts_top, tf_sub_ltr_counts_top$ltr_position), function(z) {
      motif_totals <- aggregate(
        n_hervs ~ motif_alt_id,
        data = z,
        FUN = sum
      )
      
      motif_totals <- motif_totals[
        order(motif_totals$n_hervs, decreasing = TRUE),
      ]
      
      motif_totals <- head(motif_totals, top_n_motifs)
      motif_totals$ltr_position <- unique(z$ltr_position)
      
      motif_totals
    })
  )
  
  selected_motif_ltr_pairs <- unique(
    top_motifs_by_ltr[, c("motif_alt_id", "ltr_position")]
  )
  
  tf_sub_ltr_counts_top <- merge(
    tf_sub_ltr_counts_top,
    selected_motif_ltr_pairs,
    by = c("motif_alt_id", "ltr_position")
  )
  
  # Dot plot: subfamily x motif, faceted by 5'/3' LTR
  p <- ggplot2::ggplot(
    tf_sub_ltr_counts_top,
    ggplot2::aes(
      x = subfamily,
      y = motif_alt_id,
      size = percent_hervs,
      alpha = n_hervs
    )
  ) +
    ggplot2::geom_point() +
    ggplot2::facet_wrap(~ ltr_position, scales = "free_y") +
    ggplot2::labs(
      title = "Top TF motifs by HERV subfamily and LTR position",
      subtitle = "Dot size indicates percent of HERVs in each subfamily carrying the motif",
      x = "HERV subfamily",
      y = "TF motif",
      size = "% HERVs",
      alpha = "n HERVs"
    ) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
    )
  
  ggplot2::ggsave(
    file.path(output_dir, "tfbm_top_motifs_by_subfamily_and_ltr_dotplot.png"),
    p,
    width = 12,
    height = 9,
    dpi = 300
  )
  
  # Bar plot: top motif-subfamily-LTR combinations
  tf_sub_ltr_bar <- tf_sub_ltr_counts_top[
    order(tf_sub_ltr_counts_top$n_hervs, decreasing = TRUE),
  ]
  
  tf_sub_ltr_bar$label <- paste0(
    tf_sub_ltr_bar$subfamily,
    " | ",
    tf_sub_ltr_bar$motif_alt_id
  )
  
  tf_sub_ltr_bar <- do.call(
    rbind,
    lapply(split(tf_sub_ltr_bar, tf_sub_ltr_bar$ltr_position), function(z) {
      head(z, top_n_motifs)
    })
  )
  
  p <- ggplot2::ggplot(
    tf_sub_ltr_bar,
    ggplot2::aes(
      x = reorder(label, n_hervs),
      y = n_hervs
    )
  ) +
    ggplot2::geom_col() +
    ggplot2::coord_flip() +
    ggplot2::facet_wrap(~ ltr_position, scales = "free_y") +
    ggplot2::labs(
      title = "Top subfamily-specific TF motifs in 5' and 3' LTRs",
      x = NULL,
      y = "Number of HERVs carrying motif"
    ) +
    ggplot2::theme_bw()
  
  ggplot2::ggsave(
    file.path(output_dir, "tfbm_top_subfamily_motifs_ltr5_ltr3_barplot.png"),
    p,
    width = 12,
    height = 9,
    dpi = 300
  )
  
  write.table(
    tf_sub_ltr_counts,
    file = file.path(output_dir, "plot_table_tfbm_by_subfamily_and_ltr.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  # ------------------------------------------------------------
  # 8. Write plot-supporting tables
  # ------------------------------------------------------------
  
  write.table(
    motif_herv_counts,
    file = file.path(output_dir, "plot_table_top_tfbm_motifs_by_hervs.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  write.table(
    motif_ltr_counts,
    file = file.path(output_dir, "plot_table_top_tfbm_motifs_by_ltr_position.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  write.table(
    sub_ltr,
    file = file.path(output_dir, "plot_table_subfamily_ltr_tfbm_burden.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  write.table(
    tf_sub_counts,
    file = file.path(output_dir, "plot_table_tfbm_motif_by_subfamily.tsv"),
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )
  
  message("LTR TFBM plots written to: ", output_dir)
  
  invisible(output_dir)
}