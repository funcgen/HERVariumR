# ------------------------------------------------------------
# HERVariumR concise interactive dashboards
# ------------------------------------------------------------

.dashboard_as_bool <- function(x, n = NULL) {
  if (is.null(x)) {
    if (is.null(n)) return(logical(0))
    return(rep(FALSE, n))
  }

  if (is.logical(x)) {
    x[is.na(x)] <- FALSE
    return(x)
  }

  z <- tolower(as.character(x))
  out <- z %in% c("true", "t", "1", "yes")
  out[is.na(out)] <- FALSE
  out
}

.dashboard_numeric <- function(x, n = NULL) {
  if (is.null(x)) {
    if (is.null(n)) return(numeric(0))
    return(rep(0, n))
  }
  out <- suppressWarnings(as.numeric(x))
  out[is.na(out)] <- 0
  out
}

.dashboard_safe_dt <- function(x,
                               page_length = 10,
                               filter = "top") {
  if (is.null(x) || !is.data.frame(x) || nrow(x) == 0) {
    return(htmltools::HTML("<p><em>No results available for this section.</em></p>"))
  }

  x_display <- as.data.frame(x, stringsAsFactors = FALSE)
  numeric_cols <- vapply(x_display, is.numeric, logical(1))
  if (any(numeric_cols)) {
    x_display[numeric_cols] <- lapply(
      x_display[numeric_cols],
      function(z) round(z, 2)
    )
  }

  DT::datatable(
    x_display,
    filter = filter,
    extensions = "Buttons",
    rownames = FALSE,
    options = list(
      pageLength = page_length,
      scrollX = TRUE,
      autoWidth = TRUE,
      dom = "Bfrtip",
      buttons = c("copy", "csv", "excel")
    )
  )
}

.dashboard_cards <- function(cards) {
  if (is.null(cards) || nrow(cards) == 0) return(NULL)

  card_nodes <- lapply(seq_len(nrow(cards)), function(i) {
    htmltools::div(
      class = "herv-card",
      htmltools::div(class = "herv-card-value", as.character(cards$value[i])),
      htmltools::div(class = "herv-card-label", as.character(cards$label[i]))
    )
  })

  htmltools::div(class = "herv-card-grid", card_nodes)
}

.dashboard_wrap_label <- function(x, width = 42) {
  vapply(
    as.character(x),
    function(z) paste(strwrap(z, width = width), collapse = "<br>"),
    character(1)
  )
}

.dashboard_bar <- function(df,
                           label_col,
                           value_col,
                           axis_title,
                           group_col = NULL,
                           hover_col = NULL,
                           colors = NULL,
                           max_height = 560) {
  if (is.null(df) || nrow(df) == 0) {
    return(htmltools::HTML("<p><em>No data available for this plot.</em></p>"))
  }

  keep <- !is.na(df[[value_col]])
  df <- df[keep, , drop = FALSE]
  if (nrow(df) == 0) {
    return(htmltools::HTML("<p><em>No finite values available for this plot.</em></p>"))
  }

  plot_df <- data.frame(
    label = as.character(df[[label_col]]),
    value = suppressWarnings(as.numeric(df[[value_col]])),
    stringsAsFactors = FALSE
  )

  plot_df$hover <- if (is.null(hover_col)) {
    paste0(plot_df$label, "<br>", axis_title, ": ", signif(plot_df$value, 4))
  } else {
    as.character(df[[hover_col]])
  }

  grouped <- !is.null(group_col) && group_col %in% colnames(df)
  if (grouped) plot_df$group <- as.character(df[[group_col]])

  label_order <- unique(rev(plot_df$label))
  n_labels <- length(unique(plot_df$label))
  plot_height <- min(max(250, 28 * n_labels + 95), max_height)

  if (grouped) {
    present_groups <- unique(plot_df$group)
    selected_colors <- if (!is.null(colors)) colors[present_groups] else NULL

    p <- plotly::plot_ly(
      data = plot_df,
      height = plot_height,
      x = ~value,
      y = ~label,
      type = "bar",
      orientation = "h",
      color = ~group,
      colors = selected_colors,
      hovertext = ~hover,
      hoverinfo = "text",
      textposition = "none",
      marker = list(line = list(width = 0))
    )
  } else {
    p <- plotly::plot_ly(
      data = plot_df,
      height = plot_height,
      x = ~value,
      y = ~label,
      type = "bar",
      orientation = "h",
      hovertext = ~hover,
      hoverinfo = "text",
      textposition = "none",
      marker = list(line = list(width = 0)),
      showlegend = FALSE
    )
  }

  p <- plotly::layout(
    p,
    autosize = TRUE,
    height = plot_height,
    barmode = if (grouped) "group" else "overlay",
    xaxis = list(title = axis_title, automargin = TRUE, zeroline = TRUE),
    yaxis = list(
      title = "",
      categoryorder = "array",
      categoryarray = label_order,
      automargin = TRUE
    ),
    legend = list(
      orientation = "h",
      x = 0,
      y = -0.10,
      xanchor = "left",
      yanchor = "top"
    ),
    margin = list(l = 245, r = 35, t = 10, b = if (grouped) 65 else 45)
  )

  plotly::config(p, responsive = TRUE, displaylogo = FALSE)
}

.dashboard_binary_plot <- function(df,
                                   top_n = 15,
                                   padj_cutoff = 0.05) {
  if (is.null(df) || nrow(df) == 0) {
    return(htmltools::HTML("<p><em>No binary enrichment results available.</em></p>"))
  }

  df <- df[!is.na(df$p_value), , drop = FALSE]
  if (nrow(df) == 0) {
    return(htmltools::HTML("<p><em>No valid binary enrichment p-values.</em></p>"))
  }

  a <- df$foreground_yes + 0.5
  b <- df$foreground_no + 0.5
  c <- df$background_yes + 0.5
  d <- df$background_no + 0.5
  df$effect <- log2((a / b) / (c / d))

  df$category <- "Not significant"
  sig <- !is.na(df$padj) & df$padj < padj_cutoff
  df$category[sig & df$effect > 0] <- "Significant enriched"
  df$category[sig & df$effect < 0] <- "Significant depleted"

  df <- df[order(df$p_value), , drop = FALSE]
  df <- head(df, top_n)
  df$label <- .dashboard_wrap_label(.clean_result_feature_label(df$feature))
  df$hover <- paste0(
    "Feature: ", .clean_result_feature_label(df$feature),
    "<br>Foreground: ", df$foreground_yes, "/", df$foreground_total,
    " (", round(df$foreground_percent, 1), "%)",
    "<br>Background: ", df$background_yes, "/", df$background_total,
    " (", round(df$background_percent, 1), "%)",
    "<br>log2 OR: ", signif(df$effect, 4),
    "<br>p-value: ", signif(df$p_value, 4),
    "<br>FDR: ", signif(df$padj, 4)
  )

  palette <- c(
    "Significant enriched" = "#D55E00",
    "Significant depleted" = "#0072B2",
    "Not significant" = "#BDBDBD"
  )

  .dashboard_bar(
    df,
    label_col = "label",
    value_col = "effect",
    axis_title = "log2 odds ratio",
    group_col = "category",
    hover_col = "hover",
    colors = palette
  )
}

.dashboard_numeric_plot <- function(df,
                                    top_n = 15,
                                    padj_cutoff = 0.05) {
  if (is.null(df) || nrow(df) == 0) {
    return(htmltools::HTML("<p><em>No numeric comparison results available.</em></p>"))
  }

  df <- df[!is.na(df$p_value), , drop = FALSE]
  if (nrow(df) == 0) {
    return(htmltools::HTML("<p><em>No valid numeric comparison p-values.</em></p>"))
  }

  if (all(is.na(df$delta_median) | abs(df$delta_median) < 1e-9)) {
    return(htmltools::HTML("<p><em>No non-zero median shifts were observed.</em></p>"))
  }

  df$category <- "Not significant"
  sig <- !is.na(df$padj) & df$padj < padj_cutoff
  df$category[sig & df$delta_median > 0] <- "Significant enriched"
  df$category[sig & df$delta_median < 0] <- "Significant depleted"

  df <- df[order(df$p_value), , drop = FALSE]
  df <- head(df, top_n)
  df$label <- .dashboard_wrap_label(.clean_result_feature_label(df$feature))
  df$hover <- paste0(
    "Feature: ", .clean_result_feature_label(df$feature),
    "<br>Foreground median: ", signif(df$foreground_median, 4),
    "<br>Background median: ", signif(df$background_median, 4),
    "<br>Delta median: ", signif(df$delta_median, 4),
    "<br>p-value: ", signif(df$p_value, 4),
    "<br>FDR: ", signif(df$padj, 4)
  )

  palette <- c(
    "Significant enriched" = "#D55E00",
    "Significant depleted" = "#0072B2",
    "Not significant" = "#BDBDBD"
  )

  .dashboard_bar(
    df,
    label_col = "label",
    value_col = "delta_median",
    axis_title = "Delta median",
    group_col = "category",
    hover_col = "hover",
    colors = palette
  )
}

.dashboard_key_plot <- function(df,
                                top_n = 15,
                                axis_title = NULL) {
  if (is.null(df) || nrow(df) == 0) {
    return(htmltools::HTML("<p><em>No key results available for this plot.</em></p>"))
  }

  df <- df[!is.na(df$effect_size), , drop = FALSE]
  if (nrow(df) == 0) {
    return(htmltools::HTML("<p><em>No finite effects available for this plot.</em></p>"))
  }
  if (all(abs(df$effect_size) < 1e-9)) {
    return(htmltools::HTML("<p><em>No non-zero effects were observed for this section.</em></p>"))
  }

  df <- df[order(df$padj, df$p_value, na.last = TRUE), , drop = FALSE]
  df <- head(df, top_n)
  df$label <- .dashboard_wrap_label(df$feature_label)
  df$category <- "Not significant"
  df$category[df$significant & df$effect_size > 0] <- "Significant enriched"
  df$category[df$significant & df$effect_size < 0] <- "Significant depleted"
  df$hover <- paste0(
    "Feature: ", df$feature_label,
    "<br>Layer: ", df$layer,
    "<br>Effect: ", signif(df$effect_size, 4),
    "<br>Foreground: ", signif(df$foreground_value, 4),
    "<br>Background: ", signif(df$background_value, 4),
    "<br>p-value: ", signif(df$p_value, 4),
    "<br>FDR: ", signif(df$padj, 4)
  )

  palette <- c(
    "Significant enriched" = "#D55E00",
    "Significant depleted" = "#0072B2",
    "Not significant" = "#BDBDBD"
  )

  if (is.null(axis_title)) {
    axis_title <- if (length(unique(df$effect_type)) == 1) {
      unique(df$effect_type)
    } else {
      "Effect size"
    }
  }

  .dashboard_bar(
    df,
    label_col = "label",
    value_col = "effect_size",
    axis_title = axis_title,
    group_col = "category",
    hover_col = "hover",
    colors = palette
  )
}

.dashboard_profile_cards <- function(anno) {
  f <- anno$features
  n <- nrow(f)

  count_col <- function(col) {
    if (!col %in% colnames(f)) return(NA_integer_)
    sum(.dashboard_as_bool(f[[col]]), na.rm = TRUE)
  }

  cards <- data.frame(
    label = c(
      "Input IDs", "Matched HERVs", "Missing IDs", "With domains",
      "With 5\u2032 LTR", "With 3\u2032 LTR", "With both LTRs",
      "Internal external evidence", "5\u2032/3\u2032 LTR external evidence",
      "5\u2032/3\u2032 LTR TF ChIP-seq"
    ),
    value = c(
      if (!is.null(anno$stats$n_input)) anno$stats$n_input[1] else n,
      n,
      length(anno$missing_ids),
      count_col("has_domain"),
      count_col("has_ltr5"),
      count_col("has_ltr3"),
      count_col("has_both_ltrs"),
      count_col("has_internal_encode_evidence"),
      sum(
        .dashboard_as_bool(f$ltr5_encode_evidence_including_ccre, n) |
          .dashboard_as_bool(f$ltr3_encode_evidence_including_ccre, n)
      ),
      sum(
        .dashboard_as_bool(f$has_ltr5_tfchip_overlap, n) |
          .dashboard_as_bool(f$has_ltr3_tfchip_overlap, n)
      )
    ),
    stringsAsFactors = FALSE
  )

  .dashboard_cards(cards)
}

.dashboard_profile_candidate_table <- function(anno, top_n = 20) {
  f <- anno$features
  n <- nrow(f)
  if (n == 0) return(.dashboard_safe_dt(data.frame()))

  ncol0 <- function(col) {
    if (!col %in% colnames(f)) return(rep(0, n))
    .dashboard_numeric(f[[col]], n)
  }
  bcol0 <- function(col) {
    if (!col %in% colnames(f)) return(rep(FALSE, n))
    .dashboard_as_bool(f[[col]], n)
  }

  regulatory_score <-
    ncol0("internal_encode_detected_layer_count") +
    ncol0("ltr5_encode_detected_layer_count") +
    ncol0("ltr3_encode_detected_layer_count") +
    bcol0("ltr5_encode_ccre_overlap") +
    bcol0("ltr3_encode_ccre_overlap") +
    bcol0("has_ltr5_tfchip_overlap") +
    bcol0("has_ltr3_tfchip_overlap")

  structural_score <- pmin(ncol0("domain_count"), 3) +
    bcol0("has_both_ltrs")

  f$candidate_score <- regulatory_score + 0.5 * structural_score
  f <- f[order(f$candidate_score, ncol0("domain_count"), decreasing = TRUE), , drop = FALSE]
  f <- head(f, top_n)

  selected <- c(
    "HERV_id", "subfamily", "candidate_score", "domain_count", "domains_type",
    "has_ltr5", "has_ltr3", "ltr5_tfbm_burden", "ltr3_tfbm_burden",
    "rbp_burden", "internal_encode_detected_layers",
    "ltr5_encode_detected_layers", "ltr3_encode_detected_layers",
    "ltr5_tfchip_tf_list", "ltr3_tfchip_tf_list"
  )
  selected <- selected[selected %in% colnames(f)]
  .dashboard_safe_dt(f[, selected, drop = FALSE], page_length = min(10, top_n))
}

.dashboard_profile_subfamilies <- function(anno, top_n = 15) {
  f <- anno$features
  z <- as.data.frame(sort(table(f$subfamily), decreasing = TRUE), stringsAsFactors = FALSE)
  colnames(z) <- c("label", "value")
  z <- z[z$label != "" & z$label != ".", , drop = FALSE]
  z <- head(z, top_n)
  z$hover <- paste0(z$label, "<br>HERVs: ", z$value)
  .dashboard_bar(z, "label", "value", "Number of HERVs", hover_col = "hover")
}

.dashboard_parse_domain_types <- function(x) {
  out <- character(0)
  for (value in as.character(x)) {
    if (is.na(value) || value == "" || value == ".") next
    entries <- unlist(strsplit(value, ";", fixed = TRUE))
    entries <- entries[entries != "" & entries != "."]
    if (length(entries) == 0) next

    types <- vapply(entries, function(entry) {
      z <- sub("^[^|]*\\|", "", entry)
      z <- sub(":.*$", "", z)
      z <- sub("\\|.*$", "", z)
      z
    }, character(1))

    out <- c(out, unique(types))
  }

  out[out != "" & out != "."]
}

.dashboard_profile_domains <- function(anno) {
  f <- anno$features
  types <- .dashboard_parse_domain_types(f$domains_type)

  if (length(types) == 0) {
    return(htmltools::HTML("<p><em>No retained internal domains in this list.</em></p>"))
  }

  clean <- c(
    GAG = "Gag", ENV = "Env", AP = "Protease", DUT = "dUTPase",
    INT = "Integrase", RT = "RT", RH = "RNaseH"
  )
  replace <- types %in% names(clean)
  types[replace] <- unname(clean[types[replace]])

  z <- as.data.frame(sort(table(types), decreasing = TRUE), stringsAsFactors = FALSE)
  colnames(z) <- c("label", "value")
  z$hover <- paste0(z$label, "<br>HERVs carrying domain: ", z$value)
  .dashboard_bar(z, "label", "value", "Number of HERVs", hover_col = "hover")
}

.dashboard_profile_ltr_structure <- function(anno) {
  f <- anno$features
  n <- nrow(f)
  z <- data.frame(
    label = c("5\u2032 LTR present", "3\u2032 LTR present", "Both LTRs present"),
    value = c(
      sum(.dashboard_as_bool(f$has_ltr5, n)),
      sum(.dashboard_as_bool(f$has_ltr3, n)),
      sum(.dashboard_as_bool(f$has_both_ltrs, n))
    ),
    stringsAsFactors = FALSE
  )
  z$hover <- paste0(z$label, "<br>HERVs: ", z$value)
  .dashboard_bar(z, "label", "value", "Number of HERVs", hover_col = "hover")
}

.dashboard_profile_context <- function(anno) {
  f <- anno$features
  n <- nrow(f)
  mapping <- c(
    Intergenic = "is_intergenic",
    Exon = "overlaps_exon",
    Intron = "overlaps_intron",
    CDS = "overlaps_cds",
    UTR = "overlaps_utr",
    lncRNA = "overlaps_lncRNA",
    `Protein-coding` = "overlaps_protein_coding"
  )

  z <- data.frame(
    label = names(mapping),
    value = vapply(mapping, function(col) {
      if (!col %in% colnames(f)) return(0)
      sum(.dashboard_as_bool(f[[col]], n))
    }, numeric(1)),
    stringsAsFactors = FALSE
  )
  z$hover <- paste0(z$label, "<br>HERVs: ", z$value)
  .dashboard_bar(z, "label", "value", "Number of HERVs", hover_col = "hover")
}

.dashboard_detected_layer <- function(x, layer) {
  if (is.null(x)) return(logical(0))
  grepl(paste0("(^|;)", layer, "($|;)"), as.character(x))
}

.dashboard_profile_regulatory <- function(anno) {
  f <- anno$features
  n <- nrow(f)
  if (n == 0) return(htmltools::HTML("<p><em>No matched HERVs.</em></p>"))

  evidence_cols <- c(
    "internal_encode_detected_layers",
    "ltr5_encode_detected_layers",
    "ltr3_encode_detected_layers",
    "ltr5_encode_ccre_overlap",
    "ltr3_encode_ccre_overlap",
    "has_ltr5_tfchip_overlap",
    "has_ltr3_tfchip_overlap"
  )
  if (!any(evidence_cols %in% colnames(f))) {
    return(htmltools::HTML("<p><em>External regulatory-evidence layers were not enabled for this profile.</em></p>"))
  }

  rows <- list()
  add_row <- function(region, evidence, values) {
    values <- .dashboard_as_bool(values, n)
    rows[[length(rows) + 1]] <<- data.frame(
      region = region,
      evidence = evidence,
      count = sum(values),
      percent = 100 * sum(values) / n,
      stringsAsFactors = FALSE
    )
  }

  region_specs <- list(
    Internal = list(prefix = "internal_encode", layers_col = "internal_encode_detected_layers"),
    "5\u2032 LTR" = list(prefix = "ltr5_encode", layers_col = "ltr5_encode_detected_layers"),
    "3\u2032 LTR" = list(prefix = "ltr3_encode", layers_col = "ltr3_encode_detected_layers")
  )

  for (region in names(region_specs)) {
    spec <- region_specs[[region]]
    layer_values <- if (spec$layers_col %in% colnames(f)) f[[spec$layers_col]] else rep(".", n)

    add_row(region, "DNase", .dashboard_detected_layer(layer_values, "DNase"))
    add_row(region, "H3K27ac", .dashboard_detected_layer(layer_values, "H3K27ac"))
    add_row(region, "H3K4me3", .dashboard_detected_layer(layer_values, "H3K4me3"))
    add_row(region, "Transcription", .dashboard_detected_layer(layer_values, "transcription"))

    same_col <- paste0(spec$prefix, "_same_strand_transcription_detected")
    opposite_col <- paste0(spec$prefix, "_opposite_strand_transcription_detected")
    bidirectional_col <- paste0(spec$prefix, "_bidirectional_transcription_detected")

    add_row(region, "Same-strand transcription", if (same_col %in% colnames(f)) f[[same_col]] else rep(FALSE, n))
    add_row(region, "Opposite-strand transcription", if (opposite_col %in% colnames(f)) f[[opposite_col]] else rep(FALSE, n))
    add_row(region, "Bidirectional transcription", if (bidirectional_col %in% colnames(f)) f[[bidirectional_col]] else rep(FALSE, n))

    if (region == "5\u2032 LTR") {
      add_row(region, "cCRE overlap", if ("ltr5_encode_ccre_overlap" %in% colnames(f)) f$ltr5_encode_ccre_overlap else rep(FALSE, n))
      add_row(region, "TF ChIP-seq overlap", if ("has_ltr5_tfchip_overlap" %in% colnames(f)) f$has_ltr5_tfchip_overlap else rep(FALSE, n))
    }
    if (region == "3\u2032 LTR") {
      add_row(region, "cCRE overlap", if ("ltr3_encode_ccre_overlap" %in% colnames(f)) f$ltr3_encode_ccre_overlap else rep(FALSE, n))
      add_row(region, "TF ChIP-seq overlap", if ("has_ltr3_tfchip_overlap" %in% colnames(f)) f$has_ltr3_tfchip_overlap else rep(FALSE, n))
    }
  }

  z <- do.call(rbind, rows)
  z$label <- .dashboard_wrap_label(z$evidence, width = 34)
  z$hover <- paste0(
    "Region: ", z$region,
    "<br>Evidence: ", z$evidence,
    "<br>HERVs: ", z$count, "/", n,
    "<br>Percent: ", round(z$percent, 1), "%"
  )

  .dashboard_bar(
    z,
    label_col = "label",
    value_col = "percent",
    axis_title = "Percent of matched HERVs",
    group_col = "region",
    hover_col = "hover",
    colors = c(Internal = "#4C78A8", "5\u2032 LTR" = "#F58518", "3\u2032 LTR" = "#54A24B"),
    max_height = 620
  )
}

.dashboard_profile_regulatory_table <- function(anno, top_n = 20) {
  f <- anno$features
  n <- nrow(f)
  if (n == 0) return(.dashboard_safe_dt(data.frame()))

  ncol0 <- function(col) {
    if (!col %in% colnames(f)) return(rep(0, n))
    .dashboard_numeric(f[[col]], n)
  }
  bcol0 <- function(col) {
    if (!col %in% colnames(f)) return(rep(FALSE, n))
    .dashboard_as_bool(f[[col]], n)
  }

  f$regulatory_evidence_score <-
    ncol0("internal_encode_detected_layer_count") +
    ncol0("ltr5_encode_detected_layer_count") +
    ncol0("ltr3_encode_detected_layer_count") +
    bcol0("ltr5_encode_ccre_overlap") +
    bcol0("ltr3_encode_ccre_overlap") +
    bcol0("has_ltr5_tfchip_overlap") +
    bcol0("has_ltr3_tfchip_overlap")

  f <- f[order(f$regulatory_evidence_score, decreasing = TRUE), , drop = FALSE]
  f <- head(f, top_n)
  selected <- c(
    "HERV_id", "subfamily", "regulatory_evidence_score",
    "internal_encode_detected_layers", "ltr5_encode_detected_layers",
    "ltr3_encode_detected_layers", "ltr5_encode_ccre_overlap",
    "ltr3_encode_ccre_overlap", "ltr5_tfchip_n_tfs", "ltr3_tfchip_n_tfs",
    "ltr5_tfchip_tf_list", "ltr3_tfchip_tf_list"
  )
  selected <- selected[selected %in% colnames(f)]
  .dashboard_safe_dt(f[, selected, drop = FALSE], page_length = min(10, top_n))
}

.dashboard_split_tfchip_names <- function(x) {
  x <- as.character(x)
  if (length(x) == 0 || is.na(x) || x == "" || x == ".") {
    return(character(0))
  }

  z <- unlist(strsplit(x, "[,;]", perl = TRUE))
  z <- trimws(z)
  unique(z[!is.na(z) & z != "" & z != "."])
}

.dashboard_profile_tfchip_data <- function(anno) {
  f <- anno$features
  if (is.null(f) || nrow(f) == 0) return(data.frame())

  rows <- list()
  add_ltr <- function(list_col, ltr_position) {
    if (!list_col %in% colnames(f)) return(NULL)

    for (i in seq_len(nrow(f))) {
      tfs <- .dashboard_split_tfchip_names(f[[list_col]][i])
      if (length(tfs) == 0) next

      rows[[length(rows) + 1]] <<- data.frame(
        HERV_id = as.character(f$HERV_id[i]),
        ltr_position = ltr_position,
        tf = tfs,
        stringsAsFactors = FALSE
      )
    }
  }

  add_ltr("ltr5_tfchip_tf_list", "5\u2032 LTR")
  add_ltr("ltr3_tfchip_tf_list", "3\u2032 LTR")

  if (length(rows) == 0) return(data.frame())
  unique(do.call(rbind, rows))
}

.dashboard_profile_tfchip_plot <- function(anno, top_n = 15) {
  long <- .dashboard_profile_tfchip_data(anno)
  if (nrow(long) == 0) {
    return(htmltools::HTML("<p><em>No LTR TF ChIP-seq TF overlaps were found for this profile.</em></p>"))
  }

  any_ltr <- unique(long[, c("HERV_id", "tf"), drop = FALSE])
  totals <- as.data.frame(table(any_ltr$tf), stringsAsFactors = FALSE)
  colnames(totals) <- c("tf", "total_hervs")
  totals <- totals[order(totals$total_hervs, decreasing = TRUE), , drop = FALSE]
  selected_tfs <- head(as.character(totals$tf), top_n)

  z <- as.data.frame(
    table(long$tf, long$ltr_position),
    stringsAsFactors = FALSE
  )
  colnames(z) <- c("tf", "ltr_position", "n_hervs")
  z <- z[z$tf %in% selected_tfs & z$n_hervs > 0, , drop = FALSE]
  z$tf <- factor(z$tf, levels = rev(selected_tfs))
  z$label <- .dashboard_wrap_label(as.character(z$tf), width = 34)
  z$hover <- paste0(
    "TF: ", as.character(z$tf),
    "<br>LTR: ", z$ltr_position,
    "<br>HERVs with TF ChIP-seq overlap: ", z$n_hervs
  )

  .dashboard_bar(
    z,
    label_col = "label",
    value_col = "n_hervs",
    axis_title = "Number of HERVs",
    group_col = "ltr_position",
    hover_col = "hover",
    colors = c("5\u2032 LTR" = "#F58518", "3\u2032 LTR" = "#54A24B"),
    max_height = 620
  )
}

.dashboard_profile_tfchip_table <- function(anno, top_n = 20) {
  f <- anno$features
  if (is.null(f) || nrow(f) == 0) return(.dashboard_safe_dt(data.frame()))

  n <- nrow(f)
  ncol0 <- function(col) {
    if (!col %in% colnames(f)) return(rep(0, n))
    .dashboard_numeric(f[[col]], n)
  }

  f$tfchip_total_n_tfs <- ncol0("ltr5_tfchip_n_tfs") + ncol0("ltr3_tfchip_n_tfs")
  f <- f[order(f$tfchip_total_n_tfs, decreasing = TRUE), , drop = FALSE]
  f <- head(f, top_n)

  selected <- c(
    "HERV_id", "subfamily", "tfchip_total_n_tfs",
    "ltr5_tfchip_n_tfs", "ltr5_tfchip_tf_list",
    "ltr3_tfchip_n_tfs", "ltr3_tfchip_tf_list"
  )
  selected <- selected[selected %in% colnames(f)]
  .dashboard_safe_dt(f[, selected, drop = FALSE], page_length = min(10, top_n))
}

.dashboard_profile_tfbm_plot <- function(anno, top_n = 15) {
  z <- anno$top_tfbm_motifs
  if (is.null(z) || nrow(z) == 0) {
    return(htmltools::HTML("<p><em>No detailed TFBM results available.</em></p>"))
  }

  value_col <- if ("n_hervs" %in% colnames(z)) "n_hervs" else "n_hits"
  label_col <- if ("motif_alt_id" %in% colnames(z)) "motif_alt_id" else colnames(z)[1]
  z <- z[order(z[[value_col]], decreasing = TRUE), , drop = FALSE]
  z <- head(z, top_n)
  z$label <- .dashboard_wrap_label(z[[label_col]])
  z$hover <- paste0(z[[label_col]], "<br>HERVs/hits: ", z[[value_col]])
  .dashboard_bar(z, "label", value_col, "Number of HERVs", hover_col = "hover")
}

.dashboard_profile_rbp_plot <- function(anno, top_n = 15) {
  hits <- anno$rbp_hits
  if (is.null(hits) || nrow(hits) == 0 ||
      !all(c("HERV_id", "motif_alt_id") %in% colnames(hits))) {
    return(htmltools::HTML("<p><em>No detailed RBP results available.</em></p>"))
  }

  presence <- unique(hits[, c("HERV_id", "motif_alt_id")])
  z <- as.data.frame(sort(table(presence$motif_alt_id), decreasing = TRUE), stringsAsFactors = FALSE)
  colnames(z) <- c("label", "value")
  z <- head(z, top_n)
  z$hover <- paste0(z$label, "<br>HERVs: ", z$value)
  .dashboard_bar(z, "label", "value", "Number of HERVs", hover_col = "hover")
}

.dashboard_comparison_cards <- function(cmp) {
  input <- cmp$input_summary
  fg_input <- input$n_input_ids[input$group == cmp$foreground_name][1]
  bg_input <- input$n_input_ids[input$group == cmp$background_name][1]
  fg_match <- input$n_matched_ids[input$group == cmp$foreground_name][1]
  bg_match <- input$n_matched_ids[input$group == cmp$background_name][1]
  missing_n <- if (!is.null(cmp$annotation$missing_ids)) length(cmp$annotation$missing_ids) else 0

  cards <- data.frame(
    label = c(
      "Foreground input", "Background input", "Matched foreground",
      "Matched background", "Missing IDs", "Binary features tested",
      "Numeric features tested", "FDR-significant key results"
    ),
    value = c(
      fg_input, bg_input, fg_match, bg_match, missing_n,
      if (is.null(cmp$binary_features)) 0 else nrow(cmp$binary_features),
      if (is.null(cmp$numeric_features)) 0 else nrow(cmp$numeric_features),
      if (is.null(cmp$key_results)) 0 else sum(cmp$key_results$significant, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )

  .dashboard_cards(cards)
}

.dashboard_comparison_key_table <- function(cmp, top_n = 30) {
  x <- cmp$key_results
  if (is.null(x) || nrow(x) == 0) return(.dashboard_safe_dt(data.frame()))
  x <- head(x, top_n)
  .dashboard_safe_dt(x, page_length = min(15, top_n))
}

.dashboard_comparison_regulatory_table <- function(cmp, top_n = 30) {
  x <- cmp$key_results
  if (is.null(x) || nrow(x) == 0) return(.dashboard_safe_dt(data.frame()))
  x <- x[x$layer %in% c("external_evidence", "tfchip"), , drop = FALSE]
  x <- head(x, top_n)
  .dashboard_safe_dt(x, page_length = min(15, top_n))
}

.dashboard_comparison_features_preview <- function(cmp, top_n = 20) {
  f <- cmp$features
  selected <- c(
    "HERV_id", "comparison_group", "subfamily", "domain_count", "domains_type",
    "has_ltr5", "has_ltr3", "ltr5_tfbm_burden", "ltr3_tfbm_burden",
    "rbp_burden", "internal_encode_detected_layers",
    "ltr5_encode_detected_layers", "ltr3_encode_detected_layers",
    "ltr5_tfchip_tf_list", "ltr3_tfchip_tf_list"
  )
  selected <- selected[selected %in% colnames(f)]
  .dashboard_safe_dt(head(f[, selected, drop = FALSE], top_n), page_length = min(10, top_n))
}

.dashboard_output_files <- function(x) {
  files <- x$output_files
  if (is.null(files)) return(.dashboard_safe_dt(data.frame()))

  values <- unlist(files, use.names = TRUE)
  values <- values[!is.na(values) & values != ""]
  if (length(values) == 0) return(.dashboard_safe_dt(data.frame()))

  .dashboard_safe_dt(
    data.frame(output = names(values), path = unname(values), stringsAsFactors = FALSE),
    page_length = length(values),
    filter = "none"
  )
}

.dashboard_has_rows <- function(x) {
  is.data.frame(x) && nrow(x) > 0
}

#' Generate an interactive HERVariumR dashboard
#'
#' Generates a concise HTML dashboard from a `HERVarium_annotation` or
#' `HERVarium_comparison` object. Temporary R Markdown source files are rendered
#' outside the user output directory and removed automatically.
#'
#' @param x A HERVarium_annotation or HERVarium_comparison object.
#' @param output_dir Output directory.
#' @param report_file Name of the HTML report.
#' @param title Report title.
#' @param top_n Number of top features shown in dashboard plots.
#' @param self_contained Whether to produce a self-contained HTML file.
#'
#' @return Invisibly returns the generated HTML path.
#' @export
generate_hervarium_dashboard <- function(x,
                                         output_dir = NULL,
                                         report_file = NULL,
                                         title = NULL,
                                         top_n = 25,
                                         self_contained = TRUE) {
  required <- c("rmarkdown", "DT", "plotly", "htmltools")
  missing_packages <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop(
      "The following packages are required for dashboards: ",
      paste(missing_packages, collapse = ", ")
    )
  }

  if (!inherits(x, "HERVarium_comparison") &&
      !inherits(x, "HERVarium_annotation")) {
    stop("x must be a HERVarium_comparison or HERVarium_annotation object.")
  }

  if (is.null(output_dir)) output_dir <- x$output_dir
  if (is.null(output_dir)) stop("Please provide output_dir or use an object with output_dir.")
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

  object_type <- if (inherits(x, "HERVarium_comparison")) "comparison" else "annotation"

  # Backward-compatible dashboard regeneration for comparison objects created
  # before the unified key-results table was introduced.
  if (object_type == "comparison" &&
      (is.null(x$key_results) || !is.data.frame(x$key_results))) {
    x$key_results <- build_comparison_key_results(
      binary_results = x$binary_features,
      numeric_results = x$numeric_features,
      subfamily_results = x$subfamily_enrichment,
      domain_type_results = x$domain_type_enrichment,
      tfbm_ltr5_results = x$ltr5_tfbm_motif_enrichment,
      tfbm_ltr3_results = x$ltr3_tfbm_motif_enrichment,
      rbp_motif_results = x$rbp_motif_enrichment,
      tfchip_ltr5_results = x$ltr5_tfchip_tf_enrichment,
      tfchip_ltr3_results = x$ltr3_tfchip_tf_enrichment,
      detailed_top_n = min(top_n, 15)
    )
  }

  if (is.null(report_file)) {
    report_file <- if (object_type == "comparison") {
      "hervarium_comparison_dashboard.html"
    } else {
      "hervarium_annotation_dashboard.html"
    }
  }

  if (is.null(title)) {
    title <- if (object_type == "comparison") {
      paste0("HERVariumR comparison: ", x$foreground_name, " vs ", x$background_name)
    } else {
      "HERVariumR annotation dashboard"
    }
  }

  work_dir <- tempfile("hervarium_dashboard_")
  dir.create(work_dir, recursive = TRUE)
  on.exit(unlink(work_dir, recursive = TRUE, force = TRUE), add = TRUE)

  object_rds <- file.path(work_dir, "hervarium_dashboard_object.rds")
  rmd_file <- file.path(work_dir, "hervarium_dashboard.Rmd")
  saveRDS(x, object_rds)

  if (object_type == "comparison") {
    write_hervarium_comparison_dashboard_rmd(
      rmd_file = rmd_file,
      title = title,
      top_n = top_n,
      self_contained = self_contained,
      x = x
    )
  } else {
    write_hervarium_annotation_dashboard_rmd(
      rmd_file = rmd_file,
      title = title,
      top_n = top_n,
      self_contained = self_contained,
      x = x
    )
  }

  render_env <- new.env(parent = environment(generate_hervarium_dashboard))
  output_path <- rmarkdown::render(
    input = rmd_file,
    output_file = report_file,
    output_dir = normalizePath(output_dir, mustWork = TRUE),
    quiet = TRUE,
    clean = TRUE,
    knit_root_dir = work_dir,
    envir = render_env
  )

  message("HERVariumR dashboard written to: ", output_path)
  invisible(output_path)
}

.dashboard_yaml <- function(title, self_contained = TRUE) {
  safe_title <- gsub('"', "'", title, fixed = TRUE)
  c(
    "---",
    paste0('title: "', safe_title, '"'),
    "output:",
    "  html_document:",
    "    toc: true",
    "    toc_float:",
    "      collapsed: false",
    "      smooth_scroll: true",
    "    theme: flatly",
    paste0("    self_contained: ", ifelse(self_contained, "true", "false")),
    "---",
    ""
  )
}

.dashboard_setup_lines <- function(top_n) {
  c(
    "```{r setup, include=FALSE}",
    "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE, out.width = '100%')",
    "obj <- readRDS('hervarium_dashboard_object.rds')",
    paste0("top_n <- ", as.integer(top_n)),
    "plot_top_n <- min(top_n, 15)",
    "```",
    "",
    "<style>",
    ".main-container { max-width: 1500px !important; width: 96% !important; }",
    "body { font-size: 15px; }",
    "h1 { margin-top: 1.2em; }",
    "h2 { margin-top: 1.0em; }",
    ".tocify { width: 250px !important; }",
    ".toc-content { padding-left: 20px !important; }",
    ".dataTables_wrapper { font-size: 13px; }",
    ".plotly.html-widget, .js-plotly-plot { width: 100% !important; margin-bottom: 0.25em !important; }",
    ".html-widget { margin-bottom: 0.25em !important; }",
    ".section.level2, .section.level3 { margin-bottom: 0.55em !important; }",
    ".herv-card-grid { display: flex; flex-wrap: wrap; gap: 12px; margin: 10px 0 18px 0; }",
    ".herv-card { flex: 1 1 160px; min-width: 145px; padding: 14px 16px; border: 1px solid #E1E5EA; border-radius: 8px; background: #F8FAFC; }",
    ".herv-card-value { font-size: 24px; font-weight: 700; line-height: 1.1; }",
    ".herv-card-label { margin-top: 6px; color: #555; font-size: 13px; }",
    ".interpretation-note { padding: 11px 13px; border-left: 5px solid #4C78A8; background: #F3F7FB; margin: 8px 0 14px 0; }",
    "</style>",
    ""
  )
}

write_hervarium_annotation_dashboard_rmd <- function(rmd_file,
                                                     title,
                                                     top_n = 25,
                                                     self_contained = TRUE,
                                                     x = NULL) {
  lines <- c(
    .dashboard_yaml(title, self_contained),
    .dashboard_setup_lines(top_n),
    "# Overview",
    "",
    "```{r overview-cards}",
    ".dashboard_profile_cards(obj)",
    "```",
    "",
    "<div class='interpretation-note'><strong>Candidate ranking:</strong> the preview below uses a transparent heuristic based on detected regulatory layers, LTR cCRE/TF ChIP-seq evidence, domains and LTR preservation. It is an inspection aid, not a statistical score.</div>",
    "",
    "## Top candidate HERVs",
    "",
    "```{r candidate-table}",
    ".dashboard_profile_candidate_table(obj, top_n = min(top_n, 20))",
    "```",
    "",
    "# HERV identity and structure",
    "",
    "## Top subfamilies",
    "",
    "```{r subfamilies}",
    ".dashboard_profile_subfamilies(obj, top_n = plot_top_n)",
    "```",
    "",
    "## Internal domain types",
    "",
    "```{r domains}",
    ".dashboard_profile_domains(obj)",
    "```",
    "",
    "## LTR structure",
    "",
    "```{r ltr-structure}",
    ".dashboard_profile_ltr_structure(obj)",
    "```",
    "",
    "# Transcript context",
    "",
    "```{r transcript-context}",
    ".dashboard_profile_context(obj)",
    "```",
    "",
    "# Regulatory context",
    "",
    "<div class='interpretation-note'>ENCODE/UCSC signal, cCRE overlap and TF ChIP-seq peak-cluster overlap are contextual evidence. They do not by themselves demonstrate autonomous HERV transcription or target-gene regulation.</div>",
    "",
    "## Evidence overview",
    "",
    "```{r regulatory-overview}",
    ".dashboard_profile_regulatory(obj)",
    "```",
    "",
    "## Top HERVs by regulatory evidence",
    "",
    "```{r regulatory-table}",
    ".dashboard_profile_regulatory_table(obj, top_n = min(top_n, 20))",
    "```",
    "",
    "# LTR TF ChIP-seq evidence",
    "",
    "<div class='interpretation-note'>This layer summarizes experimental TF ChIP-seq peak-cluster overlaps with the associated 5\u2032 and 3\u2032 LTR intervals. It is complementary to sequence-based TFBM predictions and does not by itself prove direct TF recruitment by the HERV sequence.</div>",
    "",
    "## Most frequent overlapping TFs",
    "",
    "```{r tfchip-profile-plot}",
    ".dashboard_profile_tfchip_plot(obj, top_n = plot_top_n)",
    "```",
    "",
    "## HERV-level TF ChIP-seq lists",
    "",
    "```{r tfchip-profile-table}",
    ".dashboard_profile_tfchip_table(obj, top_n = min(top_n, 20))",
    "```",
    ""
  )

  if (!is.null(x) && .dashboard_has_rows(x$top_tfbm_motifs)) {
    lines <- c(
      lines,
      "# Detailed LTR TFBM layer",
      "",
      "```{r detailed-tfbm}",
      ".dashboard_profile_tfbm_plot(obj, top_n = plot_top_n)",
      "```",
      "",
      "```{r detailed-tfbm-table}",
      "tfbm_cols <- c('HERV_id', 'ltr5_n_unique_tfbm_tf_names_detailed', 'ltr5_tfbm_tf_names_all', 'ltr3_n_unique_tfbm_tf_names_detailed', 'ltr3_tfbm_tf_names_all')",
      "tfbm_cols <- tfbm_cols[tfbm_cols %in% colnames(obj$features)]",
      ".dashboard_safe_dt(obj$features[, tfbm_cols, drop = FALSE], page_length = 10)",
      "```",
      ""
    )
  }

  if (!is.null(x) && .dashboard_has_rows(x$rbp_hits)) {
    lines <- c(
      lines,
      "# Detailed RBP layer",
      "",
      "```{r detailed-rbp}",
      ".dashboard_profile_rbp_plot(obj, top_n = plot_top_n)",
      "```",
      "",
      "```{r detailed-rbp-table}",
      "rbp_cols <- c('HERV_id', 'n_unique_rbp_names_detailed', 'rbp_names_all', 'top_rbp_names_preview')",
      "rbp_cols <- rbp_cols[rbp_cols %in% colnames(obj$features)]",
      ".dashboard_safe_dt(obj$features[, rbp_cols, drop = FALSE], page_length = 10)",
      "```",
      ""
    )
  }

  lines <- c(
    lines,
    "# Main output files",
    "",
    "```{r output-files}",
    ".dashboard_output_files(obj)",
    "```"
  )

  writeLines(lines, rmd_file)
}

write_hervarium_comparison_dashboard_rmd <- function(rmd_file,
                                                     title,
                                                     top_n = 25,
                                                     self_contained = TRUE,
                                                     x = NULL) {
  has_detailed <- !is.null(x) && !is.null(x$key_results) &&
    any(x$key_results$layer %in% c("tfchip_tf", "tfbm", "rbp"))

  lines <- c(
    .dashboard_yaml(title, self_contained),
    .dashboard_setup_lines(top_n),
    "# Overview",
    "",
    "```{r overview-cards}",
    ".dashboard_comparison_cards(obj)",
    "```",
    "",
    "## Top key results",
    "",
    "```{r key-results-table}",
    ".dashboard_comparison_key_table(obj, top_n = 30)",
    "```",
    "",
    "# Main feature differences",
    "",
    "## Top binary enrichments",
    "",
    "```{r binary-enrichment}",
    ".dashboard_binary_plot(obj$binary_features, top_n = plot_top_n)",
    "```",
    "",
    "## Top numeric shifts",
    "",
    "```{r numeric-shifts}",
    ".dashboard_numeric_plot(obj$numeric_features, top_n = plot_top_n)",
    "```",
    "",
    "# Regulatory evidence",
    "",
    "<div class='interpretation-note'>Positive effects indicate higher prevalence or burden in the foreground group. External regulatory annotations provide contextual evidence rather than proof of autonomous HERV activity.</div>",
    "",
    "## Binary regulatory enrichments",
    "",
    "```{r regulatory-binary}",
    "reg_binary <- obj$key_results[obj$key_results$layer %in% c('external_evidence', 'tfchip') & obj$key_results$effect_type == 'log2_odds_ratio', , drop = FALSE]",
    ".dashboard_key_plot(reg_binary, top_n = plot_top_n, axis_title = 'log2 odds ratio')",
    "```",
    "",
    "## Numeric regulatory shifts",
    "",
    "```{r regulatory-numeric}",
    "reg_numeric <- obj$key_results[obj$key_results$layer %in% c('external_evidence', 'tfchip') & obj$key_results$effect_type == 'delta_median', , drop = FALSE]",
    ".dashboard_key_plot(reg_numeric, top_n = plot_top_n, axis_title = 'Delta median')",
    "```",
    "",
    "## Regulatory key results",
    "",
    "```{r regulatory-table}",
    ".dashboard_comparison_regulatory_table(obj, top_n = 30)",
    "```",
    "",
    "# HERV identity and domains",
    "",
    "```{r identity-domains}",
    "identity_results <- obj$key_results[obj$key_results$layer %in% c('subfamily', 'internal_domains', 'ltr_structure', 'transcript_context', 'terminal_exon', 'ifn_ltr') & obj$key_results$effect_type == 'log2_odds_ratio', , drop = FALSE]",
    ".dashboard_key_plot(identity_results, top_n = plot_top_n, axis_title = 'log2 odds ratio')",
    "```",
    ""
  )

  if (has_detailed) {
    lines <- c(
      lines,
      "# Detailed regulatory enrichments",
      "",
      "```{r detailed-regulatory}",
      "detail_results <- obj$key_results[obj$key_results$layer %in% c('tfchip_tf', 'tfbm', 'rbp'), , drop = FALSE]",
      ".dashboard_key_plot(detail_results, top_n = plot_top_n, axis_title = 'log2 odds ratio')",
      "```",
      "",
      "```{r detailed-regulatory-table}",
      ".dashboard_safe_dt(head(detail_results, 50), page_length = 15)",
      "```",
      ""
    )
  }

  lines <- c(
    lines,
    "# Matched HERV preview",
    "",
    "```{r feature-preview}",
    ".dashboard_comparison_features_preview(obj, top_n = 20)",
    "```",
    "",
    "# Main output files",
    "",
    "```{r output-files}",
    ".dashboard_output_files(obj)",
    "```"
  )

  writeLines(lines, rmd_file)
}
