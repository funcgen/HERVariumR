#' Profile a single list of HERV IDs
#'
#' Annotates and summarizes a single list of HERV internal-region IDs using
#' HERVariumR annotation layers. The function can add internal domain
#' annotations, LTR information, IFN-related LTR features, terminal-exon domain
#' annotations, optional detailed LTR TFBM hits, optional detailed RBP motif
#' hits, summary plots, and an interactive dashboard.
#'
#' @param herv_ids Character vector of HERV IDs to annotate.
#' @param annotation_file Optional path to the main transcript-context annotation file. If `NULL`, the bundled HERVariumR annotation file is used.
#' @param output_dir Directory where output tables, plots, and reports will be written.
#' @param ifn_stat1_file Optional path to the STAT1 LTR summary file. If `NULL`, the bundled HERVariumR file is used.
#' @param ifn_stat1stat2_irf_file Optional path to the STAT1/STAT2/IRF LTR summary file. If `NULL`, the bundled HERVariumR file is used.
#' @param last_exon_file Optional path to the terminal-exon domain annotation file. If `NULL`, the bundled HERVariumR file is used.
#' @param domain_coverage_cutoff Numeric. Minimum HMM profile coverage required for an internal retroviral domain to be considered valid.
#' @param keep_raw_domain_columns Logical. Whether to keep unfiltered domain annotation columns with `_raw` suffixes.
#' @param add_tfbm_details Logical. Whether to add detailed LTR TFBM annotations from a FIMO file.
#' @param fimo_file Path to the detailed LTR TFBM FIMO file. This large file is not bundled with the package and must be provided by the user when `add_tfbm_details = TRUE`.
#' @param qvalue_cutoff Numeric. Maximum q-value allowed for detailed LTR TFBM hits.
#' @param add_rbp_details Logical. Whether to add detailed RBP motif annotations from a FIMO-like file.
#' @param rbp_file Path to the detailed RBP FIMO file. This large file is not bundled with the package and must be provided by the user when `add_rbp_details = TRUE`.
#' @param rbp_qvalue_cutoff Numeric. Maximum q-value allowed for detailed RBP motif hits.
#' @param use_awk Logical. Whether to use awk for fast filtering of large FIMO files.
#' @param make_plots Logical. Whether to generate summary plots.
#' @param make_ltr_tfbm_plots Logical. Whether to generate detailed LTR TFBM plots when `add_tfbm_details = TRUE`.
#' @param make_internal_domain_plot Logical. Whether to generate the internal-domain conservation bubble plot.
#' @param top_n Integer. Number of top categories to show in general plots.
#' @param top_n_motifs Integer. Number of top TF/RBP motifs to show in motif plots.
#' @param top_n_subfamilies Integer. Number of top HERV subfamilies to show in subfamily-level plots.
#' @param internal_domain_size_by Character. Bubble-size mode for the internal domain plot. One of `"fixed"` or `"domain_count"`.
#' @param verbose Logical. Whether to print progress messages.
#' @param show_progress Logical. Whether to show progress bars for long-running steps.
#'
#' @return A `HERVarium_annotation` object containing annotation tables, feature summaries, optional regulatory layers, and output metadata.
#' @export

profile_hervs <- function(herv_ids,
                          annotation_file = NULL,
                          output_dir = "HERVariumR_results",
                          
                          # Optional annotation layers
                          ifn_stat1_file = NULL,
                          ifn_stat1stat2_irf_file = NULL,
                          last_exon_file = NULL,
                          
                          # Domain filtering
                          domain_coverage_cutoff = 0.40,
                          keep_raw_domain_columns = TRUE,
                          
                          # Optional heavy TFBM layer
                          add_tfbm_details = FALSE,
                          fimo_file = NULL,
                          qvalue_cutoff = 1,
                          
                          # Optional detailed RBP layer
                          add_rbp_details = FALSE,
                          rbp_file = NULL,
                          rbp_qvalue_cutoff = 1,
                          
                          use_awk = TRUE,
                          
                          # Plotting options
                          make_plots = TRUE,
                          make_ltr_tfbm_plots = TRUE,
                          make_internal_domain_plot = TRUE,
                          top_n = 25,
                          top_n_motifs = 25,
                          top_n_subfamilies = 20,
                          internal_domain_size_by = c("fixed", "domain_count"),
                          
                          # Verbosity
                          verbose = TRUE,
                          show_progress = TRUE) {
  
  internal_domain_size_by <- match.arg(internal_domain_size_by)
  
  if (missing(herv_ids) || length(herv_ids) == 0) {
    stop("Please provide a vector of HERV IDs.")
  }
  
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  if (verbose) {
    message("Starting HERVariumR single-list profiling")
    message("Output directory: ", output_dir)
  }
  
  # ------------------------------------------------------------
  # 1. Core annotation
  # ------------------------------------------------------------
  
  if (verbose) {
    message("\n[1/5] Annotating HERV list")
  }
  
  anno <- annotate_hervs(
    herv_ids = herv_ids,
    annotation_file = annotation_file,
    ifn_stat1_file = ifn_stat1_file,
    ifn_stat1stat2_irf_file = ifn_stat1stat2_irf_file,
    last_exon_file = last_exon_file,
    domain_coverage_cutoff = domain_coverage_cutoff,
    keep_raw_domain_columns = keep_raw_domain_columns,
    output_dir = output_dir,
    verbose = verbose
  )
  
  # ------------------------------------------------------------
  # 2. Basic plots
  # ------------------------------------------------------------
  
  if (make_plots) {
    if (verbose) {
      message("\n[2/5] Generating summary annotation plots")
    }
    
    plot_herv_annotation(
      anno,
      output_dir = output_dir,
      top_n = top_n
    )
  }
  
  # ------------------------------------------------------------
  # 3. Internal domain conservation plot
  # ------------------------------------------------------------
  
  if (make_plots && make_internal_domain_plot) {
    if (verbose) {
      message("\n[3/5] Generating internal-domain conservation plot")
    }
    
    plot_internal_domain_bubbles(
      anno,
      output_dir = output_dir,
      top_n_subfamilies = top_n_subfamilies,
      top_n_domains = top_n,
      size_by = internal_domain_size_by
    )
  }
  
  # ------------------------------------------------------------
  # 4. Optional detailed LTR TFBM layer
  # ------------------------------------------------------------
  
  if (add_tfbm_details) {
    
    if (is.null(fimo_file)) {
      stop(
        "add_tfbm_details = TRUE requires a fimo_file. ",
        "Provide the path to fimo_parsed_v4.tsv or equivalent."
      )
    }
    
    if (verbose) {
      message("\n[4/5] Adding detailed LTR TFBM layer")
    }
    
    anno <- add_ltr_tfbm_details(
      anno,
      fimo_file = fimo_file,
      output_dir = output_dir,
      qvalue_cutoff = qvalue_cutoff,
      top_n_motifs = top_n_motifs,
      use_awk = use_awk,
      show_progress = show_progress
    )
    
    if (make_plots && make_ltr_tfbm_plots) {
      if (verbose) {
        message("\nGenerating detailed LTR TFBM plots")
      }
      
      plot_ltr_tfbm(
        anno,
        output_dir = output_dir,
        top_n_motifs = top_n_motifs,
        top_n_subfamilies = top_n_subfamilies
      )
    }
    
  } else {
    
    if (verbose) {
      message("\n[4/5] Skipping detailed LTR TFBM layer")
      message("Set add_tfbm_details = TRUE and provide fimo_file to enable it.")
    }
  }
  
  # ------------------------------------------------------------
  # Optional detailed RBP layer
  # ------------------------------------------------------------
  
  if (add_rbp_details) {
    
    if (is.null(rbp_file)) {
      stop(
        "add_rbp_details = TRUE requires a rbp_file. ",
        "Provide the path to RBP_fimo.tsv or equivalent."
      )
    }
    
    if (!exists("add_rbp_details", mode = "function")) {
      stop(
        "Function add_rbp_details() is not available. ",
        "Please source R/add_rbp_layer.R before running profile_hervs()."
      )
    }
    
    if (verbose) {
      message("\n[5/5]Adding detailed RBP layer")
    }
    
    anno <- add_rbp_details(
      anno,
      rbp_file = rbp_file,
      output_dir = output_dir,
      qvalue_cutoff = rbp_qvalue_cutoff,
      use_awk = use_awk,
      show_progress = show_progress
    )
  }
  
  # ------------------------------------------------------------
  # 5. Final message
  # ------------------------------------------------------------
  
  if (verbose) {
    message("\nHERVariumR profiling finished.")
    message("Results written to: ", output_dir)
  }
  
  return(anno)
}