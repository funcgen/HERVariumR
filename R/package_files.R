#' Get path to a bundled HERVariumR annotation file
#'
#' Returns the installed path to a file bundled inside `inst/extdata`.
#'
#' @param filename Name of a file inside the package `extdata` directory.
#'
#' @return Full path to the requested bundled file.
#' @export
hervarium_file <- function(filename) {
  path <- system.file("extdata", filename, package = "HERVariumR")
  
  if (path == "") {
    stop(
      "Could not find bundled HERVariumR file: ", filename,
      call. = FALSE
    )
  }
  
  path
}


.get_default_annotation_file <- function(annotation_file = NULL) {
  if (!is.null(annotation_file)) {
    return(annotation_file)
  }
  
  hervarium_file("transcript_context.with_herv_id.tsv.gz")
}


.get_default_ifn_stat1_file <- function(ifn_stat1_file = NULL) {
  if (!is.null(ifn_stat1_file)) {
    return(ifn_stat1_file)
  }
  
  hervarium_file("LTR_IFN_STAT1_summary.tsv")
}


.get_default_ifn_stat1stat2_irf_file <- function(ifn_stat1stat2_irf_file = NULL) {
  if (!is.null(ifn_stat1stat2_irf_file)) {
    return(ifn_stat1stat2_irf_file)
  }
  
  hervarium_file("LTR_IFN_STAT1STAT2_IRF_summary.tsv")
}


.get_default_last_exon_file <- function(last_exon_file = NULL) {
  if (!is.null(last_exon_file)) {
    return(last_exon_file)
  }
  
  hervarium_file("HERV_domains_transcript_context_last_exon.xlsx")
}