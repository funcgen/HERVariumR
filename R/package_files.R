#' Get path to a bundled HERVariumR file
#'
#' Returns the installed path to a file bundled inside the HERVariumR package.
#'
#' @param filename Name of the bundled file.
#' @param subdir Subdirectory inside `inst/`. Defaults to `"extdata"`.
#'
#' @return Full path to the requested bundled file.
#' @export
hervarium_file <- function(filename, subdir = "extdata") {
  path <- system.file(subdir, filename, package = "HERVariumR")

  if (path == "") {
    stop(
      "Could not find bundled HERVariumR file: ", filename,
      " in inst/", subdir, "/",
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

.get_default_internal_encode_file <- function(internal_encode_file = NULL) {
  if (!is.null(internal_encode_file)) {
    return(internal_encode_file)
  }

  hervarium_file("internal_encode_44organs_compact.tsv.gz")
}


.get_default_ltr_encode_file <- function(ltr_encode_file = NULL) {
  if (!is.null(ltr_encode_file)) {
    return(ltr_encode_file)
  }

  hervarium_file("ltr_encode_44organs_compact.tsv.gz")
}
.get_default_tfchip_file <- function(tfchip_file = NULL) {
  if (!is.null(tfchip_file)) {
    return(tfchip_file)
  }

  hervarium_file("ltr_tfchip_compact.tsv.gz")
}


.clean_herv_id_vector <- function(x) {
  z <- trimws(as.character(x))

  invalid <- is.na(z) |
    z == "" |
    z == "." |
    toupper(z) == "NA"

  out <- unique(z[!invalid])

  attr(out, "n_invalid_removed") <- sum(invalid)
  attr(out, "n_duplicates_removed") <- sum(!invalid) - length(out)

  out
}
