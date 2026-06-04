load_transcript_context <- function(annotation_file = NULL) {
  if (is.null(annotation_file)) {
    annotation_file <- system.file(
      "extdata",
      "transcript_context.with_herv_id.tsv.gz",
      package = "HERVariumR"
    )
  }
  
  if (annotation_file == "") {
    stop("Annotation file not found. Check that transcript_context.with_herv_id.tsv.gz exists in inst/extdata/")
  }
  
  read.delim(
    annotation_file,
    sep = "\t",
    header = TRUE,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}