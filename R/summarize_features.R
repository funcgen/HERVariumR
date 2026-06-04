# R/summarize_features.R

summarize_herv_features <- function(x) {
  if (inherits(x, "HERVarium_annotation")) {
    df <- x$summary
  } else {
    df <- x
  }
  
  data.frame(
    HERV_id = df$HERV_id,
    subfamily = sub("-int_.*$", "", df$HERV_id),
    locid = df$locid,
    chrom = df$chrom,
    start = df$start,
    end = df$end,
    strand = df$strand,
    
    has_domain = df$domain_count > 0,
    domain_count = df$domain_count,
    max_domain_cov = df$max_domain_cov,
    domains_type = df$domains_type,
    
    
    has_gag = grepl("gag", df$domains_type, ignore.case = TRUE),
    has_pol = grepl("pol", df$domains_type, ignore.case = TRUE),
    has_env = grepl("env", df$domains_type, ignore.case = TRUE),
    has_accessory = grepl("accessory", df$domains_type, ignore.case = TRUE),
    
    has_complete_gag_pol_env =
      grepl("gag", df$domains_type, ignore.case = TRUE) &
      grepl("pol", df$domains_type, ignore.case = TRUE) &
      grepl("env", df$domains_type, ignore.case = TRUE),
    
    has_ltr5 = df$has_ltr5 == 1,
    has_ltr3 = df$has_ltr3 == 1,
    has_both_ltrs = df$has_ltr5 == 1 & df$has_ltr3 == 1,
    
    ltr5_name = df$ltr5_name,
    ltr3_name = df$ltr3_name,
    ltr5_tfbm_burden = df$ltr5_tfbm_burden,
    ltr3_tfbm_burden = df$ltr3_tfbm_burden,
    total_ltr_tfbm_burden = df$ltr5_tfbm_burden + df$ltr3_tfbm_burden,
    
    rbp_burden = df$rbp_burden,
    rbp_unique = df$rbp_unique,
    
    feature_overlap = df$feature_overlap,
    ov_gene_types = df$ov_gene_types,
    
    is_intergenic = grepl("intergenic", df$feature_overlap, ignore.case = TRUE),
    is_intergenic_same_strand = grepl("intergenic", df$feature_overlap_same_strand, ignore.case = TRUE),
    
    overlaps_gene = grepl("gene", df$feature_overlap, ignore.case = TRUE),
    overlaps_exon = grepl("exon", df$feature_overlap, ignore.case = TRUE),
    overlaps_intron = grepl("intron", df$feature_overlap, ignore.case = TRUE),
    overlaps_cds = grepl("cds", df$feature_overlap, ignore.case = TRUE),
    overlaps_utr = grepl("utr", df$feature_overlap, ignore.case = TRUE),
    
    overlaps_lncRNA = grepl("lncRNA", df$ov_gene_types, ignore.case = TRUE),
    overlaps_protein_coding = grepl("protein_coding", df$ov_gene_types, ignore.case = TRUE),
    
    stringsAsFactors = FALSE
  )
}