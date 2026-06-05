# HERVariumR

**HERVariumR** is an R package for fast functional annotation, profiling, and comparison of human endogenous retrovirus (HERV) locus lists.

The package is designed for researchers who quantify HERV expression at the locus level and obtain lists of HERVs of interest, such as differentially expressed HERVs, condition-specific HERVs, interferon-responsive HERVs, or HERVs associated with a given cell type or disease state.

HERVariumR converts those lists into interpretable biological summaries by integrating multiple annotation layers:

* internal retroviral protein-domain annotations;
* Gag, Pol, Env and accessory-domain classification;
* HERV subfamily information;
* 5′ and 3′ LTR annotation;
* LTR transcription-factor-binding motif burden;
* interferon-related STAT1, STAT1::STAT2 and IRF motif summaries;
* optional detailed LTR TFBM motif enrichment;
* optional RNA-binding protein motif annotation;
* transcript-context annotation, including intronic, exonic, protein-coding and lncRNA overlap;
* terminal-exon domain annotations;
* static plots and interactive HTML dashboards.

The main philosophy of HERVariumR is simple: **start from a list of HERV IDs and obtain a compact, publication-ready biological interpretation**.

---

## Installation

HERVariumR can be installed directly from GitHub:

```r
install.packages("remotes")
remotes::install_github("funcgen/HERVariumR")
```

Load the package:

```r
library(HERVariumR)
packageVersion("HERVariumR")
```

---

## Main use cases

HERVariumR currently supports two common workflows.

### 1. Profile a single HERV list

Use this when you have one list of HERVs and want to know what kind of loci they are.

Examples:

* HERVs upregulated in one condition;
* HERVs expressed in one cell type;
* HERVs overlapping accessible chromatin;
* HERVs associated with a phenotype;
* candidate HERV loci selected from another analysis.

The main function is:

```r
profile_hervs()
```

This returns a `HERVarium_annotation` object.

---

### 2. Compare a foreground HERV list against a background

Use this when you want to know whether a HERV list is enriched for specific features compared with a relevant background universe.

Examples:

* upregulated HERVs versus all tested HERVs;
* downregulated HERVs versus all tested HERVs;
* disease-associated HERVs versus detected HERVs;
* HERVs in one cluster versus all expressed HERVs.

The main function is:

```r
compare_herv_lists()
```

This returns a `HERVarium_comparison` object.

The foreground IDs are removed from the background internally, so the final comparison is performed between mutually exclusive groups.

---

## Bundled annotation resources

HERVariumR includes several annotation files in `inst/extdata`. You can retrieve their installed paths with:

```r
hervarium_file("transcript_context.with_herv_id.tsv.gz")
hervarium_file("LTR_IFN_STAT1_summary.tsv")
hervarium_file("LTR_IFN_STAT1STAT2_IRF_summary.tsv")
hervarium_file("HERV_domains_transcript_context_last_exon.xlsx")
```

These bundled resources are used by default by the main functions.

Large motif-level files, such as genome-wide FIMO results for LTR transcription-factor-binding motifs or RBP motifs, are not bundled by default because of their size. They can be supplied manually when running the detailed regulatory layers.

---

## Quick start: profile a list of HERVs

The package includes an example differential-expression table in `sample_data`.

```r
library(HERVariumR)

deg_file <- system.file(
  "sample_data",
  "HERV_DESeq2_DEG_all.tsv",
  package = "HERVariumR"
)

deg <- read.delim(deg_file, stringsAsFactors = FALSE)

head(deg)
colnames(deg)
```

For example, define a list of significantly upregulated HERVs:

```r
test_ids <- deg$transcript[
  deg$log2FoldChange > 0 &
    deg$padj < 0.05
]
```

Run a simple HERVariumR profile:

```r
res_cALD_BBB <- profile_hervs(
  herv_ids = test_ids,
  output_dir = "hervarium_test_results/simple",
  internal_domain_size_by = "domain_count"
)
```

Generate an interactive HTML dashboard:

```r
generate_hervarium_dashboard(
  res_cALD_BBB,
  top_n = 25
)
```

This produces a `HERVarium_annotation` object containing the matched annotation table, feature summaries, missing IDs, summary statistics, output directory information and optional regulatory layers.

---

## Full profiling workflow with detailed TFBM and RBP layers

If detailed motif-level files are available, HERVariumR can add LTR TFBM and RBP motif summaries.

```r
res_cALD_BBB_full <- profile_hervs(
  herv_ids = test_ids,
  domain_coverage_cutoff = 0,
  output_dir = "hervarium_test_results/with_tfbm_and_rbp",

  # Detailed regulatory layers
  add_tfbm_details = TRUE,
  fimo_file = "HERVariumR_large_data/fimo_parsed_v4.tsv",
  add_rbp_details = TRUE,
  rbp_file = "HERVariumR_large_data/RBP_fimo.tsv",
  qvalue_cutoff = 1,
  rbp_qvalue_cutoff = 1,
  use_awk = TRUE,

  # Plotting
  internal_domain_size_by = "domain_count"
)

generate_hervarium_dashboard(
  res_cALD_BBB_full,
  top_n = 25
)
```

The `use_awk = TRUE` option is recommended for large FIMO files because it avoids loading the full file into memory before filtering.

---

## Compare a foreground HERV list against a background

A typical transcriptomic use case is to compare significantly upregulated HERVs against all HERVs tested in the differential-expression analysis.

```r
test_ids <- deg$transcript[
  deg$log2FoldChange > 0 &
    deg$padj < 0.05
]

test_ids_background <- deg$transcript
```

Run the comparison:

```r
cmp_cALD_BBB_full <- compare_herv_lists(
  foreground_ids = test_ids,
  background_ids = test_ids_background,
  foreground_name = "cALD_BBB_sig_up",
  background_name = "tested_non_up",
  output_dir = "hervarium_test_results/comparison_with_tfbm_and_rbp",

  # Domain validity threshold
  domain_coverage_cutoff = 0,

  # Detailed regulatory layers
  add_tfbm_details = TRUE,
  fimo_file = "HERVariumR_large_data/fimo_parsed_v4.tsv",
  add_rbp_details = TRUE,
  rbp_file = "HERVariumR_large_data/RBP_fimo.tsv",
  qvalue_cutoff = 1,
  rbp_qvalue_cutoff = 1,

  # Plots
  make_plots = TRUE,
  plot_top_n = 25,
  plot_all_binary_features = TRUE
)

generate_hervarium_dashboard(
  cmp_cALD_BBB_full,
  top_n = 25
)
```

Although the background is supplied as all tested HERVs, HERVariumR removes the foreground IDs from the background internally. Therefore, the actual comparison is:

```text
significantly upregulated HERVs
versus
tested but not significantly upregulated HERVs
```

This is the recommended setup for enrichment-style analyses of differential-expression results.

---

## Output objects

### `HERVarium_annotation`

Returned by:

```r
profile_hervs()
annotate_hervs()
```

Typical elements include:

```r
names(res_cALD_BBB_full)
```

Common slots include:

* `summary`: matched HERV-level annotation table;
* `features`: compact feature matrix used for summaries and comparisons;
* `missing_ids`: input IDs not found in the annotation;
* `stats`: summary statistics for the input list;
* `ltr_tfbm_hits`: optional detailed LTR TFBM hits;
* `herv_tfbm_summary`: optional HERV-level TFBM summaries;
* `top_tfbm_motifs`: most frequent TFBM motifs;
* `rbp_hits`: optional detailed RBP motif hits;
* `herv_rbp_summary`: optional HERV-level RBP summaries;
* `output_dir`: output directory used by the run;
* `domain_coverage_cutoff`: domain coverage threshold used.

---

### `HERVarium_comparison`

Returned by:

```r
compare_herv_lists()
```

Typical elements include:

```r
names(cmp_cALD_BBB_full)
```

Common slots include:

* `foreground_name`;
* `background_name`;
* `foreground_ids`;
* `background_ids`;
* `features`;
* `input_summary`;
* `binary_features`;
* `numeric_features`;
* `subfamily_enrichment`;
* `domain_hits`;
* `domain_type_enrichment`;
* `ifn_ltr_features`;
* `ltr_tfbm_burden`;
* `rbp_burden`;
* `ltr5_tfbm_motif_enrichment`;
* `ltr3_tfbm_motif_enrichment`;
* `any_ltr_tfbm_motif_enrichment`;
* `rbp_motif_enrichment`;
* `output_dir`.

---

## Interpreting the main comparison tables

### Binary feature enrichment

The `binary_features` table tests whether binary annotations are enriched in the foreground.

Examples of binary features include:

* `has_domain`;
* `has_gag`;
* `has_pol`;
* `has_env`;
* `has_ltr5`;
* `has_ltr3`;
* `has_both_ltrs`;
* `has_any_ifn_related_ltr_motif`;
* `has_terminal_exon_domain_protein_coding`;
* `overlaps_exon`;
* `overlaps_intron`;
* `overlaps_lncRNA`;
* `overlaps_protein_coding`.

The table reports foreground and background counts, percentages, odds ratio, p-value, direction and adjusted p-value.

Positive enrichment means that the feature is more common in the foreground than in the background.

---

### Numeric feature shifts

The `numeric_features` table compares continuous or count-like annotations between foreground and background.

Examples include:

* `domain_count`;
* `max_domain_cov`;
* `ltr5_tfbm_burden`;
* `ltr3_tfbm_burden`;
* `total_ltr_tfbm_burden`;
* `rbp_burden`;
* `rbp_unique`;
* `n_terminal_domain_hits_protein_coding`.

The table reports foreground and background medians, delta median, means, p-value, direction and adjusted p-value.

A positive delta median means that the feature is higher in the foreground.

---

### Subfamily enrichment

The `subfamily_enrichment` table tests whether specific HERV subfamilies are overrepresented in the foreground.

This can reveal whether a HERV expression signature is dominated by a particular family or lineage.

---

### Domain type enrichment

The `domain_type_enrichment` table tests whether specific internal retroviral domain types are enriched in the foreground.

Examples include:

* Env;
* Gag;
* Protease;
* RT;
* RNaseH;
* Integrase;
* dUTPase.

This layer is useful for distinguishing generic HERV activation from activation of loci with specific retroviral protein-domain annotations.

---

### LTR TFBM burden and motif enrichment

HERVariumR separates 5′ LTR and 3′ LTR annotations when possible.

The `ltr_tfbm_burden` table compares the global number of predicted transcription-factor-binding motif hits in:

* 5′ LTRs;
* 3′ LTRs;
* all LTRs combined.

When detailed FIMO files are provided, HERVariumR also performs motif-level enrichment separately for:

* `ltr5_tfbm_motif_enrichment`;
* `ltr3_tfbm_motif_enrichment`;
* `any_ltr_tfbm_motif_enrichment`.

This allows users to ask whether a HERV list is enriched for specific predicted transcription-factor-binding motifs in 5′ LTRs, 3′ LTRs, or either LTR.

---

### RBP burden and motif enrichment

When an RBP motif file is provided, HERVariumR compares:

* total RBP motif burden;
* number of unique RBP motifs;
* individual RBP motif enrichment.

Relevant output tables include:

* `rbp_burden`;
* `rbp_motif_enrichment`;
* `herv_rbp_summary`;
* `rbp_hits`.

This layer helps explore whether HERV RNAs in a given list may differ in predicted post-transcriptional regulatory potential.

---

## Example biological application: cALD blood–brain barrier model

HERVariumR can be applied to HERV differential-expression results from disease-relevant transcriptomic datasets.

As an example, the package includes HERV differential-expression results from a reanalysis of the GSE108012 dataset, a public bulk RNA-seq dataset of iPSC-derived brain microvascular endothelial cells used as a blood–brain barrier model of childhood cerebral adrenoleukodystrophy.

In this use case, HERVariumR compares significantly upregulated HERVs in ccALD-derived iBMECs against all HERVs tested in the differential-expression analysis.

The goal is not to prove causality, but to show how HERVariumR converts a flat HERV differential-expression list into a multilayer biological profile.

For example, HERVariumR can reveal whether upregulated HERVs are enriched for:

* domain-containing loci;
* Env-containing loci;
* specific retroviral domain types;
* 5′ or 3′ LTR motif burden;
* interferon-related LTR motif features;
* specific transcription-factor-binding motifs;
* RBP motif burden;
* specific RBP motifs;
* particular HERV subfamilies;
* protein-coding or lncRNA transcript contexts.

---

## Recommended foreground/background design

For differential-expression analyses, the recommended setup is:

```text
foreground = significant HERVs of interest
background = all HERVs tested in the same analysis
```

Examples:

```text
foreground = significantly upregulated HERVs
background = all tested HERVs
```

or:

```text
foreground = significantly downregulated HERVs
background = all tested HERVs
```

HERVariumR removes foreground IDs from the background internally, so the comparison is made against tested non-foreground HERVs.

Avoid using all annotated HERVs in the genome as background unless all of them could realistically have been detected and tested in the experiment. For transcriptomic analyses, the best background is usually the set of HERVs that passed expression filtering and entered the differential-expression model.

---

## Domain coverage cutoff

Internal retroviral domains are filtered using `domain_coverage_cutoff`.

For example:

```r
domain_coverage_cutoff = 0.40
```

keeps only domain hits covering at least 40% of the HMM profile.

Using:

```r
domain_coverage_cutoff = 0
```

keeps all annotated domain hits.

A stricter cutoff is useful when the goal is to focus on more conserved domain annotations. A permissive cutoff can be useful for exploratory analyses or when the user wants to retain weak or partial domain evidence.

---

## Large external files

Some optional layers require large external files.

### LTR TFBM file

Used when:

```r
add_tfbm_details = TRUE
```

Required argument:

```r
fimo_file = "path/to/fimo_parsed_v4.tsv"
```

Expected content: FIMO-like transcription-factor-binding motif results mapped to LTR sequence names.

---

### RBP motif file

Used when:

```r
add_rbp_details = TRUE
```

Required argument:

```r
rbp_file = "path/to/RBP_fimo.tsv"
```

Expected content: FIMO-like RBP motif results mapped to HERV internal-region sequence names.

---

## Generated files

Depending on the selected options, HERVariumR writes output files to `output_dir`.

Common outputs include:

* input summaries;
* matched annotation tables;
* feature matrices;
* binary feature enrichment tables;
* numeric feature comparison tables;
* subfamily enrichment tables;
* domain type enrichment tables;
* LTR TFBM burden tables;
* motif enrichment tables;
* RBP enrichment tables;
* PNG plots;
* interactive HTML dashboards.

The dashboard can be generated with:

```r
generate_hervarium_dashboard(object)
```

where `object` is either a `HERVarium_annotation` or `HERVarium_comparison` object.

---

## Minimal example

```r
library(HERVariumR)

deg_file <- system.file(
  "sample_data",
  "HERV_DESeq2_DEG_all.tsv",
  package = "HERVariumR"
)

deg <- read.delim(deg_file, stringsAsFactors = FALSE)

foreground <- deg$transcript[
  deg$log2FoldChange > 0 &
    deg$padj < 0.05
]

background <- deg$transcript

cmp <- compare_herv_lists(
  foreground_ids = foreground,
  background_ids = background,
  foreground_name = "significant_up",
  background_name = "tested_non_up",
  output_dir = "HERVariumR_example",
  make_plots = TRUE
)

generate_hervarium_dashboard(cmp)
```

---

## Citation

If you use HERVariumR, please cite the associated manuscript when available.

For now, please cite the GitHub repository:

```text
Montserrat Ayuso T. HERVariumR: fast functional annotation and comparison of HERV lists.
GitHub: https://github.com/funcgen/HERVariumR
```

---

## Notes and limitations

HERVariumR is designed for functional interpretation and hypothesis generation.

Important points:

* HERV expression and HERV feature enrichment are correlative.
* Predicted TFBM and RBP motifs do not prove physical binding.
* Internal domain annotations indicate sequence similarity and conservation, not necessarily protein expression or function.
* Enrichment results depend strongly on the background universe.
* For transcriptomic studies, the most appropriate background is usually all HERVs tested in the same differential-expression analysis.
* Large FIMO files can be memory-intensive; `use_awk = TRUE` is recommended when available.
* Results should be interpreted together with experimental design, expression quantification strategy, mappability and biological context.

---

## License

Add license information here.

---

## Contact

For questions, suggestions or issues, please use the GitHub issue tracker:

```text
https://github.com/funcgen/HERVariumR/issues
```
