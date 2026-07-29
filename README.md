# HERVariumR

<p align="center">
  <img src="man/figures/HERVariumR_logo.png" width="350">
</p>

**HERVariumR** is an R package for fast functional annotation, profiling, and comparison of human endogenous retrovirus (HERV) locus lists.

The package is designed for researchers who quantify HERV expression at the locus level and obtain lists of HERVs of interest, such as differentially expressed HERVs, condition-specific HERVs, interferon-responsive HERVs, or HERVs associated with a given cell type or disease state.

The main philosophy is simple: **start from a list of HERV IDs and obtain a compact biological interpretation**.

---

## What HERVariumR integrates

HERVariumR converts HERV lists into multilayer summaries by integrating:

* internal retroviral protein-domain annotations;
* Gag, Pol, Env and accessory-domain classification;
* HERV subfamily information;
* 5′ and 3′ LTR annotation;
* LTR transcription-factor-binding motif burden;
* interferon-related STAT1, STAT1::STAT2 and IRF motif summaries;
* compact ENCODE/UCSC 44-organ regulatory evidence for internal HERV regions and LTRs;
* compact cCRE overlap summaries;
* compact LTR TF ChIP-seq peak-cluster overlap evidence;
* transcript-context annotation, including intronic, exonic, protein-coding and lncRNA overlap;
* terminal-exon domain annotations;
* optional detailed LTR TFBM motif enrichment from large FIMO files;
* optional RNA-binding protein motif annotation from large FIMO-like files;
* compact standard outputs, optional extended static plots, and interactive HTML dashboards.

Two regulatory concepts are intentionally kept separate:

```text
TFBM      = predicted TF motif in the LTR sequence
TF ChIP   = experimental external TF ChIP-seq peak-cluster overlap with the LTR interval
```

Both can be useful, but they answer different questions.

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

HERVariumR supports two common workflows.

### 1. Profile a single HERV list

Use this when you have one list of HERVs and want to know what kind of loci they are.

Examples:

* HERVs upregulated in one condition;
* HERVs expressed in one cell type;
* HERVs overlapping accessible chromatin;
* HERVs associated with a phenotype;
* candidate HERV loci selected from another analysis.

Main function:

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

Main function:

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
hervarium_file("internal_encode_44organs_compact.tsv.gz")
hervarium_file("ltr_encode_44organs_compact.tsv.gz")
hervarium_file("ltr_tfchip_compact.tsv.gz")
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

deg <- read.delim(
  deg_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)
```

Define a list of significantly upregulated HERVs:

```r
test_ids <- unique(deg$transcript[
  !is.na(deg$log2FoldChange) &
    !is.na(deg$padj) &
    deg$log2FoldChange > 0 &
    deg$padj < 0.05
])
```

Run a standard HERVariumR profile:

```r
res_cALD_BBB <- profile_hervs(
  herv_ids = test_ids,
  output_dir = "hervarium_test_results/profile",
  output_level = "standard"
)
```

By default, `profile_hervs()` adds the bundled compact external-evidence layers:

```r
add_external_evidence = TRUE
add_tfchip = TRUE
```

The standard workflow automatically writes a compact HERV-level table, an RDS object, an input summary, and a concise interactive HTML dashboard. Static plots are disabled by default in standard mode.

`generate_hervarium_dashboard()` remains available when a dashboard needs to be regenerated, for example with a different `top_n` value:

```r
generate_hervarium_dashboard(
  res_cALD_BBB,
  top_n = 25
)
```

`profile_hervs()` returns a `HERVarium_annotation` object containing the matched annotation table, compact feature summaries, missing IDs, summary statistics, output metadata, and any optional regulatory layers.

## Compare a foreground HERV list against a background

A typical transcriptomic use case is to compare significantly upregulated HERVs against all HERVs tested in the same differential-expression analysis.

```r
foreground <- unique(deg$transcript[
  !is.na(deg$log2FoldChange) &
    !is.na(deg$padj) &
    deg$log2FoldChange > 0 &
    deg$padj < 0.05
])

background <- unique(deg$transcript)
```

Run the standard comparison:

```r
cmp_cALD_BBB <- compare_herv_lists(
  foreground_ids = foreground,
  background_ids = background,
  foreground_name = "cALD_BBB_sig_up",
  background_name = "tested_non_up",
  output_dir = "hervarium_test_results/comparison",
  output_level = "standard",
  plot_top_n = 25
)
```

The standard workflow automatically writes the compact comparison table, the unified key-results table, an RDS object, an input summary, and a concise interactive HTML dashboard.

Although the supplied background contains all tested HERVs, HERVariumR removes foreground IDs from the background internally. Therefore, the actual comparison is:

```text
significantly upregulated HERVs
versus
tested but not significantly upregulated HERVs
```

This is the recommended design for enrichment-style analyses of differential-expression results.

## Full profiling workflow with detailed TFBM and RBP layers

Detailed motif-level files are optional and are not bundled with the package. When these files are available, HERVariumR can add LTR TFBM and RBP motif summaries.

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

  # Preserve detailed tables and static plots
  output_level = "extended",
  make_plots = TRUE
)
```

The `use_awk = TRUE` option is recommended for large FIMO files because it filters the file before importing the selected records into R.

Detailed TFBM and RBP layers remain predictions derived from motif matches. They should not be interpreted as validated TF or RBP binding events.

## Output objects

### `HERVarium_annotation`

Returned by:

```r
profile_hervs()
annotate_hervs()
```

Core components include:

* `summary`: matched HERV-level annotation table;
* `features`: compact user-facing HERV-level feature table;
* `missing_ids`: input IDs not found in the annotation;
* `stats`: summary statistics for the input list;
* `domain_coverage_cutoff`: domain coverage threshold used;
* `ltr_tfbm_hits`: optional detailed LTR TFBM hits;
* `herv_tfbm_summary`: optional HERV-level TFBM summaries;
* `top_tfbm_motifs`: optional motif-frequency table;
* `rbp_hits`: optional detailed RBP motif hits;
* `herv_rbp_summary`: optional HERV-level RBP summaries.

Objects returned by `profile_hervs()` additionally contain:

* `output_level`: `"standard"` or `"extended"`;
* `output_files`: paths to the principal generated files and optional output directories;
* `output_dir`: output directory used by the run.

Internal wide feature matrices may be used transiently while annotations, plots, and statistical summaries are computed, but the principal returned and written feature tables are compact.

---

### `HERVarium_comparison`

Returned by:

```r
compare_herv_lists()
```

Common components include:

* `foreground_name`;
* `background_name`;
* `foreground_ids`;
* `background_ids`;
* `annotation`: the annotated HERV set used for the comparison;
* `features`: compact user-facing feature table with `comparison_group`;
* `input_summary`;
* `key_results`: unified main comparison-results table;
* `binary_features`;
* `numeric_features`;
* `external_evidence_binary_features`;
* `external_evidence_numeric_features`;
* `subfamily_enrichment`;
* `domain_hits`;
* `domain_type_enrichment`;
* `ifn_ltr_features`;
* `ltr_tfbm_burden`;
* `tfchip_burden`;
* `ltr5_tfchip_tf_enrichment`;
* `ltr3_tfchip_tf_enrichment`;
* `any_ltr_tfchip_tf_enrichment`;
* `rbp_burden`;
* `ltr5_tfbm_motif_enrichment`;
* `ltr3_tfbm_motif_enrichment`;
* `any_ltr_tfbm_motif_enrichment`;
* `rbp_motif_enrichment`;
* `domain_coverage_cutoff`;
* `output_level`;
* `output_files`;
* `output_dir`.

## Default compact feature tables

HERVariumR writes compact feature tables by default. Wider intermediate feature matrices may be created during annotation, plotting, and comparison, but they are pruned before the principal user-facing object and TSV outputs are returned.

The compact feature tables retain the main biological information: HERV identity, domain preservation, 5′/3′ LTR identity, LTR-specific TFBM burden, RBP burden, LTR-specific STAT1/STAT1::STAT2 motif flags, terminal-exon domain context, transcript context, compact internal ENCODE/UCSC evidence, LTR-specific ENCODE/UCSC plus cCRE evidence, and LTR-specific TF ChIP-seq overlap summaries.

For external evidence, the compact table intentionally avoids broad `any_ltr_*`, `both_ltr_*`, and `total_ltr_*` summaries when separate 5′ and 3′ LTR columns are more interpretable. Internal-region cCRE columns are also omitted from the compact table; cCRE is exposed for associated LTRs through columns such as:

```text
ltr5_encode_ccre_overlap
ltr3_encode_ccre_overlap
ltr5_encode_ccre_n_overlaps
ltr3_encode_ccre_n_overlaps
```

Layer identities are shown with detected-layer columns, for example:

```text
internal_encode_detected_layers = DNase;H3K27ac;transcription
ltr5_encode_detected_layers     = DNase;H3K4me3
ltr3_encode_detected_layers     = .
```

Detected transcription is also exposed by strand direction when available:

```text
internal_encode_same_strand_transcription_detected
internal_encode_opposite_strand_transcription_detected
internal_encode_bidirectional_transcription_detected
ltr5_encode_same_strand_transcription_detected
ltr5_encode_opposite_strand_transcription_detected
ltr5_encode_bidirectional_transcription_detected
ltr3_encode_same_strand_transcription_detected
ltr3_encode_opposite_strand_transcription_detected
ltr3_encode_bidirectional_transcription_detected
```

When the bundled TF ChIP-seq layer is enabled, the compact table retains the canonical, non-truncated LTR-specific TF lists:

```text
ltr5_tfchip_tf_list
ltr3_tfchip_tf_list
```

It also retains LTR-specific burden fields such as:

```text
has_ltr5_tfchip_overlap
has_ltr3_tfchip_overlap
ltr5_tfchip_n_peak_clusters
ltr3_tfchip_n_peak_clusters
ltr5_tfchip_n_tfs
ltr3_tfchip_n_tfs
```

Detailed TFBM TF-name columns are added only when `add_tfbm_details = TRUE` and a detailed FIMO file is supplied.

## Interpreting the main comparison tables

### Unified key-results table

`03_comparison_key_results.tsv` is the principal user-facing comparison-results table in standard mode. It combines key binary, numeric, subfamily, domain, external-regulatory, TF ChIP-seq, and optional motif-enrichment results using a shared schema:

```text
rank
layer
feature
feature_label
test_type
effect_type
foreground_value
background_value
effect_size
p_value
padj
direction
significant
```

The complete component-specific result tables remain available inside the returned `HERVarium_comparison` object. In extended mode, they are also written under `details/`.

### Binary feature enrichment

The `binary_features` component tests whether binary annotations are enriched in the foreground.

Examples include:

* `has_domain`;
* `has_gag`;
* `has_pol`;
* `has_env`;
* `has_ltr5`;
* `has_ltr3`;
* `has_both_ltrs`;
* `has_ltr5_stat1_motif`;
* `has_ltr3_stat1_motif`;
* `has_ltr5_stat1stat2_motif`;
* `has_ltr3_stat1stat2_motif`;
* `has_terminal_exon_domain_protein_coding`;
* `overlaps_exon`;
* `overlaps_intron`;
* `overlaps_lncRNA`;
* `overlaps_protein_coding`;
* `has_internal_encode_evidence`;
* `ltr5_encode_evidence_including_ccre`;
* `ltr3_encode_evidence_including_ccre`;
* `ltr5_encode_ccre_overlap`;
* `ltr3_encode_ccre_overlap`;
* `has_ltr5_tfchip_overlap`;
* `has_ltr3_tfchip_overlap`.

The table reports foreground and background counts and percentages, odds ratio, p-value, adjusted p-value, and direction. Positive enrichment means that the feature is more common in the foreground.

### Numeric feature shifts

The `numeric_features` component compares continuous or count-like annotations between foreground and background.

Examples include:

* `domain_count`;
* `max_domain_cov`;
* `ltr5_tfbm_burden`;
* `ltr3_tfbm_burden`;
* `rbp_burden`;
* `rbp_unique`;
* `max_terminal_domain_coverage_protein_coding`;
* `internal_encode_detected_layer_count`;
* `ltr5_encode_detected_layer_count`;
* `ltr3_encode_detected_layer_count`;
* `ltr5_encode_ccre_n_overlaps`;
* `ltr3_encode_ccre_n_overlaps`;
* `ltr5_tfchip_n_peak_clusters`;
* `ltr3_tfchip_n_peak_clusters`;
* `ltr5_tfchip_n_tfs`;
* `ltr3_tfchip_n_tfs`.

The table reports foreground and background medians and means, delta median, p-value, adjusted p-value, and direction. A positive delta median means that the feature is higher in the foreground.

### External evidence tables

The external-evidence components collect ENCODE/UCSC regulatory evidence, cCRE overlaps, and selected LTR TF ChIP-seq features. In extended mode they are written as:

```text
details/16_comparison_external_evidence_binary_features.tsv
details/17_comparison_external_evidence_numeric_features.tsv
```

These tables separate external regulatory context from core HERV identity, internal-domain, and transcript-context annotations.

### LTR TF ChIP-seq burden and TF enrichment

The compact LTR TF ChIP-seq layer asks whether associated HERV LTR intervals overlap TF ChIP-seq peak clusters.

In extended mode, the dedicated tables are written as:

```text
details/18_comparison_tfchip_burden.tsv
details/19_comparison_ltr5_tfchip_tf_enrichment.tsv
details/20_comparison_ltr3_tfchip_tf_enrichment.tsv
details/21_comparison_any_ltr_tfchip_tf_enrichment.tsv
```

The burden table compares peak-cluster counts, numbers of TFs, total overlap base pairs, and related overlap summaries. The per-TF tables test whether foreground HERVs are enriched for positional overlap with ChIP-seq evidence for particular TFs.

TF ChIP-seq overlap is experimental positional context. It does not prove autonomous TF recruitment by the HERV sequence.

### Subfamily enrichment

The `subfamily_enrichment` component tests whether specific HERV subfamilies are overrepresented in the foreground.

### Domain type enrichment

The `domain_type_enrichment` component tests whether specific internal retroviral domain types are enriched in the foreground, including Env, Gag, Protease, RT, RNaseH, Integrase, and dUTPase annotations.

Internal-domain annotation indicates sequence similarity and preservation. It does not establish protein expression or function.

### LTR TFBM burden and motif enrichment

HERVariumR keeps 5′ and 3′ LTR annotations separate when possible.

The `ltr_tfbm_burden` component compares the global number of predicted transcription-factor-binding motif hits in 5′ LTRs, 3′ LTRs, and both positions combined for the dedicated burden analysis.

When a detailed FIMO file is supplied, motif-level enrichment is retained separately in:

* `ltr5_tfbm_motif_enrichment`;
* `ltr3_tfbm_motif_enrichment`;
* `any_ltr_tfbm_motif_enrichment`.

TFBM results are sequence-based predictions, not validated binding events.

### RBP burden and motif enrichment

When an RBP motif file is supplied, HERVariumR compares total RBP motif burden, the number of unique RBP motifs, and individual RBP motif enrichment.

Relevant object components include:

* `rbp_burden`;
* `rbp_motif_enrichment`;
* `herv_rbp_summary`;
* `rbp_hits`.

RBP motif matches represent predicted sequence features, not validated RBP binding events.

## Interpreting external regulatory evidence

ENCODE/UCSC signal, cCRE overlap and TF ChIP-seq overlap provide **external regulatory context** overlapping HERV intervals.

They should not be interpreted as direct proof that:

* the HERV autonomously drives regulatory activity;
* the HERV is autonomously expressed;
* a nearby gene is regulated by the HERV;
* a TF physically binds because a motif is present.

Recommended interpretation:

```text
TFBM motif present       = sequence-based regulatory potential
TF ChIP-seq overlap      = experimental TF ChIP-seq peak-cluster context
cCRE overlap             = external candidate regulatory-element overlap
44-organ signal          = organ-level external epigenomic/transcriptional context
```

For ENCODE/UCSC signal thresholds, `detected`, `robust` and `high` should be read hierarchically. The compact default feature table reports detected-layer identities, while robust/high support is still used internally for selected summary features such as robust multilayer evidence:

```text
high ⊂ robust ⊂ detected
```

So a `*_detected` column means **at least detected signal**, not “detected-only excluding robust or high”.

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
* external ENCODE/UCSC regulatory evidence;
* LTR-specific cCRE overlap;
* LTR TF ChIP-seq peak-cluster overlap;
* specific TF ChIP-seq TF enrichments;
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

Both main workflows support:

```r
output_level = c("standard", "extended")
```

The default is `"standard"`.

### Standard single-list output

```text
herv_features_compact.tsv
input_summary.tsv
hervarium_annotation_object.rds
hervarium_annotation_dashboard.html
missing_herv_ids.txt          # only when needed
```

`herv_features_compact.tsv` is the principal per-locus output table.

### Standard comparison output

```text
01_input_summary.tsv
02_features_with_comparison_group.tsv
03_comparison_key_results.tsv
hervarium_comparison_object.rds
hervarium_comparison_dashboard.html
missing_herv_ids.txt          # only when needed
```

`03_comparison_key_results.tsv` is the principal statistical-results table.

### Extended output

Using:

```r
output_level = "extended"
```

retains all standard outputs and additionally writes:

```text
details/    # detailed annotation and comparison tables
plots/      # static plots when enabled
```

When `make_plots = NULL`, static plots default to disabled in standard mode and enabled in extended mode. They can be requested explicitly with `make_plots = TRUE`.

The interactive dashboard is generated automatically unless:

```r
make_dashboard = FALSE
```

A dashboard can also be regenerated from a returned object:

```r
generate_hervarium_dashboard(
  object,
  top_n = 15
)
```

Temporary R Markdown source files and rendering objects are created outside the user output directory and removed automatically.

## Minimal example

```r
library(HERVariumR)

deg_file <- system.file(
  "sample_data",
  "HERV_DESeq2_DEG_all.tsv",
  package = "HERVariumR"
)

deg <- read.delim(
  deg_file,
  stringsAsFactors = FALSE,
  check.names = FALSE
)

foreground <- unique(deg$transcript[
  !is.na(deg$log2FoldChange) &
    !is.na(deg$padj) &
    deg$log2FoldChange > 0 &
    deg$padj < 0.05
])

background <- unique(deg$transcript)

cmp <- compare_herv_lists(
  foreground_ids = foreground,
  background_ids = background,
  foreground_name = "significant_up",
  background_name = "tested_non_up",
  output_dir = "HERVariumR_example",
  output_level = "standard"
)
```

The compact comparison table, unified key-results table, RDS object, input summary, and interactive dashboard are generated automatically.

## Citation

If you use HERVariumR, please cite the associated manuscript when available.

For now, please cite the GitHub repository:

```text
Montserrat-Ayuso T. HERVariumR: fast functional annotation and comparison of HERV lists.
GitHub: https://github.com/funcgen/HERVariumR
```

---

## Notes and limitations

HERVariumR is designed for functional interpretation and hypothesis generation.

Important points:

* HERV expression and HERV feature enrichment are correlative.
* Predicted TFBM and RBP motifs do not prove physical binding.
* TF ChIP-seq overlap is external evidence of peak-cluster overlap, not proof that the HERV sequence drives TF recruitment.
* cCRE overlap is positional regulatory context, not direct target-gene evidence.
* Internal domain annotations indicate sequence similarity and conservation, not necessarily protein expression or function.
* Enrichment results depend strongly on the background universe.
* For transcriptomic studies, the most appropriate background is usually all HERVs tested in the same differential-expression analysis.
* Large FIMO files can be memory-intensive; `use_awk = TRUE` is recommended when available.
* Results should be interpreted together with experimental design, expression quantification strategy, mappability and biological context.

---

## License

MIT License
Copyright (c) 2026 Tomàs Montserrat-Ayuso

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

---

## Contact

For questions, suggestions or issues, please use the GitHub issue tracker:

```text
https://github.com/funcgen/HERVariumR/issues
```
