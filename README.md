# scRNA-seq analysis pipeline (Seurat v5)

**Author:** Brian Wray

A config-driven, multi-report pipeline for a standard single-cell RNA-seq
analysis: QC → filtering → integration → clustering → marker genes → summary.
Each stage is an R Markdown report that renders to PDF, so every run leaves a
readable record of what was done and what the data looked like at each step.

Built for Northwestern's Quest HPC (SLURM + `renv`), but the analysis itself is
portable — see [Running elsewhere](#running-elsewhere).

---

## Pipeline at a glance

| Report | Reads | Does | Writes |
|---|---|---|---|
| `scRNA_part1_QC.Rmd` | CellRanger output dirs | Loads matrices (optional SoupX ambient-RNA correction), computes per-cell QC metrics, plots raw distributions | `_0_raw_seurat_object.RDS`, QC tables, raw QC PNGs |
| `scRNA_part2a_integration.Rmd` | part 1 RDS | Applies QC cutoffs, optional DoubletFinder, SCTransform → PCA → UMAP → RPCA integration → builds the `int_snn` graph → **subsample-stability resolution sweep** | integrated RDS (graph + upstream-param stamp), `proposed_resolution.txt`, filtering/attrition CSVs, stability + clustree PNGs |
| `scRNA_part2b_clustering.Rmd` | part 2a RDS | Staleness-checks the handoff, clusters at the resolved resolution, runs `FindAllMarkers`, draws per-cluster marker plots | clustered RDS, marker CSVs, cluster PNGs, `run_parameters.csv` |
| `scRNA_part3_summary.Rmd` | figures + tables from all three stages | Stitches everything into one client-facing summary PDF. Computes nothing new. | summary PDF |
| `scRNA_part4_group_comparisons.Rmd` | part 2b RDS | **Template/stub.** Between-group DE (e.g. by sex, treatment). Edit before use — see [Part 4 is a template](#part-4-is-a-template). | group marker CSVs |

`functions_for_report.R` holds all helpers — stage-directory resolution, 10x
loading, metric tables, plotting, the stability selector. Sourced by every
report.

**The intended loop:** render part 1 → read the QC PDF → set your cutoffs in
`config.yaml` → render 2a + 2b (2a proposes a clustering resolution, 2b clusters
and finds markers) → render part 3 for the summary. Re-clustering later is cheap:
set one config value and re-run 2b alone (see [Running](#running)).

---

## Output layout

Each stage writes into a directory nested inside the stage it depends on, so a
directory's path *is* its provenance — a part 2b result can't be paired with the
wrong part 2a, because it lives inside it.

```
<rds-file-path>/
  part1[_<part1_suffix>]/
    _0_raw_seurat_object.RDS
    part1_tables_for_report.RData
    ggplot/                          raw QC distributions
    part2a[_<part2a_suffix>]/
      _1_seurat5_obj_merged.RDS
      _2_seurat5_obj_integrated.RDS
      part2_tables_for_report.RData
      filtering_stats_table.csv
      qc_attrition_table.csv
      proposed_resolution.txt
      ggplot/                        before/after QC, merged + integrated UMAPs,
                                     stability sweep, clustree
      res_<resolution>/
        _3_seurat5_obj_clustered.RDS
        cells_perCluster_perSample.csv
        run_parameters.csv
        marker_genes/                per-cluster and combined marker CSVs
        ggplot/                      cluster UMAP, PCA, per-cluster marker plots
```

Three consequences worth knowing:

- **Part 1's object is shared.** It doesn't depend on QC cutoffs, so re-running
  2a with different cutoffs never re-runs the expensive load.
- **Resolutions don't collide.** Clustering at 0.4 and then 0.6 produces
  `res_0.4/` and `res_0.6/` side by side, each with its own markers and figures.
- **`res_<resolution>/` is the client deliverable.** The part 3 PDF and
  `marker_genes/` both live there, so it can be zipped and sent as-is.

Suffixes name directories, not files. `part1_suffix: soupx_off` gives
`part1_soupx_off/`; leaving it `''` gives `part1/`. A separator is inserted for
you, so `alt` and `_alt` behave identically.

---

## Setup on Quest

**1. Clone**

RStudio on an analytics node: *File → New Project → Version Control → Git*, with
repository URL
`https://github.com/brianFSM/scRNA-seq_seurat_v5_stability.git`

Or from the command line:

```bash
git clone https://github.com/brianFSM/scRNA-seq_seurat_v5_stability.git
cd scRNA-seq_seurat_v5_stability
module purge all
module load R/4.5.1
```

**2. Restore the environment**

The `renv.lock` committed in this repo is the source of truth, and includes the
GitHub-only packages (`scclusteval`, `presto`, `DoubletFinder`). From R, in the
project directory:

```r
renv::init(bare = TRUE)
renv::restore()
```

Packages symlink from the `renv` cache if you've built them before; otherwise
they compile fresh (slow the first time).

If a GitHub package fails to restore — most often after an R version change —
install it by hand and re-snapshot:

```r
renv::install("remotes")
remotes::install_github("crazyhottommy/scclusteval")     # REQUIRED by part 2a
remotes::install_github("chris-mcginnis-ucsf/DoubletFinder", force = TRUE)
remotes::install_github("immunogenomics/presto")
renv::snapshot()
```

`presto` depends on `ragg`, which sometimes won't build on Quest. If it fails,
pin an older source version first:

```r
install.packages(
  "https://cran.r-project.org/src/contrib/Archive/ragg/ragg_1.4.0.tar.gz",
  repos = NULL, type = "source")
```

**3. PDF rendering (TinyTeX)** — the reports build to PDF via xelatex. First time
only:

```r
install.packages("tinytex")
tinytex::install_tinytex()   # sets up ~/.TinyTeX
tinytex::is_tinytex()        # should return TRUE

tinytex::reinstall_tinytex() # if you move to a new R version
```

**4. Optional packages.** Everything below is optional; the pipeline degrades
gracefully without it.

| Package | Effect if missing |
|---|---|
| `clustree` + `ggraph` | Part 2a's clustree chunk skips silently |
| `presto` | `FindAllMarkers` in part 2b still runs, but much more slowly |

---

## Configure the run

```bash
cp config_template.yaml config.yaml
```

`config.yaml` is **input only** — nothing in the pipeline writes to it, so your
comments and manual edits survive. `config_template.yaml` is the annotated
reference; read it alongside this section.

`config.yaml` is gitignored so your real paths stay out of version control.

### Project + data

```yaml
project:
  analyst-name: Your Name
  project-name: 'Test'
  project-description: "ToDo"
  organism: mouse            # 'mouse' or 'human'
  samples:
    - Sample1
    - Sample2

data:
  parent_directory_path: path/to/CellRanger
  rds-file-path: ./results   # output root; stages nest beneath it
```

`samples` are directory names under `parent_directory_path`, each holding
CellRanger output. `organism` selects the mito/ribo/hemoglobin patterns and the
Y-chromosome gene list; anything other than `human` or `mouse` stops with an
explanatory error.

### Part 1

```yaml
analysis:
  node_type: analytics

  part1:
    part1_rds_save_filename: _0_raw_seurat_object.RDS
    part1_report_tables_filename: part1_tables_for_report.RData
    part1_suffix: ''           # names part1_<suffix>/
    generate_metrics_tables: 'TRUE'
    run_SoupX: 'TRUE'          # ambient RNA correction; needs the raw matrix
    ggplot_dir: ggplot
```

Use `part1_suffix` only when **part 1 itself** varies — SoupX on versus off, say.
Different QC cutoffs are a part 2a concern and get `part2a_suffix` instead.

### QC cutoffs (part 2a) — set these *after* reading the part 1 PDF

```yaml
  part2:
    part2a_suffix: ''          # names part2a_<suffix>/ inside the part 1 dir

    mito_ceiling: 5            # max % mitochondrial
    RNA_count_floor: 500
    RNA_count_ceiling: 25000
    feature_count_floor: 300
    feature_count_ceiling: 6000
    rbc_ceiling: 100           # max % hemoglobin
    ribo_floor: 0              # see the caveat in Methods summary
```

Change a cutoff and you should change `part2a_suffix` too, so the new run lands
beside the old one instead of overwriting it. Part 2a stamps its parameters into
the integrated object and part 2b refuses a stale handoff, so a forgotten suffix
surfaces as an error rather than a silent mismatch.

### Clustering resolution (part 2a proposes, part 2b applies)

```yaml
    cluster_res_min: 0.1
    cluster_res_max: 1.0
    cluster_res_step: 0.1
    stability_reps: 5
    stability_subsample_frac: 0.8
    clustering_resolution_value: auto
```

Part 2a sweeps the grid, scores each resolution for subsample stability, and
writes its pick to `proposed_resolution.txt` in its own directory. With
`clustering_resolution_value: auto`, part 2b uses that. Set a number instead to
override; part 2a never overwrites your value.

`stability_reps` is the main cost knob — the sweep runs that many subsampled
re-clusterings per resolution.

### Output filenames and toggles

```yaml
    merged_seurat_obj_file: _1_seurat5_obj_merged.RDS
    integrated_seurat_obj_file: _2_seurat5_obj_integrated.RDS
    clustered_seurat_obj_file: _3_seurat5_obj_clustered.RDS
    part2_report_tables_filename: part2_tables_for_report.RData
    marker_gene_dir: marker_genes

    run_doubletFinder: 'TRUE'
    doublet_rate_per_1k: 0.008   # expected doublets per 1,000 cells recovered
    doublet_rate_cap: 0.20
    marker_plot_types:           # any of feat, dot, vln
      - feat
      - dot
```

Filenames are relative to whichever stage directory owns them. The doublet rate
is applied per sample and scales with that sample's cell count, so a small
library isn't charged a large library's rate.

---

## Running

Submit each report as a SLURM job. `run_templates.sh` renders **any number of
Rmds in order**, in one job:

```bash
sbatch run_templates.sh scRNA_part1_QC.Rmd
# read scRNA_part1_QC.pdf, set your cutoffs in config.yaml, then:
sbatch run_templates.sh scRNA_part2a_integration.Rmd scRNA_part2b_clustering.Rmd
sbatch run_templates.sh scRNA_part3_summary.Rmd
```

Rendering 2a and 2b in one job is deliberate: they run sequentially in the same
session, and `render()` stops on error, so if 2a fails, 2b never runs (fail-fast
ordering without a SLURM dependency). Because they run sequentially, peak memory
is whichever stage is larger, not the sum.

**Re-cluster at a different resolution** (fast — skips integration entirely):
set `clustering_resolution_value` to a number, then re-run 2b alone. The results
land in a new `res_<resolution>/` beside the existing one.

```bash
sbatch run_templates.sh scRNA_part2b_clustering.Rmd scRNA_part3_summary.Rmd
```

**Alternate config file** — select it with the `CONFIG_FILE` env var (default
`config.yaml`):

```bash
sbatch --export=ALL,CONFIG_FILE=config_test.yaml run_templates.sh scRNA_part1_QC.Rmd
```

Edit the `#SBATCH --mail-user` line in `run_templates.sh` to your address; the
defaults (`--mem=80G`, `-t 48:00:00`, account `b1042`, genomics partition) suit a
large multi-sample run — trim for small data.

> **Caveat:** output PDFs are named after the *Rmd* and land in the project
> directory, not in the stage directory with everything else. So a second run —
> a different `part2a_suffix`, or `CONFIG_FILE=config_test.yaml` — overwrites the
> previous PDF even though its objects, tables, and figures are kept separate.
> Copy the PDF into its stage directory before re-running, or rename it.

**On rendering:** hitting *Knit* in RStudio on an analytics node has been known
to hang forever; submitting via SLURM works. The reports run under a sequential
`future` plan for this reason.

### Running elsewhere

Nothing about the analysis is Quest-specific — only the module loads and the
SLURM header are. To run locally, restore with `renv`, then render directly:

```r
rmarkdown::render("scRNA_part1_QC.Rmd",
                  params = list(config.args = "config.yaml"),
                  envir  = new.env(parent = globalenv()))
```

The `CONFIG_FILE` mechanism isn't needed here — pass the config path in `params`
directly, or render the 2a/2b pair with two `render()` calls (fresh env each),
which is exactly what the SLURM path does.

---

## Methods summary

The defaults, and why they're defensible for a basic pipeline:

- **Loading:** the raw 10x matrix is preferred; `min.features = 100` at object
  creation strips empty droplets only. The analytical feature floor is applied in
  part 2a, so the part-1 QC plots show honest raw distributions.
- **Ambient RNA:** optional SoupX. `autoEstCont()` estimates each sample's
  contamination fraction independently; part 1 prints the diagnostic plot per
  sample with its estimate.
- **Doublets:** optional DoubletFinder per sample, before merging. The expected
  rate scales with cells recovered (`doublet_rate_per_1k` × cells ÷ 1,000, capped
  at `doublet_rate_cap`), then is reduced by the estimated homotypic proportion.
- **Normalization:** `SCTransform` (v2), applied to the merged object.
- **Integration:** `RPCAIntegration` on the SCT assay, 30 dims, with `k.weight`
  lowered automatically for small samples to avoid the anchor-weighting error.
- **Clustering:** Louvain on the SNN graph (`int_snn`) built from the *integrated*
  embedding. The resolution is auto-selected by **subsample stability** rather
  than silhouette: part 2a resamples 80% of cells several times, re-clusters each
  subsample, and measures how reproducibly each cluster reappears (Jaccard, via
  `scclusteval`). Every resolution is scored on the same subsamples, so the
  comparison between them is paired. It then picks the finest resolution whose
  median per-cluster stability clears a **data-derived** threshold (a
  chooseR-style decision rule). Part 2b clusters at the chosen — or your
  overridden — resolution.
- **Markers:** `FindAllMarkers` on the SCT assay (after `PrepSCTFindMarkers`),
  Wilcoxon rank-sum, `only.pos = TRUE`, Bonferroni-corrected `p_val_adj`. These
  are cluster *identifiers*, not a formal differential-expression analysis.

**Two choices a reader should question — surfaced, not buried:**

1. **`ribo_floor` defaults to 0 (off), deliberately.** When set, it removes cells
   with *less* ribosomal content than the threshold — a non-standard QC step that
   can delete legitimate low-ribosomal cell types. Turn it on only with a reason.
2. **The auto-selected resolution is a starting point, not a verdict.** Stability
   selection resists under-clustering, but it can also suppress small, real
   populations (rare clusters are inherently less reproducible under resampling).
   Inspect `stability_vs_resolution.png` — and `clustree.png`, if generated — in
   the summary; if the pick doesn't match the biology you expect, set
   `clustering_resolution_value` by hand and re-run part 2b.

### Part 4 is a template

`scRNA_part4_group_comparisons.Rmd` is a starting point, not a finished analysis.
It demonstrates a between-group contrast (e.g. by sex) and **must be edited for
your design**. Two things to fix before using it for anything real:

1. The group assignment is a placeholder that infers sex from a sample-name
   pattern. Replace it with your actual design, and check that every sample
   landed in the group you expect.
2. Per-cell Wilcoxon across sample groups is **pseudoreplication** — it treats
   cells as independent replicates when the experimental unit is the sample. For
   a genuine group comparison, aggregate to **pseudobulk** per sample and test
   with DESeq2/edgeR.

---

## Notes / known limitations

- **Reproducibility is the point.** Use `renv`; don't install stray packages into
  the project library. If you add a dependency, `renv::snapshot()` and commit the
  lockfile.
- **Freshness is checked one step back, not all the way up.** Part 2b refuses an
  integrated object built under different QC or integration parameters. Nothing
  yet detects a *part 1* re-run invalidating the 2a and 2b results beneath it —
  if you re-run part 1, treat everything under it as stale.
- **Suffix discipline is manual.** Nothing stops you from changing a cutoff
  without changing `part2a_suffix`; the staleness check will catch the resulting
  mismatch at part 2b, but the directory name will no longer describe its
  contents. Every report prints its resolved paths and parameters, and part 3
  includes a run-parameters table, so a mistake is traceable after the fact.
- If something breaks, it's usually the config YAML, a missing module, or a
  package that didn't restore (`scclusteval`, `presto`, TinyTeX) — check those
  first.

## License
