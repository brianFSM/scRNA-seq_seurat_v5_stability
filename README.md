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
| `scRNA_part1_QC.Rmd` | CellRanger output dirs | Loads matrices (optional SoupX), computes per-cell QC metrics, plots raw distributions with filter thresholds shaded | `_0_raw_seurat_object.RDS`, QC tables, PNGs |
| `scRNA_part2a_integration.Rmd` | part 1 RDS | Applies QC cutoffs, optional DoubletFinder, SCTransform → PCA → UMAP → RPCA integration → builds the `int_snn` graph → **subsample-stability resolution sweep** | integrated RDS (with graph + upstream-param stamp), stability/clustree PNGs, **writes a proposed resolution to config** |
| `scRNA_part2b_clustering.Rmd` | part 2a RDS | Staleness-checks the handoff, clusters at the config resolution, runs `FindAllMarkers`, generates per-cluster feature/dot plots | clustered RDS, marker CSVs, PNGs, **writes the derived resolution label + cluster count to config** |
| `scRNA_part3_summary.Rmd` | PNGs + CSVs from parts 1–2 | Stitches everything into one client-facing summary PDF. Computes nothing new. | summary PDF |
| `scRNA_part4_group_comparisons.Rmd` | part 2b RDS | **Template/stub.** Between-group DE (e.g. by sex, treatment). Edit before use — see [Part 4 is a template](#part-4-is-a-template). | group marker CSVs |

`functions_for_report.R` holds all helpers (path handling, 10x loading, metric
tables, plotting, the stability selector). Sourced by every report.

**The intended loop:** render part 1 → read the QC PDF → set your cutoffs in
`config.yaml` → render 2a + 2b (2a proposes a clustering resolution, 2b clusters
and discovers markers) → render part 3 for the summary. Re-clustering later is
cheap: edit one config value and re-run 2b alone (see [Running](#running)).

---

## Setup on Quest

**1. Clone**

RStudio on an analytics node: *File → New Project → Version Control → Git*, with
repository URL `https://github.com/brianFSM/scRNA-seq_seurat_V5.git`

Or from the command line:

```bash
git clone https://github.com/brianFSM/scRNA-seq_seurat_V5.git
cd scRNA-seq_seurat_V5
module purge all
module load R/4.4.0
```

**2. Restore the environment**

The `renv.lock` committed in this repo is the source of truth. From R, in the
project directory:

```r
renv::init(bare = TRUE)
renv::restore()
```

Packages symlink from the `renv` cache if you've built them before; otherwise
they compile fresh (slow the first time).

**3. PDF rendering (TinyTeX)** — the reports build to PDF via xelatex. First time
only:

```r
install.packages("tinytex")
tinytex::install_tinytex()   # sets up ~/.TinyTeX
tinytex::is_tinytex()        # should return TRUE
```

**4. Two packages worth confirming are present.** If `renv::restore()` didn't
bring them (GitHub-only or optional), install and re-snapshot:

```r
# REQUIRED by part 2a: subsample-stability resolution selection
renv::install("crazyhottommy/scclusteval")

# OPTIONAL: clustree override plot in part 2a (chunk auto-skips if absent)
renv::install("clustree")

# RECOMMENDED: dramatically speeds up FindAllMarkers in part 2b
# (ragg often won't build under R/4.4.0, so pin an older one first)
install.packages(
  "https://cran.r-project.org/src/contrib/Archive/ragg/ragg_1.4.0.tar.gz",
  repos = NULL, type = "source")
devtools::install_github("immunogenomics/presto")

renv::snapshot()
```

> **Note for QDSC:** if you want the exact shared lockfile the core standardizes
> on, copy it in *before* `renv::restore()`
> (`cp /projects/b1197/.../renv.lock .`). This is optional and only works if you
> have access to that project space — the repo's own `renv.lock` is sufficient
> for anyone else.

---

## Configure the run

```bash
cp config_template.yaml config.yaml
```

Edit `config.yaml` to point at your data. `config.yaml` is disposable and gets
*rewritten* by the pipeline (parts 2a/2b write chosen-resolution values back into
it, which flattens comments) — keep `config_template.yaml` pristine as your
reference.

### Project + data

```yaml
project:
  analyst-name: Your Name
  project-name: 'MyProject'
  project-description: "One-line description; appears in the reports."
  organism: mouse            # 'mouse' or 'human' — sets mito/ribo/hb/Y gene patterns
  run_cell_cycle: 'TRUE'
  samples:                   # a real YAML sequence; dashes/spaces in names are fine
    - Sample1
    - Sample2

data:
  parent_directory_path: /path/to/CellRanger   # one sub-directory per sample
  rds-file-path: /path/to/results              # all output objects/plots land here
```

**Input directory assumption:** for each entry in `samples`,
`parent_directory_path/<sample>/` must contain a CellRanger `outs/` (or the
matrices directly), i.e. a `raw_feature_bc_matrix/` (preferred; `filtered_...` is
used as a fallback) plus `metrics_summary.csv` for the sequencing-QC table. Raw
is required if you enable SoupX. Samples with a missing `metrics_summary.csv` are
warned about and simply omitted from the metrics table rather than crashing.

### QC cutoffs (part 2a) — set these *after* reading the part 1 PDF

```yaml
part2:
  mito_ceiling: 5           # max % mitochondrial
  RNA_count_floor: 500
  RNA_count_ceiling: 25000
  feature_count_floor: 300
  feature_count_ceiling: 6000
  rbc_ceiling: 100          # max % hemoglobin
  ribo_floor: 0             # min % ribosomal — see the warning below
```

### Clustering resolution (part 2a proposes, part 2b applies)

```yaml
part2:
  # resolution grid swept for stability
  cluster_res_min: 0.1
  cluster_res_max: 1.0
  cluster_res_step: 0.1
  # subsample-stability sweep controls
  stability_reps: 5             # subsamples per resolution (more = steadier, slower)
  stability_subsample_frac: 0.8
  stability_threshold: 0.75     # VESTIGIAL — see note
  # THE knob you edit to re-cluster:
  clustering_resolution_value: 0.4
  # derived by part 2b, do not hand-edit:
  clustering_resolution: None   # e.g. "int_snn_res.0.4"
  num_clusters: 0
```

- `clustering_resolution_value` is the numeric resolution — the one value you
  edit by hand to re-cluster. Part 2a writes its proposed value here.
- `clustering_resolution` and `num_clusters` are *derived outputs* written by
  part 2b (the metadata column name and the actual cluster count). Don't edit
  them; 2b overwrites them every run.
- `stability_threshold` is **no longer used** — the current selector derives its
  threshold from the data (the highest lower-confidence-bound across
  resolutions). The key is retained only for backward compatibility and can be
  ignored or removed.

### Toggles

```yaml
part1:
  run_SoupX: 'TRUE'              # ambient RNA correction; needs the raw matrix
  generate_metrics_tables: 'TRUE'
  part1_suffix: ''                # suffixes part-1 output filenames (alt runs)
part2:
  run_doubletFinder: 'TRUE'      # per-sample doublet detection before merging
  part2_suffix: ''                # suffixes part-2 output filenames (alt runs)
```

`run_SoupX` and `run_doubletFinder` are on by default; if you turn them off, say
so when reporting results. The `*_suffix` keys let you keep alternate runs
(different cutoffs) side by side by suffixing every output filename — parts 2b
and 3 stay in sync via a shared filename helper.

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
edit `clustering_resolution_value` in `config.yaml`, then re-run 2b alone:

```bash
sbatch run_templates.sh scRNA_part2b_clustering.Rmd
```

**Alternate config file** — select it with the `CONFIG_FILE` env var (default
`config.yaml`):

```bash
sbatch --export=ALL,CONFIG_FILE=config_test.yaml run_templates.sh scRNA_part1_QC.Rmd
```

Output PDFs are named after the Rmd (`scRNA_part1_QC.pdf`, etc.); objects,
tables, and plots go to `rds-file-path`. Edit the `#SBATCH --mail-user` line in
`run_templates.sh` to your address; the defaults (`--mem=80G`, `-t 48:00:00`,
account `b1042`, genomics partition) suit a large multi-sample run — trim for
small data.

> **Caveat:** the output PDF name comes from the *Rmd*, not the config, so a
> `CONFIG_FILE=config_test.yaml` run of part 1 overwrites `scRNA_part1_QC.pdf`.
> If you're running a real and a test config side by side, that's the collision
> to watch — the RDS/plot outputs are separated by `*_suffix`, but the PDF is not.

**On rendering:** hitting *Knit* in RStudio on an analytics node has been known
to hang forever; submitting via SLURM works. The reports run under a sequential
`future` plan for this reason.

### Running elsewhere

Nothing about the analysis is Quest-specific — only the module loads, the SLURM
header, and the reference GTF paths in part 1 are. To run locally, restore with
`renv`, then render directly:

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
  creation strips empty droplets only. The analytical feature floor is applied
  in part 2a, so the part-1 QC plots show honest raw distributions.
- **Normalization:** `SCTransform` (v2), applied to the merged object.
- **Integration:** `RPCAIntegration` on the SCT assay, 30 dims, with `k.weight`
  set dynamically for small samples to avoid the anchor-weighting error.
- **Clustering:** Louvain on an SNN graph (`int_snn`) built from the *integrated*
  embedding. The resolution is auto-selected by **subsample stability** rather
  than silhouette: part 2a repeatedly resamples 80% of cells, re-clusters, and
  measures how reproducibly each cluster reappears (Jaccard, via `scclusteval`);
  it then picks the finest resolution whose median per-cluster stability clears a
  **data-derived** threshold (a chooseR-style decision rule). Part 2b clusters at
  the chosen (or your overridden) resolution.
- **Markers:** `FindAllMarkers` on the SCT assay (after `PrepSCTFindMarkers`),
  Wilcoxon rank-sum, `only.pos = TRUE`, Bonferroni-corrected `p_val_adj`. These
  are cluster *identifiers*, not a formal differential-expression analysis.

**Two choices a reader should question — surfaced, not buried:**

1. **`ribo_floor` defaults to 0 (off), deliberately.** When set, it removes cells
   with *less* ribosomal content than the threshold — a non-standard QC step that
   can delete legitimate low-ribosomal cell types. Turn it on only with a reason.
2. **The auto-selected resolution is a starting point, not a verdict.** Stability
   selection resists under-clustering, but it can also suppress small, real
   populations (rare clusters are inherently less reproducible under
   resampling). Inspect `stability_vs_resolution.png` (and `clustree.png`, if
   generated) in the summary; if the pick doesn't match the biology you expect,
   set `clustering_resolution_value` by hand and re-run part 2b.

### Part 4 is a template

`scRNA_part4_group_comparisons.Rmd` is a starting point, not a finished analysis.
It demonstrates a between-group contrast (e.g. by sex) and **must be edited for
your design**. Two things to fix before using it for anything real:

1. Re-run `PrepSCTFindMarkers` after changing the grouping/idents.
2. Per-cell Wilcoxon across sample groups is **pseudoreplication** — it treats
   cells as independent replicates when the experimental unit is the sample. For
   a genuine group comparison, aggregate to **pseudobulk** per sample and test
   with DESeq2/edgeR.

---

## Notes / known limitations

- **Reproducibility is the point.** Use `renv`; don't install stray packages into
  the project library. `scclusteval` is a real dependency now (part 2a) — make
  sure it's in the lockfile.
- Alternate runs (different cutoffs) can coexist via `part1_suffix` /
  `part2_suffix`, which suffix every output filename.
- If something breaks, it's usually the config YAML, a missing module, or a
  package that didn't restore (`scclusteval`, `presto`, TinyTeX) — check those
  first.

## License

See `LICENSE`.
