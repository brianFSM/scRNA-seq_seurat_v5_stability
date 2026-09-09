# functions_for_report.R


# Helper function that formats filenames. Allows for adding
# suffix to config for alternate versions of analysis (e.g. 
# different clustering resolutions, filters, etc)
suffix_path <- function(base_dir, filename, suffix = "") { 

	root <- tools::file_path_sans_ext(filename) 
	ext  <- tools::file_ext(filename) 
	fn   <- if (nzchar(suffix)) paste0(root, suffix, if (nzchar(ext)) paste0(".", ext) else "") else filename 
	file.path(base_dir, fn)
}

# Helper function to make the sample names pretty and workable.
#
# `samples` in the config is now a proper YAML sequence, so this normally
# receives a character vector and just needs light cleanup. The legacy path
# (a single dash/quote-delimited string) is still handled for back-compat, so
# old config files don't break.
sanitize_sample_names <- function(formatted.samples){

  ids <- unlist(formatted.samples, use.names = FALSE)

  # Legacy single-string format, e.g.  -"Sample1" -"Sample2"
  if (length(ids) == 1 && grepl("[[:space:]\"]|^-", ids)) {
    ids <- sub("^-", "", ids)          # strip only the leading formatting dash
    ids <- gsub(" -", " ", ids)        # collapse the " -" separators
    ids <- gsub("\"", "", ids)         # drop quotes
    ids <- unlist(strsplit(ids, "[[:space:]]+"))
  }

  # Common cleanup for both paths: trim whitespace, drop quotes/empties.
  ids <- trimws(gsub("\"", "", ids))
  ids <- ids[nzchar(ids)]

  if (length(ids) == 0) stop("No sample names could be parsed from the config 'samples' entry.")

  ids
}


# This got complicated as I was learning about how ggsave saves images (i.e. 
# be default it saves huge images) and how the pdf rendering presents those
# images (i.e. it makes them huge if they're save huge)
save_png_plot <- function(p, filename, ggplot.dir,
                          width = 7, height = 5, dpi = 150) {

	if (!dir.exists(ggplot.dir)) dir.create(ggplot.dir, recursive = TRUE, showWarnings = FALSE) 
	out <- file.path(ggplot.dir,  filename) 
	print(out)
	
	# works for ggplot or any grid grob via cowplot::ggdraw 
	if (inherits(p, "ggplot")) { 
		ggplot2::ggsave(out, plot = p, width = width, height = height, 
				units = "in", dpi = dpi, limitsize = TRUE) 
	} else if (!is.null(p$gtable)) { # pheatmap object 
		png(out, width = width, height = height, units = "in", res = dpi) 
		grid::grid.newpage(); grid::grid.draw(p$gtable); dev.off() 
	} else if (inherits(p, "grob") || inherits(p, "gTree")) { 
		png(out, width = width, height = height, units = "in", res = dpi) 
		grid::grid.newpage(); grid::grid.draw(p); dev.off() 
	} else { 
		# last resort: try plotting it 
		png(out, width = width, height = height, units = "in", res = dpi) 
		print(p); dev.off() 
	} 
	invisible(out) 
}


# Filter out the metrics we don't want, and format everything nicely
make_qc_table <- function(ids, dataset.path){
  
	d10x.metrics <- lapply(ids, function(sample.name){ 
				       metrics.path.cleaned <- file.path(dataset.path, sample.name, "metrics_summary.csv") 
				       metrics.path.raw     <- file.path(dataset.path, sample.name, "outs/metrics_summary.csv") 
				       tryCatch({ 
					       if (file.exists(metrics.path.cleaned)) { 
						       read.csv(metrics.path.cleaned, colClasses = "character") 
					       } else if (file.exists(metrics.path.raw)) { 
						       read.csv(metrics.path.raw, colClasses = "character") 
					       } else { 
						       data.frame() 
					       } 
				       }, error = function(cond){ 
					       message(paste0("Error loading the sample ", sample.name, ": ", conditionMessage(cond))) 
					       data.frame() 
				       }) 
				})

	# --- Make them a data.frame, keyed by sample --- 
	# Bind on a Sample column rather than rownames: rownames(df) <- ids throws
	# when a sample's metrics_summary.csv is missing (empty data.frame) or when
	# column sets differ between samples. Tag each frame with its sample first,
	# then bind by name so mismatched/absent samples degrade gracefully.
	names(d10x.metrics) <- ids
	d10x.metrics <- Filter(function(x) nrow(x) > 0, d10x.metrics)
	if (length(d10x.metrics) == 0) {
		warning("No metrics_summary.csv files could be read for any sample; returning empty QC table.")
		return(data.frame(metric = character(0)))
	}
	dropped <- setdiff(ids, names(d10x.metrics))
	if (length(dropped) > 0) {
		message("No metrics_summary.csv found for: ", paste(dropped, collapse = ", "),
			" — these samples are omitted from the QC metrics table.")
	}
	for (s in names(d10x.metrics)) d10x.metrics[[s]]$.sample <- s
	df <- dplyr::bind_rows(d10x.metrics)   # tolerant of differing columns
	rownames(df) <- df$.sample
	df$.sample <- NULL

	# --- Keep only desired metrics (only those actually present) --- 
	keep <- c("Estimated.Number.of.Cells", 
		  "Mean.Reads.per.Cell", 
		  "Median.Genes.per.Cell", 
		  "Valid.Barcodes") 
	keep_present <- intersect(keep, colnames(df))
	missing_cols <- setdiff(keep, colnames(df))
	if (length(missing_cols) > 0) {
		message("QC metrics not found in metrics_summary.csv (CellRanger version drift?): ",
			paste(missing_cols, collapse = ", "))
	}
	if (length(keep_present) == 0) {
		warning("None of the expected QC metric columns were found; returning empty QC table.")
		return(data.frame(metric = character(0)))
	}
	df_sub <- df[, keep_present, drop = FALSE] 
	
	# --- Prettify metric names --- 
	pretty_names <- function(x) { 
		s <- gsub("\\.", " ", x) 
		s <- tolower(s) 
		substr(s, 1, 1) <- toupper(substr(s, 1, 1)) 
		s 
	} 
	colnames(df_sub) <- pretty_names(colnames(df_sub))
  
 
	# --- Transpose; put 'metric' first --- 
	qc.table <- t(df_sub) 
	qc.table <- as.data.frame(qc.table, stringsAsFactors = FALSE) 
	qc.table$metric <- rownames(qc.table) 
	row.names(qc.table) <- NULL 
	qc.table <- qc.table[, c(ncol(qc.table), 1:(ncol(qc.table)-1))]
  
  	# --- Order metrics (optional) --- 
	metric_order <- c("Estimated number of cells", 
			  "Mean reads per cell", 
			  "Median genes per cell", 
			  "Valid barcodes") 
	qc.table <- qc.table %>% 
		mutate(metric = factor(metric, levels = metric_order)) %>% 
		arrange(metric) %>% 
		mutate(metric = as.character(metric))
	# return 
	qc.table
}


# This is mostly used for putting the sample names in the right order so things 
# are consistent in the plots
make_metric_table_from_list <- function(seurat.list, 
					metric, 
					sample_levels = NULL, 
					sort_alpha = TRUE) { 
	df <- do.call(rbind, lapply(names(seurat.list), function(s) { 
					    so <- seurat.list[[s]] 
					    stopifnot(metric %in% colnames(so[[]])) 
					    data.frame(Sample = s, value = so[[metric, drop = TRUE]], 
						       row.names = colnames(so), check.names = FALSE) 
				       })) 
	names(df)[2] <- metric 

	if (!is.null(sample_levels)) { 
		df$Sample <- factor(df$Sample, levels = sample_levels) 
	} else if (sort_alpha) { 
		df$Sample <- factor(df$Sample, levels = sort(unique(df$Sample))) 
	} else { 
		df$Sample <- factor(df$Sample, levels = unique(df$Sample))
  }
  
  df
}


make_metric_table <- function(obj, sample_col, metric, sample_levels = NULL, sort_alpha = TRUE) {
  stopifnot(metric %in% colnames(obj[[]]))
  out <- data.frame(
    Sample = obj[[sample_col, drop = TRUE]],
    value  = obj[[metric, drop = TRUE]],
    row.names = colnames(obj),
    check.names = FALSE
  )
  names(out)[2] <- metric

  if (!is.null(sample_levels)) {
    out$Sample <- factor(as.character(out$Sample), levels = sample_levels)
  } else if (sort_alpha) {
    out$Sample <- factor(as.character(out$Sample), levels = sort(unique(as.character(out$Sample))))
  } else {
    out$Sample <- factor(as.character(out$Sample), levels = unique(as.character(out$Sample)))
  }
  out
}


################################################################################
# Takes a DF and prints out a nice formatted table to the rendered PDF
################################################################################

render_formatted_table <- function(input.df) {
  if (!"metric" %in% names(input.df))
    stop("render_formatted_table() expects a 'metric' column ",
         "(it formats the CellRanger metrics table). Use render_plain_table() ",
         "for general data frames.")
  # ---- format numbers (no manual LaTeX) ----
  num_only <- function(x) as.numeric(gsub("[^0-9.]", "", x))
  
  fmt_one_cell <- function(metric, val_chr){
    if (is.na(val_chr) || val_chr == "") return(NA_character_)
    if (identical(metric, "Valid barcodes")){
      # CellRanger stores Valid Barcodes as a percentage integer (e.g. "93") so 
      # dividing by 100 converts it to a proportion for percent().
      v <- num_only(val_chr) / 100
      if (is.na(v)) return(NA_character_)
      percent(v, accuracy = 0.1)                # e.g. "93.8%"
    } else {
      v <- num_only(val_chr)
      if (is.na(v)) return(NA_character_)
      comma(v)                                   # e.g. "24,047"
    }
  }
  
  for (j in 2:ncol(input.df)) {
    input.df[[j]] <- mapply(fmt_one_cell, input.df$metric, input.df[[j]])
  }
  input.df[is.na(input.df)] <- ""
  input.df[] <- lapply(input.df, as.character)
  
  # ---- render (escape TRUE) ----
  input.df %>%
    rename(Metric = metric) %>%
    kable(format   = "latex",
          booktabs = TRUE,
          align    = c("l", rep("r", ncol(input.df)-1)),
          caption  = "Sequencing / QC metrics by sample",
          # escape = TRUE by default; keep it that way!
          longtable = FALSE) %>%
    kable_styling(
      latex_options = c("striped", "hold_position", "scale_down"),
      font_size = 12                                # body font
    ) %>%
    column_spec(1, bold = TRUE) %>%                 # bold Metric column
    add_header_above(c(" " = 1, "Samples" = ncol(input.df)-1)) %>%
    row_spec(0, angle = 45, font_size = 8)          # rotate/shrink header row (sample names)
}


# helper to locate 10x paths
read_10x_any <- function(sample.name, dataset.path, run.soupx, type = "auto") {
  # Define all possible paths
  raw_paths  <- file.path(dataset.path, sample.name, c("outs/raw_feature_bc_matrix", "raw_feature_bc_matrix"))
  filt_paths <- file.path(dataset.path, sample.name, c("outs/filtered_feature_bc_matrix", "filtered_feature_bc_matrix"))
  selected_path <- NULL
  # --- Logic Gate ---
  # 1. User specifically wants RAW
  if (type == "raw") {
    selected_path <- raw_paths[dir.exists(raw_paths)][1]
    if (is.na(selected_path)) stop("Raw matrix requested but not found for: ", sample.name)
    # 2. User specifically wants FILTERED
  } else if (type == "filtered") {
    selected_path <- filt_paths[dir.exists(filt_paths)][1]
    if (is.na(selected_path)) stop("Filtered matrix requested but not found for: ", sample.name)
    # 3. Default "auto" behavior: Try raw, then filtered
  } else {
    if (any(dir.exists(raw_paths))) {
      selected_path <- raw_paths[dir.exists(raw_paths)][1]
      message("Auto-detected: Using RAW matrix for ", sample.name)
    } else if (any(dir.exists(filt_paths))) {
      selected_path <- filt_paths[dir.exists(filt_paths)][1]
      message("Auto-detected: Using FILTERED matrix for ", sample.name)
    }
  }
  # --- Final Execution ---
  if (is.null(selected_path) || is.na(selected_path)) {
    stop("Could not locate any 10x matrix folders for sample: ", sample.name)
  }
  message("Loading data from: ", selected_path)
  
  # --- SoupX ---
  if (run.soupx) {
    # SoupX requires both raw and filtered matrices
    raw_path <- raw_paths[dir.exists(raw_paths)][1]
    filt_path <- filt_paths[dir.exists(filt_paths)][1]
    if (is.na(raw_path))  stop("SoupX requires a raw matrix, but none was found for: ",      sample.name, " at ", raw_path)
    if (is.na(filt_path)) stop("SoupX requires a filtered matrix, but none was found for: ", sample.name)
    message("Running SoupX for: ", sample.name)
    
    tod <- Read10X(data.dir = raw_path)   # Table of Droplets (all barcodes)
    toc <- Read10X(data.dir = filt_path)  # Table of Counts  (cell barcodes only)
    sc  <- SoupChannel(tod, toc)
    
    rm(tod); gc()  # free the large raw matrix immediately
    
    # downsample to keep memory manageable
    max.cells.for.clustering <- 5000
    if (ncol(toc) > max.cells.for.clustering) {
      set.seed(42)
      cells.to.use <- sample(colnames(toc), max.cells.for.clustering)
    } else {
      cells.to.use <- colnames(toc)
    }
    
    # Cluster cells within SoupX using a basic Seurat workflow, 
    # as SoupX needs cluster labels to estimate the contamination fraction
    tmp <- CreateSeuratObject(counts = toc)
    tmp <- NormalizeData(tmp,                          verbose = FALSE)
    tmp <- FindVariableFeatures(tmp, nfeatures = 1000, verbose = FALSE)
    tmp <- ScaleData(tmp,                              verbose = FALSE)
    tmp <- RunPCA(tmp, npcs = 10,                      verbose = FALSE)
    tmp <- FindNeighbors(tmp, dims = 1:10, k.param = 10, verbose = FALSE)
    tmp <- FindClusters(tmp, resolution = 0.3,         verbose = FALSE)
    sc  <- setClusters(sc, setNames(tmp$seurat_clusters, colnames(tmp)))
    
    sc  <- autoEstCont(sc, verbose = FALSE)
    return(adjustCounts(sc, roundToInt = TRUE))
  }
  
  Read10X(data.dir = selected_path)
}

vln_boxplot <- function(data, x, y, title) {
  p <- ggplot(data, aes(x = {{ x }}, y = {{ y }})) +
    geom_violin() +
    geom_boxplot(width = 0.1) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(title)

  p
}

shaded_vln_boxplot <- function(data, 
                               x, 
                               y, 
                               title, 
                               show.outliers = TRUE,
                               this_floor = -Inf, 
                               this_ceiling = Inf){
  ggplot(data, aes(x = {{x}}, y = {{y}})) +
    # shaded bands (put first so they sit behind the geoms)
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = this_ceiling, ymax = Inf,
             fill = "red", alpha = 0.08) +
    annotate("rect", xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = this_floor,
             fill = "red", alpha = 0.08) +
    geom_violin() +
    geom_boxplot(width = 0.1) +
    geom_hline(yintercept = c(this_floor, this_ceiling),
               color = "red", linetype = "dashed", linewidth = 0.6) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(title)
}


# Single source of truth for marker-plot filenames. Both the writer
# (feat_plots_top_genes, part 2) and the reader (part 3 summary) call this, so
# the two can never disagree on the suffix/assay/cluster naming scheme.
marker_plot_path <- function(ggplot.dir,
                             cluster.num,
                             kind = c("feat", "vln", "dot"),
                             assay.used = "SCT",
                             filename.suffix = "") {
  kind <- match.arg(kind)
  tag  <- switch(kind, feat = "featPlot", vln = "vlnPlot", dot = "dotPlot")
  fn   <- paste0("clustering_", filename.suffix,
                 "_marker_gene_", tag, "_cl_", cluster.num, "_", assay.used, ".png")
  file.path(ggplot.dir, fn)
}


# Assemble a list of per-gene panels into a titled grid and write it to disk.
# Shared by the feature- and violin-plot branches of feat_plots_top_genes().
save_panel_grid <- function(panels, title, path, num_cols, base_text,
                            width, height, dpi) {
  panels <- lapply(panels, function(p) {
    p + theme(
      axis.text.x = element_text(size = base_text, angle = 45, hjust = 1),
      axis.text.y = element_text(size = base_text),
      plot.title  = element_text(size = base_text + 1)
    )
  })
  grid  <- cowplot::plot_grid(plotlist = panels, ncol = num_cols)
  head  <- cowplot::ggdraw() + cowplot::draw_label(title, fontface = "bold",
                                                   size = base_text + 2)
  final <- cowplot::plot_grid(head, grid, ncol = 1, rel_heights = c(0.12, 1))
  ggplot2::ggsave(filename = path, plot = final, width = width, height = height,
                  units = "in", dpi = dpi, limitsize = TRUE)
  invisible(path)
}


# Marker plots for one cluster.
#
# `plots` selects which panels to draw; any subset of:
#   "feat" - UMAP feature plots, one panel per gene (grid)
#   "dot"  - single dot plot, genes x all clusters (mean expr + % expressing)
#   "vln"  - violin plots, one panel per gene (grid). Off by default: the dot
#            plot answers "is this gene cluster-specific" more compactly, but
#            violins show the within-cluster distribution, which the dot plot
#            hides. Turn on when you need to see bimodality or a marker driven
#            by a handful of high-expressing cells.
#
# Returns a named list of the paths actually written (invisibly), or NULL for a
# cluster with no markers.
feat_plots_top_genes <- function(cluster.num,
                                 seurat.obj,
                                 marker.genes.df,
                                 filename.suffix,
                                 ggplot.dir,
                                 plots = c("feat", "dot"),
                                 assay.used = "SCT",
                                 genes_per_page = 6,
                                 fig_width = 9,
                                 fig_height = 10,
                                 dot_height = 5,
                                 dpi = 150,
                                 num_cols = 3,
                                 umap_pt_size = 0.02,
                                 vln_pt_size = 0.05,
                                 base_text = 10,
                                 this_reduction = "umap.postint",
                                 to_raster = FALSE) {
  
  plots <- match.arg(plots, choices = c("feat", "vln", "dot"), several.ok = TRUE)
  if (!dir.exists(ggplot.dir)) dir.create(ggplot.dir, recursive = TRUE, showWarnings = FALSE)
  
  curr.markers <- marker.genes.df[marker.genes.df$cluster == cluster.num, ]
  genes <- as.character(na.omit(curr.markers$gene))
  genes <- genes[seq_len(min(length(genes), genes_per_page))]
  
  # A cluster can legitimately have no positive markers (only.pos = TRUE).
  if (length(genes) == 0) {
    message("Cluster ", cluster.num, ": no marker genes, skipping plots.")
    return(invisible(NULL))
  }
  
  out <- list()
  
  if ("feat" %in% plots) {
    out$feature_plot <- save_panel_grid(
      panels = Seurat::FeaturePlot(
        seurat.obj, features = genes, reduction = this_reduction,
        cols = c("grey85", "navy"), ncol = num_cols, pt.size = umap_pt_size,
        label.size = base_text, combine = FALSE,
        raster = to_raster, raster.dpi = c(dpi, dpi)),
      title  = paste0("Cluster ", cluster.num, " — top marker features"),
      path   = marker_plot_path(ggplot.dir, cluster.num, "feat", assay.used, filename.suffix),
      num_cols = num_cols, base_text = base_text,
      width = fig_width, height = fig_height, dpi = dpi)
  }
  
  if ("vln" %in% plots) {
    out$violin_plot <- save_panel_grid(
      panels = Seurat::VlnPlot(
        seurat.obj, features = genes, pt.size = vln_pt_size, combine = FALSE),
      title  = paste0("Cluster ", cluster.num, " — top marker violins"),
      path   = marker_plot_path(ggplot.dir, cluster.num, "vln", assay.used, filename.suffix),
      num_cols = num_cols, base_text = base_text,
      width = fig_width, height = fig_height, dpi = dpi)
  }
  
  if ("dot" %in% plots) {
    # DotPlot returns ONE ggplot covering all clusters, so no grid assembly.
    dot.plot.path <- marker_plot_path(ggplot.dir, cluster.num, "dot", assay.used, filename.suffix)
    dp <- Seurat::DotPlot(seurat.obj, features = genes, assay = assay.used) +
      ggtitle(paste0("Cluster ", cluster.num, " — top markers across all clusters")) +
      theme(axis.text.x = element_text(size = base_text, angle = 45, hjust = 1),
            axis.text.y = element_text(size = base_text),
            plot.title  = element_text(size = base_text + 1, face = "bold"))
    ggplot2::ggsave(filename = dot.plot.path, plot = dp,
                    width = fig_width, height = dot_height, units = "in",
                    dpi = dpi, limitsize = TRUE)
    out$dot_plot <- dot.plot.path
  }
  
  invisible(out)
}


  
################################################################################
# Part 2a / 2b split: handoff staleness + stability-based resolution selection
################################################################################

# The set of upstream parameters that, if changed, invalidate the integrated
# object produced by part 2a. Stamped into obj@misc by 2a and re-checked by 2b.
# Editing only the clustering resolution does NOT change this, so 2b can re-run
# freely; changing a QC cutoff DOES, forcing a 2a re-run.
upstream_params <- function(p2, run.doubletFinder, organism,
                            part1.suffix, part2.suffix) {
  list(
    mito_ceiling          = p2$mito_ceiling,
    RNA_count_floor       = p2$RNA_count_floor,
    RNA_count_ceiling     = p2$RNA_count_ceiling,
    feature_count_floor   = p2$feature_count_floor,
    feature_count_ceiling = p2$feature_count_ceiling,
    rbc_ceiling           = p2$rbc_ceiling,
    ribo_floor            = p2$ribo_floor,
    run_doubletFinder     = isTRUE(as.logical(run.doubletFinder)),
    organism              = organism,
    part1_suffix          = part1.suffix %||% "",
    part2_suffix          = part2.suffix %||% ""
  )
}

# Null-coalescing helper (rlang provides one, but keep this dependency-free).
`%||%` <- function(a, b) if (is.null(a)) b else a

# Compare the stamped upstream params to the current config; stop on mismatch.
stop_if_stale <- function(obj, current.params) {
  stamped <- obj@misc$upstream_params
  if (is.null(stamped)) {
    warning("Integrated object has no upstream_params stamp (built by an older ",
            "part 2a?). Proceeding, but re-run part 2a if results look off.")
    return(invisible(NULL))
  }
  if (!identical(stamped, current.params)) {
    changed <- names(current.params)[
      !mapply(identical, stamped[names(current.params)], current.params)]
    stop("Integrated object is STALE: upstream settings changed since part 2a ",
         "was run (", paste(changed, collapse = ", "), "). Re-run part 2a ",
         "before part 2b.")
  }
  invisible(TRUE)
}

# Subsample cluster-stability selection.
#
# Division of labour / attribution:
#   * subsampling + re-clustering loop .... this function (scclusteval expects
#     subsampling to be done outside R, in its Snakemake workflow)
#   * cluster-wise Jaccard stability metric  scclusteval::PairWiseJaccardSets()
#     (Tang et al. 2021, Bioinformatics; theory: Hennig 2007)
#   * decision rule ....................... chooseR-style (Patterson-Cross et al.
#     2021): bootstrap CI of the median per-cluster stability, threshold = the
#     highest lower-CI bound, choose the FINEST resolution whose median clears it.
#
# Per resolution:
#   * cluster the full graph once -> reference labels
#   * repeat n_subsample times: resample subsample_frac of cells, rebuild the
#     graph on the same reduction, re-cluster, and compare the re-clustering to
#     the reference labels RESTRICTED to those same cells, via
#     scclusteval::PairWiseJaccardSets(); take each reference cluster's best
#     Jaccard (row max)
#   * average across subsamples -> one robustness score PER reference cluster
#   * summarise that per-cluster distribution by its median + a percentile
#     bootstrap CI of the median
#
# Resolutions yielding <2 clusters are excluded (a 1-cluster "clustering" is
# trivially stable and would otherwise dominate the threshold).
#
# Returns list(summary = data.frame(resolution, n_clusters, median, low, high),
#              threshold = <numeric>, chosen_resolution = <numeric>).
select_resolution_by_stability <- function(obj,
                                            resolutions,
                                            graph.name = "int_snn", 
					    nn.graph.name = "int_nn", 
                                            reduction  = "integrated.rpca",
                                            dims       = 1:30,
                                            n_subsample = 5,
                                            subsample_frac = 0.8,
                                            boot_reps = 1000,
                                            conf = 0.95,
                                            seed = 1) {
  if (!requireNamespace("scclusteval", quietly = TRUE)) {
    stop("Package 'scclusteval' is required for stability selection. Install with:\n",
         "  renv::install(\"crazyhottommy/scclusteval\"); renv::snapshot()")
  }
  set.seed(seed)
  all.cells <- colnames(obj)

  # Percentile bootstrap CI of the median (dependency-free; chooseR uses a BCa
  # bootstrap via boot::, this is the lightweight equivalent).
  boot_median_ci <- function(x, conf = 0.95, R = 1000) {
    x <- x[is.finite(x)]
    if (length(x) < 2) {
      m <- if (length(x) == 1) x else NA_real_
      return(c(low = m, high = m))
    }
    meds <- vapply(seq_len(R), function(i) stats::median(sample(x, replace = TRUE)),
                   numeric(1))
    a <- (1 - conf) / 2
    q <- stats::quantile(meds, c(a, 1 - a), names = FALSE)
    c(low = q[1], high = q[2])
  }

  rows <- lapply(resolutions, function(res) {
    ref.obj <- FindClusters(obj, graph.name = graph.name,
                            resolution = res, verbose = FALSE)
    ref.lab <- stats::setNames(as.character(Idents(ref.obj)), colnames(ref.obj))
    ref.clusters <- sort(unique(ref.lab))

    # per-cluster best Jaccard accumulated across subsamples (NA where a ref
    # cluster is absent from a given subsample; averaged with na.rm)
    acc <- matrix(NA_real_, nrow = length(ref.clusters), ncol = n_subsample,
                  dimnames = list(ref.clusters, NULL))
    for (b in seq_len(n_subsample)) {
      cells.b <- sample(all.cells, floor(length(all.cells) * subsample_frac))
      sub <- subset(obj, cells = cells.b)
      sub <- FindNeighbors(sub, reduction = reduction, dims = dims,
                           graph.name = c(nn.graph.name, graph.name), verbose = FALSE)
      sub <- FindClusters(sub, graph.name = graph.name,
                          resolution = res, verbose = FALSE)
      test.lab <- stats::setNames(as.character(Idents(sub)), colnames(sub))

      # Compare on the SHARED cells: reference labels restricted to the subsample
      # vs the re-clustering. scclusteval returns a rows(ref) x cols(test) matrix
      # of pairwise Jaccard; row max = each reference cluster's best match.
      ref.sub <- ref.lab[names(test.lab)]
      jmat <- scclusteval::PairWiseJaccardSets(ref.sub, test.lab)
      best <- apply(jmat, 1, max)
      acc[names(best), b] <- best
    }
    cluster.stability <- rowMeans(acc, na.rm = TRUE)  # one score per ref cluster
    ci <- boot_median_ci(cluster.stability, conf = conf, R = boot_reps)

    data.frame(resolution = res,
               n_clusters = length(ref.clusters),
               median     = stats::median(cluster.stability, na.rm = TRUE),
               low        = unname(ci["low"]),
               high       = unname(ci["high"]),
               row.names  = NULL)
  })

  df <- do.call(rbind, rows)

  # Threshold + choice over NON-TRIVIAL resolutions only (>= 2 clusters).
  valid <- df[df$n_clusters >= 2 & is.finite(df$low), , drop = FALSE]
  if (nrow(valid) == 0) {
    return(list(summary = df, threshold = NA_real_,
                chosen_resolution = max(df$resolution)))
  }
  threshold <- max(valid$low)
  clears <- valid[is.finite(valid$median) & valid$median >= threshold, , drop = FALSE]
  chosen <- if (nrow(clears) > 0) max(clears$resolution) else
            valid$resolution[which.max(valid$median)]

  list(summary = df, threshold = threshold, chosen_resolution = chosen)
}

################################################################################
# Faceted QC plotting
#
# qc_long()          : stack the per-metric tables saved by part 1 into one
#                      long frame (Sample, metric, value [, stage])
# qc_thresholds()    : build the per-metric floor/ceiling frame from config
# qc_facet_plot()    : one violin+box figure, faceted by metric (and stage)
#
# These replace the per-metric vln_boxplot()/shaded_vln_boxplot() calls. The
# originals are kept for one-off single-metric plots elsewhere.
################################################################################

# Trim the longest common prefix off sample labels so 4+ samples stay legible
# inside narrow facets (SFN5_KO_KPC99_F_N1 -> F_N1). Cuts back to the last
# separator so tokens aren't split mid-word.
short_sample_labels <- function(x) {
  lv <- levels(factor(x))
  if (length(lv) < 2) return(factor(x))
  n <- min(nchar(lv)); i <- 0L
  while (i < n && length(unique(substr(lv, i + 1L, i + 1L))) == 1L) i <- i + 1L
  pre <- sub("[^_.\\-]*$", "", substr(lv[1], 1L, i))
  if (nchar(pre) == 0L) return(factor(x))
  k <- nchar(pre)
  factor(substr(as.character(x), k + 1L, nchar(as.character(x))),
         levels = substr(lv, k + 1L, nchar(lv)))
}

qc_long <- function(tables, stage = NULL, short.labels = TRUE) {
  stopifnot(is.list(tables), !is.null(names(tables)))
  out <- do.call(rbind, lapply(names(tables), function(nm) {
    d <- tables[[nm]]
    vcol <- setdiff(names(d), "Sample")[1]
    data.frame(Sample = as.character(d$Sample),
               metric = nm,
               value  = d[[vcol]],
               row.names = NULL, check.names = FALSE)
  }))
  out$metric <- factor(out$metric, levels = names(tables))
  out$Sample <- if (short.labels) short_sample_labels(out$Sample) else factor(out$Sample)
  if (!is.null(stage)) out$stage <- factor(stage, levels = c("raw", "filtered"))
  out
}

# spec: named list, e.g. list(nCount_RNA = c(lower = 500, upper = 25000),
#                            percent.mt  = c(upper = 5))
qc_thresholds <- function(spec, metric.levels) {
  do.call(rbind, lapply(names(spec), function(nm) {
    v <- spec[[nm]]
    grab <- function(k) if (is.na(v[k]) || is.null(v[k])) NA_real_ else as.numeric(v[k])
    data.frame(metric = factor(nm, levels = metric.levels),
               lower  = grab("lower"),
               upper  = grab("upper"))
  }))
}

qc_facet_plot <- function(dat,
                          thresholds    = NULL,
                          log.y         = FALSE,
                          facet.nrow    = 1,
                          title         = NULL,
                          drop.inactive = TRUE,
                          shade.stage   = "raw") {

  by.stage <- "stage" %in% names(dat)

  # A bound that excludes zero cells is not a filter. Blank it so red shading
  # always means "cells are removed here" (config ships rbc_ceiling = 100 and
  # ribo_floor = 0, neither of which cuts anything).
  if (!is.null(thresholds) && drop.inactive) {
    for (i in seq_len(nrow(thresholds))) {
      v <- dat$value[dat$metric == thresholds$metric[i]]
      if (!is.na(thresholds$lower[i]) && !any(v < thresholds$lower[i], na.rm = TRUE))
        thresholds$lower[i] <- NA_real_
      if (!is.na(thresholds$upper[i]) && !any(v > thresholds$upper[i], na.rm = TRUE))
        thresholds$upper[i] <- NA_real_
    }
  }

  p <- ggplot(dat, aes(x = Sample, y = value))

  if (!is.null(thresholds)) {
    lo <- thresholds[!is.na(thresholds$lower), , drop = FALSE]
    hi <- thresholds[!is.na(thresholds$upper), , drop = FALSE]
    if (by.stage) {
      st <- factor(shade.stage, levels = levels(dat$stage))
      if (nrow(lo)) lo$stage <- st
      if (nrow(hi)) hi$stage <- st
    }
    if (nrow(lo)) p <- p +
      geom_rect(data = lo, inherit.aes = FALSE,
                aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = lower),
                fill = "red", alpha = 0.08) +
      geom_hline(data = lo, inherit.aes = FALSE, aes(yintercept = lower),
                 colour = "red", linetype = "dashed", linewidth = 0.5)
    if (nrow(hi)) p <- p +
      geom_rect(data = hi, inherit.aes = FALSE,
                aes(xmin = -Inf, xmax = Inf, ymin = upper, ymax = Inf),
                fill = "red", alpha = 0.08) +
      geom_hline(data = hi, inherit.aes = FALSE, aes(yintercept = upper),
                 colour = "red", linetype = "dashed", linewidth = 0.5)
  }

  # scale = "width" is the fix for the current plots: a handful of cells at
  # percent.mt ~95 or nCount ~175k otherwise squash every violin into a sliver.
  p <- p +
    geom_violin(scale = "width", fill = NA, linewidth = 0.3) +
    geom_boxplot(width = 0.12, fill = NA, linewidth = 0.3,
                 outlier.size = 0.15, outlier.alpha = 0.25)

  p <- p + if (by.stage) {
    facet_grid(metric ~ stage, scales = "free_y", switch = "y")
  } else {
    facet_wrap(~ metric, scales = "free_y", nrow = facet.nrow)
  }

  if (log.y) p <- p + scale_y_log10(
    labels = scales::label_number(scale_cut = scales::cut_short_scale()))

  p +
    labs(title = title, x = NULL, y = NULL) +
    theme(axis.text.x  = element_text(angle = 45, hjust = 1, size = 7),
          strip.text   = element_text(size = 8),
          panel.spacing = unit(0.6, "lines"))
}

qc_attrition <- function(md, spec) {
  miss <- setdiff(names(spec), names(md))
  if (length(miss)) stop("metrics not found in metadata: ", paste(miss, collapse = ", "))

  fails <- vapply(names(spec), function(nm) {
    v <- md[[nm]]; b <- spec[[nm]]
    f <- rep(FALSE, length(v))
    # guard NA: v < bound yields NA, which propagates into colSums
    if (!is.na(b["lower"])) f <- f | (!is.na(v) & v < b["lower"])
    if (!is.na(b["upper"])) f <- f | (!is.na(v) & v > b["upper"])
    f
  }, logical(nrow(md)))

  n.crit <- rowSums(fails)
  n.fail <- colSums(fails)

  data.frame(
    criterion      = names(spec),
    # "min"/"max" instead of >=/<=: avoids LaTeX text-mode angle brackets
    bound          = vapply(spec, function(b) paste(c(
                       if (!is.na(b["lower"])) paste0("min ", b["lower"]),
                       if (!is.na(b["upper"])) paste0("max ", b["upper"])),
                       collapse = ", "), ""),
    cells_lost     = n.fail,
    lost_only_here = colSums(fails & n.crit == 1),
    pct_of_total   = round(100 * n.fail / nrow(md), 2),
    row.names = NULL, stringsAsFactors = FALSE)
}

render_plain_table <- function(df, caption = NULL, align = NULL,
                               col.names = NULL, font_size = 10) {
  if (is.null(align))     align     <- c("l", rep("r", ncol(df) - 1))
  if (is.null(col.names)) col.names <- gsub("_", " ", names(df))
  df %>%
    kable(format = "latex", booktabs = TRUE, align = align,
          caption = caption, col.names = col.names, longtable = FALSE) %>%
    kable_styling(latex_options = c("striped", "hold_position"),
                  font_size = font_size)
}
