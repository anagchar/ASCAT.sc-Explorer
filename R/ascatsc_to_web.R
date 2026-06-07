# =============================================================================
# ascatsc_to_web.R
# Single-file pipeline: ASCAT.sc .rds → JSON for the React web explorer
#
# Usage (minimal):
#   source("ascatsc_to_web.R")
#   rds_to_web("your_results.rds", "ascat_data.json")
#
# Usage (with cell type annotations):
#   source("ascatsc_to_web.R")
#   rds_to_web("your_results.rds", "ascat_data.json",
#              cell_type_file  = "final_cell_type.txt",
#              barcode_map_file = "barcodes_atac_gex.csv")
#
# Dependencies: data.table, jsonlite
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
})

# =============================================================================
# HELPER: map RNA cell type annotations to DNA barcodes
# =============================================================================

map_cell_types_to_dna <- function(cell_type_file, barcode_map_file) {
  ct   <- fread(cell_type_file,   header = TRUE)
  bmap <- fread(barcode_map_file, header = TRUE)
  setnames(ct,   c("barcode_rna", "cell_type"))
  setnames(bmap, c("barcode_dna", "barcode_rna"))
  merged <- merge(bmap, ct, by = "barcode_rna", all = FALSE)
  message(sprintf("  Mapped %d / %d RNA cell types to DNA barcodes",
                  nrow(merged), nrow(ct)))
  merged[, .(barcode = barcode_dna, cell_type)]
}

# =============================================================================
# HELPER: hierarchical clustering grouped by cell type
# =============================================================================

cluster_by_cell_type <- function(mat, cell_types,
                                 dist_method   = "manhattan",
                                 hclust_method = "ward.D2") {
  ct <- as.data.table(cell_types)
  setnames(ct, c("barcode", "cell_type"))

  available <- colnames(mat)
  ct <- ct[barcode %in% available]

  if (nrow(ct) == 0) stop("No barcodes in cell_types match the profile column names.")

  message(sprintf("  Matched %d / %d cells to cell type annotations",
                  nrow(ct), length(available)))

  # Cells with no annotation are kept as "Unknown" — not silently dropped
  unmatched <- setdiff(available, ct$barcode)
  if (length(unmatched) > 0) {
    message(sprintf("  %d cells have no annotation — labelled 'Unknown'", length(unmatched)))
    ct <- rbind(ct, data.table(barcode = unmatched, cell_type = "Unknown"))
  }

  ordered_cells <- character(0)
  dendrograms   <- list()

  for (ctype in sort(unique(ct$cell_type))) {
    group_cells <- ct[cell_type == ctype, barcode]
    if (length(group_cells) == 1) {
      ordered_cells       <- c(ordered_cells, group_cells)
      dendrograms[[ctype]] <- NULL
    } else {
      sub_mat <- mat[, group_cells, drop = FALSE]
      sub_mat[is.na(sub_mat)] <- -1
      hc <- hclust(dist(t(sub_mat), method = dist_method), method = hclust_method)
      ordered_cells       <- c(ordered_cells, labels(as.dendrogram(hc)))
      dendrograms[[ctype]] <- as.dendrogram(hc)
    }
  }

  list(
    ordered_cells = ordered_cells,
    cell_type_df  = ct[match(ordered_cells, barcode)],
    dendrograms   = dendrograms
  )
}

# =============================================================================
# HELPER: extract genomic bins from res object
# =============================================================================

extract_bins <- function(res) {
  message("Extracting genomic bins...")
  lSe_use    <- if (!is.null(res$nlSe)) res$nlSe else res$lSe
  use_plural <- !is.null(res$nlSe)

  bins <- rbindlist(lapply(names(lSe_use), function(chr) {
    data.table(
      chr   = chr,
      start = if (use_plural) lSe_use[[chr]]$starts else lSe_use[[chr]]$start,
      end   = if (use_plural) lSe_use[[chr]]$ends   else lSe_use[[chr]]$end
    )
  }))

  bins[, bin      := .I]
  bins[, end_cum  := cumsum((end - start) + 1)]
  bins[, start_cum := c(1, end_cum[.I - 1] + 1)]

  message(sprintf("  Found %d bins across %d chromosomes",
                  nrow(bins), length(unique(bins$chr))))
  bins
}

# =============================================================================
# HELPER: extract total CN profiles + raw per-bin logR-derived CN
# =============================================================================

extract_total_profiles <- function(res) {
  message("Extracting total CN profiles...")
  if (is.null(res$allProfiles)) stop("Cannot find allProfiles in result object")

  cell_names <- names(res$allProfiles)
  message(sprintf("  Processing %d cells...", length(cell_names)))

  profiles_list <- lapply(cell_names, function(cell) {
    dt <- as.data.table(res$allProfiles[[cell]])
    rep(dt$total_copy_number, as.numeric(dt$num.mark)) |> as.numeric()
  })
  profiles <- setnames(data.table(do.call(cbind, profiles_list)), cell_names)

  raw <- NULL
  if (!is.null(res$allTracks.processed) && !is.null(res$allSolutions)) {
    message("  Extracting raw per-bin CN values...")
    raw_list <- lapply(cell_names, function(cell) {
      dt <- rbindlist(res$allTracks.processed[[cell]]$lCTS)
      res$allSolutions[[cell]]$ploidy * 2^(dt$smoothed)
    })
    raw <- setnames(data.table(do.call(cbind, raw_list)), cell_names)
  } else {
    message("  Warning: allTracks.processed / allSolutions missing — no raw dots in profile plot")
  }

  list(profiles = profiles, raw = raw)
}

# =============================================================================
# HELPER: detect allele-specific mode
# =============================================================================

detect_allele_specific <- function(res) {
  if (!is.null(res$allProfiles_AS_smoothed) && length(res$allProfiles_AS_smoothed) > 0) return(TRUE)
  if (!is.null(res$allProfiles_AS)          && length(res$allProfiles_AS)          > 0) return(TRUE)
  first <- res$allProfiles[[names(res$allProfiles)[1]]]
  if (is.data.frame(first)) return(all(c("nMajor", "nMinor") %in% colnames(first)))
  FALSE
}

# =============================================================================
# HELPER: extract allele-specific profiles (nMajor / nMinor)
# =============================================================================

extract_allele_profiles <- function(res, bins) {
  message("Extracting allele-specific profiles...")
  as_profiles <- if (!is.null(res$allProfiles_AS_smoothed)) res$allProfiles_AS_smoothed else res$allProfiles_AS
  if (is.null(as_profiles)) stop("No allProfiles_AS or allProfiles_AS_smoothed found.")

  cell_names  <- names(as_profiles)
  n_bins      <- nrow(bins)
  is_smoothed <- !is.null(res$allProfiles_AS_smoothed)

  profiles_list <- lapply(cell_names, function(cell) {
    cell_data <- as_profiles[[cell]]
    prof <- if (is_smoothed) {
      if (!is.data.frame(cell_data) && !is.matrix(cell_data))
        return(list(nMajor = rep(NA_real_, n_bins), nMinor = rep(NA_real_, n_bins), total = rep(NA_real_, n_bins)))
      as.data.table(cell_data)
    } else {
      if (!is.list(cell_data) || is.null(cell_data$nprof.fixed))
        return(list(nMajor = rep(NA_real_, n_bins), nMinor = rep(NA_real_, n_bins), total = rep(NA_real_, n_bins)))
      as.data.table(cell_data$nprof.fixed)
    }
    num_mark <- as.numeric(as.data.table(res$allProfiles[[cell]])$num.mark)
    prof[is.na(nA), nA := 0]
    prof[is.na(nB), nB := total_copy_number]
    nA <- rep(as.numeric(prof$nA), num_mark)
    nB <- rep(as.numeric(prof$nB), num_mark)
    list(nMajor = pmax(nA, nB), nMinor = pmin(nA, nB), total = nA + nB)
  })
  names(profiles_list) <- cell_names

  nMajor <- setnames(data.table(do.call(cbind, lapply(profiles_list, `[[`, "nMajor"))), cell_names)
  nMinor <- setnames(data.table(do.call(cbind, lapply(profiles_list, `[[`, "nMinor"))), cell_names)
  total  <- setnames(data.table(do.call(cbind, lapply(profiles_list, `[[`, "total"))),  cell_names)

  # Drop cells where every bin is NA in both alleles
  fully_na <- sapply(nMajor, function(col) all(is.na(col))) &
              sapply(nMinor, function(col) all(is.na(col)))
  if (any(fully_na)) {
    dropped <- cell_names[fully_na]
    message(sprintf("  Dropping %d cells with no allele-specific call", length(dropped)))
    keep <- !fully_na
    nMajor <- nMajor[, keep, with = FALSE]
    nMinor <- nMinor[, keep, with = FALSE]
    total  <- total[,  keep, with = FALSE]
  }

  list(nMajor = nMajor, nMinor = nMinor, total = total)
}

# =============================================================================
# STEP 1: prepare_explorer_data()
# Processes res object into everything the exporter needs.
# =============================================================================

prepare_explorer_data <- function(res, cell_types = NULL) {
  message("=== Preparing ASCAT.sc data ===")

  bins     <- extract_bins(res)
  profiles <- extract_total_profiles(res)
  is_as    <- detect_allele_specific(res)
  message(sprintf("  Mode: %s", ifelse(is_as, "Allele-Specific", "Total CN")))

  allele_profiles <- if (is_as) extract_allele_profiles(res, bins) else NULL

  cell_names <- colnames(profiles$profiles)

  # --- Clustering ---
  message("  Computing hierarchical clustering...")
  bin_cols  <- c("chr", "start", "end", "bin", "start_cum", "end_cum")
  mat       <- as.matrix(data.table(cbind(bins, profiles$profiles))[, setdiff(names(data.table(cbind(bins, profiles$profiles))), bin_cols), with = FALSE])
  mat[is.na(mat)] <- -1

  cell_type_df <- NULL
  hc           <- NULL
  if (!is.null(cell_types)) {
    message("  Using cell type grouped clustering...")
    ct_result    <- cluster_by_cell_type(mat, cell_types)
    ordered_cells <- ct_result$ordered_cells
    cell_type_df  <- ct_result$cell_type_df
  } else {
    hc            <- hclust(dist(t(mat), method = "manhattan"), method = "ward.D2")
    ordered_cells <- hc$labels[hc$order]
  }

  # --- Quality metrics ---
  message("  Computing per-cell quality metrics (mapd_logR)...")

  # mapd_logR: within-chromosome median absolute pairwise difference on
  # GC-corrected logR. Cross-chromosome diffs excluded (they span boundaries).
  compute_mapd_logR <- function(logR_vals, chr_vals) {
    if (length(logR_vals) != length(chr_vals)) return(NA_real_)
    diffs <- unlist(tapply(logR_vals, chr_vals, function(v) {
      v <- v[!is.na(v)]
      if (length(v) < 2) return(numeric(0))
      abs(diff(v))
    }), use.names = FALSE)
    if (length(diffs) == 0) return(NA_real_)
    median(diffs, na.rm = TRUE)
  }

  cell_quality <- sapply(cell_names, function(cell_id) {
    raw    <- profiles$raw[[cell_id]]
    fitted <- profiles$profiles[[cell_id]]
    if (is.null(raw)) return(NA_real_)
    median(abs(raw - fitted), na.rm = TRUE)
  })

  cell_mapd <- sapply(cell_names, function(cell_id) {
    if (is.null(res$allTracks.processed[[cell_id]])) return(NA_real_)
    logR_vals <- rbindlist(res$allTracks.processed[[cell_id]]$lCTS)$smoothed
    if (length(logR_vals) != nrow(bins)) return(NA_real_)
    compute_mapd_logR(logR_vals, bins$chr)
  })

  quality_dt <- data.table(
    cell            = cell_names,
    median_residual = as.numeric(cell_quality),
    mapd            = as.numeric(cell_mapd)
  )

  message(sprintf("  Thresholds — residual: %.3f  mapd: %.3f",
                  median(quality_dt$median_residual, na.rm = TRUE) + 2 * mad(quality_dt$median_residual, na.rm = TRUE),
                  median(quality_dt$mapd,            na.rm = TRUE) + 2 * mad(quality_dt$mapd,            na.rm = TRUE)))
  message(sprintf("  Ready: %d cells, %d bins", length(cell_names), nrow(bins)))

  list(
    res                = res,
    bins               = bins,
    profiles           = profiles,
    allele_profiles    = allele_profiles,
    is_allele_specific = is_as,
    ordered_cells      = ordered_cells,
    quality_dt         = quality_dt,
    hc                 = hc,
    cell_type_df       = cell_type_df
  )
}

# =============================================================================
# STEP 2: export_for_web()
# Serialises app_data to JSON.
# =============================================================================

export_for_web <- function(app_data,
                           output_path = "ascat_data.json",
                           pretty      = FALSE) {
  message("=== Exporting to JSON ===")

  bins       <- app_data$bins
  profiles   <- app_data$profiles
  is_as      <- app_data$is_allele_specific
  cell_names <- colnames(profiles$profiles)
  n_cells    <- length(cell_names)
  n_bins     <- nrow(bins)

  message(sprintf("  Cells: %d  |  Bins: %d  |  Mode: %s",
                  n_cells, n_bins, ifelse(is_as, "Allele-Specific", "Total CN")))

  # 1. metadata
  metadata <- list(n_cells = n_cells, n_bins = n_bins, is_allele_specific = is_as)

  # 2. bins
  bins_out <- list(
    chr       = as.character(bins$chr),
    start     = as.integer(bins$start),
    end       = as.integer(bins$end),
    start_cum = as.numeric(bins$start_cum),
    end_cum   = as.numeric(bins$end_cum)
  )

  # 3. chr_info
  chr_info_dt <- bins[, .(start_cum = min(start_cum), end_cum = max(end_cum),
                           mid_cum = (min(start_cum) + max(end_cum)) / 2), by = chr]
  chr_info <- lapply(seq_len(nrow(chr_info_dt)), function(i)
    list(chr = chr_info_dt$chr[i], start_cum = chr_info_dt$start_cum[i],
         end_cum = chr_info_dt$end_cum[i], mid_cum = chr_info_dt$mid_cum[i]))

  # 4. profiles
  message("  Serializing CN profiles...")
  total_profiles <- lapply(as.list(profiles$profiles), function(v) as.integer(round(v)))

  # 5. raw
  raw_out <- NULL
  if (!is.null(profiles$raw)) {
    message("  Serializing raw profiles...")
    raw_out <- lapply(as.list(profiles$raw), function(v) round(as.numeric(v), 4))
  }

  # 6. quality
  message("  Serializing quality metrics...")
  qdt <- app_data$quality_dt
  res <- app_data$res
  coverage_vals <- NULL
  if (!is.null(res$allSolutions)) {
    coverage_vals <- tryCatch(
      sapply(cell_names, function(cell) {
        sol <- res$allSolutions[[cell]]
        if (!is.null(sol$coverage)) as.numeric(sol$coverage) else NA_real_
      }), error = function(e) NULL)
  }
  quality_out <- setNames(lapply(cell_names, function(cell) {
    idx <- which(qdt$cell == cell)
    q <- list(
      median_residual = if (length(idx) > 0) round(qdt$median_residual[idx[1]], 6) else NA,
      mapd            = if (length(idx) > 0) round(qdt$mapd[idx[1]],            6) else NA
    )
    if (!is.null(coverage_vals) && !is.na(coverage_vals[cell]))
      q$coverage <- round(coverage_vals[cell], 2)
    q
  }), cell_names)

  # 7. clustering_order
  clustering_order <- app_data$ordered_cells

  # 8. nMajor / nMinor
  nMajor_out <- nMinor_out <- NULL
  if (is_as && !is.null(app_data$allele_profiles)) {
    message("  Serializing allele-specific profiles...")
    nMajor_out <- lapply(as.list(app_data$allele_profiles$nMajor), as.integer)
    nMinor_out <- lapply(as.list(app_data$allele_profiles$nMinor), as.integer)
  }

  # 9. CI bands
  ci_out <- NULL
  if (is_as && !is.null(res)) {
    as_profiles <- if (!is.null(res$allProfiles_AS_smoothed)) res$allProfiles_AS_smoothed else res$allProfiles_AS
    if (!is.null(as_profiles) && !is.null(res$allProfiles)) {
      message("  Serializing CI bands...")
      is_smoothed <- !is.null(res$allProfiles_AS_smoothed)
      ci_list <- lapply(cell_names, function(cell) tryCatch({
        cell_data  <- as_profiles[[cell]]
        total_prof <- as.data.table(res$allProfiles[[cell]])
        num_mark   <- as.numeric(total_prof$num.mark)
        seg_prof   <- if (is_smoothed) {
          if (!is.data.frame(cell_data) && !is.matrix(cell_data)) return(NULL)
          as.data.table(cell_data)
        } else {
          if (!is.list(cell_data) || is.null(cell_data$nprof.fixed)) return(NULL)
          as.data.table(cell_data$nprof.fixed)
        }
        if (!all(c("q05", "q95", "total_copy_number") %in% names(seg_prof))) return(NULL)
        total_cn   <- as.numeric(seg_prof$total_copy_number)
        baf_lower  <- as.numeric(seg_prof$q05)
        baf_upper  <- as.numeric(seg_prof$q95)
        list(
          nMajor_lower = round(rep(baf_lower * total_cn,          num_mark), 4),
          nMajor_upper = round(rep(baf_upper * total_cn,          num_mark), 4),
          nMinor_lower = round(rep(total_cn - baf_upper * total_cn, num_mark), 4),
          nMinor_upper = round(rep(total_cn - baf_lower * total_cn, num_mark), 4)
        )
      }, error = function(e) NULL))
      names(ci_list) <- cell_names
      ci_list <- Filter(Negate(is.null), ci_list)
      if (length(ci_list) > 0) ci_out <- ci_list
    }
  }

  # 10. dendrogram
  dendro_out <- NULL
  if (!is.null(app_data$hc)) {
    message("  Serializing dendrogram...")
    hc <- app_data$hc
    dendro_out <- list(
      merge  = lapply(seq_len(nrow(hc$merge)), function(i) as.integer(hc$merge[i, ])),
      height = round(as.numeric(hc$height), 6)
    )
  }

  # 11. cell_types
  cell_types_out <- NULL
  if (!is.null(app_data$cell_type_df)) {
    ct <- app_data$cell_type_df
    cell_types_out <- setNames(as.list(ct$cell_type), ct$barcode)
  }

  # Assemble payload
  message("  Assembling payload...")
  payload <- list(metadata = metadata, bins = bins_out, chr_info = chr_info,
                  profiles = total_profiles, quality = quality_out,
                  clustering_order = clustering_order)
  if (!is.null(raw_out))        payload$raw        <- raw_out
  if (!is.null(nMajor_out))     payload$nMajor     <- nMajor_out
  if (!is.null(nMinor_out))     payload$nMinor     <- nMinor_out
  if (!is.null(ci_out))         payload$ci         <- ci_out
  if (!is.null(dendro_out))     payload$dendrogram <- dendro_out
  if (!is.null(cell_types_out)) payload$cell_types <- cell_types_out

  message(sprintf("  Writing → %s", output_path))
  jsonlite::write_json(payload, output_path, auto_unbox = TRUE,
                       digits = 6, pretty = pretty, null = "null", na = "null")

  size_mb <- file.info(output_path)$size / 1024^2
  message(sprintf("  Done. %.1f MB  |  %d cells  |  allele-specific: %s",
                  size_mb, n_cells, is_as))
  if (!is.null(cell_types_out))
    message(sprintf("  Cell types: %d unique labels", length(unique(unlist(cell_types_out)))))
  message("=== Export complete ===")

  invisible(output_path)
}

# =============================================================================
# TOP-LEVEL WRAPPER: rds_to_web()
# One call does everything.
# =============================================================================

rds_to_web <- function(rds_path,
                       output_path      = "ascat_data.json",
                       cell_type_file   = NULL,
                       barcode_map_file = NULL,
                       pretty           = FALSE) {

  message(sprintf("Loading %s ...", rds_path))
  res <- readRDS(rds_path)

  cell_types <- NULL
  if (!is.null(cell_type_file) && !is.null(barcode_map_file)) {
    message("Mapping cell type annotations...")
    cell_types <- map_cell_types_to_dna(cell_type_file, barcode_map_file)
  }

  app_data <- prepare_explorer_data(res, cell_types = cell_types)
  export_for_web(app_data, output_path = output_path, pretty = pretty)
}
