# ASCAT.sc Explorer

An interactive browser-based viewer for single-cell copy-number data produced for [ASCAT.sc](https://github.com/VanLoo-lab/ASCAT.sc).

**Live app:** https://anagchar.github.io/ASCAT.sc-Explorer

---

## What it shows

- Genome-wide CN heatmap across all cells (zoomable, scrollable)
- Allele-specific view (major / minor copy numbers)
- Hierarchical dendrogram with ward.D2 clustering
- Per-cell profile panel (click any row)
- QC sidebar — filter cells by MAPD and median residual
- Cell-type colouring when annotations are supplied

---

## Workflow

### Step 1 — convert your RDS to JSON (in R)

```r
# install.packages(c("data.table", "jsonlite"))   # if not already installed

source("R/ascatsc_to_web.R")

# Minimal — total CN only
rds_to_web("your_results.rds", "ascat_data.json")

# With cell-type annotations from a multiome experiment
rds_to_web(
  rds_path         = "your_results.rds",
  output_path      = "ascat_data.json",
  cell_type_file   = "final_cell_type.txt",   # two columns: barcode_rna, cell_type
  barcode_map_file = "barcodes_atac_gex.csv"  # two columns: barcode_dna, barcode_rna
)
```

The script reads directly from the `ASCAT.sc` result object and exports
everything the app needs: CN profiles, allele-specific calls (nMajor/nMinor),
CI bands, dendrogram, quality metrics, and cell-type labels.

### Step 2 — upload the JSON

Go to **https://anagchar.github.io/ASCAT.sc-Explorer**, drag-and-drop (or
click to browse) your `ascat_data.json`, and explore.

---

## Running locally (no install required)

The repo includes a pre-built `build/` folder. Just clone and open:

```bash
git clone https://github.com/anagchar/ASCAT.sc-Explorer.git
# then open build/index.html in your browser
open ASCAT.sc-Explorer/build/index.html        # macOS
xdg-open ASCAT.sc-Explorer/build/index.html   # Linux
# on Windows: double-click build\index.html
```

No npm, no server required.

---

## JSON format reference

The app accepts any JSON file with at minimum these top-level keys:

| Key | Type | Description |
|-----|------|-------------|
| `bins` | object | `chr`, `start`, `end`, `start_cum`, `end_cum` arrays (one entry per bin) |
| `chr_info` | array | `{chr, start_cum, end_cum, mid_cum}` per chromosome |
| `profiles` | object | `{ cellName: [cn_bin1, cn_bin2, …] }` — integer total CN |

Optional keys: `nMajor`, `nMinor`, `ci`, `raw`, `quality`, `dendrogram`,
`clustering_order`, `cell_types`, `metadata`.

The app will recompute clustering in-browser if `dendrogram` is absent.

---

## Deploying your own instance

```bash
npm run deploy     # builds and pushes to gh-pages branch
```

Update the `homepage` field in `package.json` to match your own GitHub Pages URL before running.

---

## Dependencies

| | |
|---|---|
| **Browser app** | React 19, D3 v7, Tailwind CSS (CDN) |
| **R export** | data.table, jsonlite |

