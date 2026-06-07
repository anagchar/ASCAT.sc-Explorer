// Color palettes matching ASCATsc_plot.R exactly.
// Total CN: get_total_cn_colors() — 12 levels (0–10, 10+)
// Allele-specific: get_allele_specific_colors() — 21 named states
// Profile alleles: Major #E03546, Minor #1b38ae

// ── Total copy number ────────────────────────────────────────────────────────
// R: colors 0→"#153570", 1→"#577aba", 2→"#c1c1c1", 3→"#e3b55f", 4→"#d6804f",
//          5→"#b3402e", 6→"#821010", 7→"#6a0936", 8→"#ab1964", 9→"#b6519f",
//         10→"#ad80b9", "10+"→"#c2a9d1"
// The web app caps display at 6+ so we map 6–10+ all to the ≥6 bucket.
// NA/unknown gets a neutral grey.
export const CN_COLORS = {
  0:  "#153570",
  1:  "#577aba",
  2:  "#c1c1c1",
  3:  "#e3b55f",
  4:  "#d6804f",
  5:  "#b3402e",
  6:  "#821010",
  na: "#BEBEBE",
};

// Light-mode uses the exact same R colours — R itself renders on white,
// so these are already designed for light backgrounds.
export const CN_COLORS_LIGHT = {
  0:  "#153570",
  1:  "#577aba",
  2:  "#888888",
  3:  "#c49631",
  4:  "#d6804f",
  5:  "#b3402e",
  6:  "#821010",
  na: "#aaaaaa",
};

// ── Allele-specific states ───────────────────────────────────────────────────
// R: get_allele_specific_colors() — exact hex values copied verbatim
export const AS_COLORS = {
  "0|0": "#153570",
  "1|0": "#345d8e",
  "1|1": "#c1c1c1",
  "2|0": "#737373",
  "2|1": "#edd78e",
  "3|0": "#c49631",
  "2|2": "#a8d9a8",
  "3|1": "#5fa85f",
  "4|0": "#2d7a2d",
  "3|2": "#f2b4a8",
  "4|1": "#cf3a3a",
  "5|0": "#7d1111",
  "3|3": "#d4b896",
  "4|2": "#a67c52",
  "5|1": "#7a5230",
  "6|0": "#4a2c12",
  "4|3": "#c7aed5",
  "5|2": "#9b6bb5",
  "6|1": "#6e3a8e",
  "7|0": "#4a1570",
  "7+":  "#3d0d5c",
  "NA":  "#BEBEBE",
};

// ── Profile plot allele colors ───────────────────────────────────────────────
// R: plot_allele_profile() — Major "#E03546", Minor "#1b38ae"
export const MAJOR_COLOR = "#E03546";
export const MINOR_COLOR = "#1b38ae";

// ── Legend groups (for CN legend bar in heatmap) ─────────────────────────────
// Matches the 0–6 range the web app displays
export const CN_LEGEND_ENTRIES = [
  [0, "0"],
  [1, "1"],
  [2, "2"],
  [3, "3"],
  [4, "4"],
  [5, "5"],
  [6, "6+"],
];

// ── Helpers ──────────────────────────────────────────────────────────────────

export function cnColor(v, light = false) {
  const pal = light ? CN_COLORS_LIGHT : CN_COLORS;
  if (v == null || isNaN(v)) return pal.na;
  const r = Math.round(v);
  return pal[r <= 0 ? 0 : r >= 6 ? 6 : r] || pal.na;
}

function hexToRGB(h) {
  return [
    parseInt(h.slice(1, 3), 16),
    parseInt(h.slice(3, 5), 16),
    parseInt(h.slice(5, 7), 16),
  ];
}

function makeCnRGB(light = false) {
  const pal = light ? CN_COLORS_LIGHT : CN_COLORS;
  const c = {};
  Object.entries(pal).forEach(([k, h]) => { c[k] = hexToRGB(h); });
  return v => {
    if (v == null || isNaN(v)) return c.na;
    const r = Math.round(v);
    return c[r <= 0 ? "0" : r >= 6 ? "6" : String(r)] || c.na;
  };
}

export const cnRGBDark  = makeCnRGB(false);
export const cnRGBLight = makeCnRGB(true);
