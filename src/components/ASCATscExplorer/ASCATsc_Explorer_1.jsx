import { useState, useEffect, useRef, useMemo, useCallback, useDeferredValue, memo, forwardRef, useImperativeHandle } from "react";
import { AS_COLORS } from "../../constants/colors";
import { generateDemoData } from "./demoData";
import { axes, chrLabel, css, esc, fit, fmt, fmtBp, fmtCount, fmtInt, FONT, hexRGB, niceTicks, quantile, saveCanvas, tip } from "./chart";
import { computeQCMetrics, fmtMetric, isOutlier, METRIC_META, metricSummaries, num } from "./qc";
import QCTab, { useWidth as useWidthOf } from "./QCTab";

/* Precomputed RGB lookup for allele-specific states (mirrors R's get_allele_specific_colors) */
function hexToRGBArr(h) {
  return [parseInt(h.slice(1,3),16), parseInt(h.slice(3,5),16), parseInt(h.slice(5,7),16)];
}
const AS_RGB = {};
for (const [state, hex] of Object.entries(AS_COLORS)) AS_RGB[state] = hexToRGBArr(hex);

/* Circlize-style G-banding stain colours */
const STAIN_COLORS = {
  gneg:    "#FFFFFF",
  gpos25:  "#C8C8C8",
  gpos50:  "#A0A0A0",
  gpos75:  "#787878",
  gpos100: "#000000",
  acen:    "#CC0000",
  gvar:    "#A0A0A0",
  stalk:   "#5C5C9C",
};

/* ═══════════════════════════════════════════════════════════════════════════
   DENDROGRAM HELPERS
   ═══════════════════════════════════════════════════════════════════════════ */
function layoutDendrogram(merge, height, nLeaves, cellMapping = null) {
  const segments = [];
  const nodeY = new Float64Array(merge.length);
  const nodeX = new Float64Array(merge.length);
  nodeY.fill(NaN);
  const getY = idx => {
    const leafIdx = idx < 0 ? -idx - 1 : null;
    if (leafIdx !== null) {
      if (!cellMapping) return leafIdx + 0.5;
      const displayY = cellMapping[leafIdx];
      return Number.isFinite(displayY) && displayY >= 0 ? displayY + 0.5 : null;
    }
    const y = nodeY[idx - 1];
    return Number.isFinite(y) ? y : null;
  };
  const getX = idx => idx < 0 ? 0 : nodeX[idx - 1];
  for (let i = 0; i < merge.length; i++) {
    const [left, right] = merge[i];
    const ly = getY(left), ry = getY(right);
    const lx = getX(left), rx = getX(right);
    const h = height[i];
    nodeX[i] = h;

    if (ly == null && ry == null) continue;
    if (ly != null && ry != null) {
      nodeY[i] = (ly + ry) / 2;
      segments.push({ x1: lx, y1: ly, x2: h, y2: ly });
      segments.push({ x1: rx, y1: ry, x2: h, y2: ry });
      segments.push({ x1: h, y1: ly, x2: h, y2: ry });
      continue;
    }

    // Pruned branch: one child is filtered out, keep the surviving branch connected.
    if (ly != null) {
      nodeY[i] = ly;
      segments.push({ x1: lx, y1: ly, x2: h, y2: ly });
    } else {
      nodeY[i] = ry;
      segments.push({ x1: rx, y1: ry, x2: h, y2: ry });
    }
  }
  let maxH = 0;
  for (let i = 0; i < height.length; i++) if (height[i] > maxH) maxH = height[i];
  return { segments, maxH };
}

function buildDendrogramFromProfiles(profiles, nMajor = null, nMinor = null) {
  const cells = Object.keys(profiles || {});
  const n = cells.length;
  if (n < 2) return null;

  // Mirror R: combined <- rbind(maj_mat, min_mat); combined[is.na(combined)] <- -1
  // Each cell's vector is [maj0..majN, min0..minN] (all major first, then all minor).
  // NA values become -1 exactly as in the R code.
  const useAS = nMajor && nMinor;
  const vectors = cells.map(c => {
    if (useAS && nMajor[c] && nMinor[c]) {
      const maj = nMajor[c], min = nMinor[c];
      const len = maj.length;
      const v = new Array(len * 2);
      for (let i = 0; i < len; i++) v[i]       = maj[i] ?? -1;
      for (let i = 0; i < len; i++) v[len + i]  = min[i] ?? -1;
      return v;
    }
    return profiles[c] || [];
  });
  const maxNodes = 2 * n - 1;
  const dist = new Float64Array(maxNodes * maxNodes);
  const active = new Uint8Array(maxNodes);
  const nodeCode = new Int32Array(maxNodes);
  const nodeSize = new Int32Array(maxNodes);

  const idx = (i, j) => i * maxNodes + j;
  const setDist = (i, j, d) => { dist[idx(i, j)] = d; dist[idx(j, i)] = d; };
  const getDist = (i, j) => dist[idx(i, j)];

  for (let i = 0; i < n; i++) {
    active[i] = 1;
    nodeCode[i] = -(i + 1);
    nodeSize[i] = 1;
  }

  // Store squared Manhattan distances — Ward.D2 Lance-Williams updates operate
  // on squared distances throughout; heights are emitted as sqrt at merge time.
  for (let i = 0; i < n; i++) {
    const a = vectors[i];
    for (let j = i + 1; j < n; j++) {
      const b = vectors[j];
      const len = Math.max(a.length, b.length);
      let manhattan = 0;
      // NA → -1 matches R: mat[is.na(mat)] <- -1
      for (let k = 0; k < len; k++) manhattan += Math.abs((a[k] ?? -1) - (b[k] ?? -1));
      setDist(i, j, manhattan * manhattan);
    }
  }

  const merge = [];
  const height = [];
  const children = new Map();
  let nextNode = n;
  let activeCount = n;

  while (activeCount > 1) {
    let minD = Infinity;
    let s = -1;
    let t = -1;
    for (let i = 0; i < nextNode; i++) {
      if (!active[i]) continue;
      for (let j = i + 1; j < nextNode; j++) {
        if (!active[j]) continue;
        const d = getDist(i, j);
        if (d < minD) { minD = d; s = i; t = j; }
      }
    }
    if (s < 0 || t < 0 || !Number.isFinite(minD)) break;

    const stepCode = merge.length + 1;
    const codeS = nodeCode[s];
    const codeT = nodeCode[t];
    merge.push([codeS, codeT]);
    // Height = sqrt of squared distance, matching R's ward.D2 convention
    height.push(Math.sqrt(minD));
    children.set(stepCode, [codeS, codeT]);

    const u = nextNode++;
    nodeCode[u] = stepCode;
    nodeSize[u] = nodeSize[s] + nodeSize[t];

    // Ward.D2 Lance-Williams update on squared distances
    const dST = getDist(s, t);
    for (let v = 0; v < u; v++) {
      if (!active[v] || v === s || v === t) continue;
      const dSV = getDist(s, v);
      const dTV = getDist(t, v);
      const nS = nodeSize[s], nT = nodeSize[t], nV = nodeSize[v];
      const tot = nS + nT + nV;
      const d2 = ((nV + nS) / tot) * dSV + ((nV + nT) / tot) * dTV - (nV / tot) * dST;
      setDist(u, v, Math.max(0, d2));
    }

    active[s] = 0;
    active[t] = 0;
    active[u] = 1;
    activeCount -= 1;
  }

  if (merge.length !== n - 1) return null;

  const rootCode = merge.length;
  const order = [];
  const stack = [rootCode];
  while (stack.length) {
    const code = stack.pop();
    if (code < 0) {
      order.push(cells[-code - 1]);
      continue;
    }
    const kids = children.get(code);
    if (!kids) continue;
    stack.push(kids[1]);
    stack.push(kids[0]);
  }

  return { merge, height, order };
}

function normalizeLoadedData(rawData) {
  const data = { ...rawData };
  const profileCells = Object.keys(data.profiles || {});

  const mergeRaw = data?.dendrogram?.merge ?? data?.merge ?? data?.hclust?.merge;
  const heightRaw = data?.dendrogram?.height ?? data?.height ?? data?.hclust?.height;
  const orderRaw = data?.clustering_order ?? data?.order ?? data?.dendrogram?.order ?? data?.hclust?.order;

  const merge = Array.isArray(mergeRaw)
    ? mergeRaw
      .map(row => Array.isArray(row) && row.length >= 2 ? [Number(row[0]), Number(row[1])] : null)
      .filter(Boolean)
    : null;
  const height = Array.isArray(heightRaw) ? heightRaw.map(v => Number(v)).filter(v => Number.isFinite(v)) : null;

  const hasValidDendrogram = !!(
    merge &&
    height &&
    merge.length > 0 &&
    height.length === merge.length
  );

  const hasAS = data.nMajor && data.nMinor && Object.keys(data.nMajor).length > 0;

  if (hasValidDendrogram) {
    data.dendrogram = {
      ...(typeof data.dendrogram === "object" && data.dendrogram ? data.dendrogram : {}),
      merge,
      height,
    };
    // Build a separate AS dendrogram from allele-specific vectors when no
    // AS dendrogram was supplied in the JSON (the uploaded data only ever has
    // one dendrogram, which was computed from total CN).
    if (hasAS && !data.dendrogram_as) {
      const computedAS = buildDendrogramFromProfiles(data.profiles, data.nMajor, data.nMinor);
      if (computedAS) {
        data.dendrogram_as = { merge: computedAS.merge, height: computedAS.height };
        data.clustering_order_as = computedAS.order;
      }
    }
  } else if (profileCells.length > 1) {
    const computed = buildDendrogramFromProfiles(data.profiles);
    if (computed) {
      data.dendrogram = {
        ...(typeof data.dendrogram === "object" && data.dendrogram ? data.dendrogram : {}),
        merge: computed.merge,
        height: computed.height,
      };
      if (!Array.isArray(data.clustering_order) || data.clustering_order.length === 0) {
        data.clustering_order = computed.order;
      }
    }
    // Build AS dendrogram and store its leaf order separately
    if (hasAS) {
      const computedAS = buildDendrogramFromProfiles(data.profiles, data.nMajor, data.nMinor);
      if (computedAS) {
        data.dendrogram_as = { merge: computedAS.merge, height: computedAS.height };
        data.clustering_order_as = computedAS.order;
      }
    }
  }

  if (Array.isArray(orderRaw) && orderRaw.length > 0 && profileCells.length > 0) {
    const numericOrder = orderRaw.map(v => Number(v));
    const allNumeric = numericOrder.every(v => Number.isFinite(v));
    let orderCells = null;

    if (allNumeric) {
      const oneBased = numericOrder.every(v => v >= 1 && v <= profileCells.length);
      const zeroBased = numericOrder.every(v => v >= 0 && v < profileCells.length);
      if (oneBased) orderCells = numericOrder.map(v => profileCells[v - 1]);
      else if (zeroBased) orderCells = numericOrder.map(v => profileCells[v]);
    }

    const candidateOrder = orderCells || orderRaw;
    const filteredOrder = candidateOrder
      .map(v => String(v))
      .filter((cell, idx, arr) => profileCells.includes(cell) && arr.indexOf(cell) === idx);

    if (filteredOrder.length > 0) data.clustering_order = filteredOrder;
  }

  if (!Array.isArray(data.clustering_order) || data.clustering_order.length === 0) {
    data.clustering_order = profileCells;
  }

  // Compute bins_with_cna from profiles for any cell that doesn't already have it
  if (data.profiles && data.quality) {
    for (const cell of profileCells) {
      if (data.quality[cell] && data.quality[cell].bins_with_cna == null) {
        const profile = data.profiles[cell];
        if (Array.isArray(profile)) {
          data.quality[cell] = {
            ...data.quality[cell],
            bins_with_cna: profile.filter(v => Math.round(v) !== 2).length,
          };
        }
      }
    }
  }

  // Compute the Gini coefficient (coverage evenness) from per-bin read counts
  // for any cell that doesn't already carry a precomputed value.
  if (data.reads && data.quality) {
    for (const cell of profileCells) {
      if (data.quality[cell] && data.quality[cell].gini == null && data.reads[cell]) {
        const gini = computeGini(data.reads[cell]);
        if (gini != null) data.quality[cell] = { ...data.quality[cell], gini };
      }
    }
  }

  return data;
}

/* ═══════════════════════════════════════════════════════════════════════════
   CYTOBANDS  (hg38 fetched lazily from UCSC public API)
   ═══════════════════════════════════════════════════════════════════════════ */
const _cytoBandsCache = { data: null };
function useCytobands() {
  const [bands, setBands] = useState(_cytoBandsCache.data);
  useEffect(() => {
    if (_cytoBandsCache.data) { setBands(_cytoBandsCache.data); return; }
    fetch("https://api.genome.ucsc.edu/getData/track?genome=hg38;track=cytoBandIdeo")
      .then(r => r.json())
      .then(json => {
        const raw = json.cytoBandIdeo || [];
        const parsed = raw.map(b => ({
          chr: b.chrom, start: +b.chromStart, end: +b.chromEnd,
          name: b.name, stain: b.gieStain,
        }));
        _cytoBandsCache.data = parsed;
        setBands(parsed);
      })
      .catch(() => {}); // silently degrade — ideogram simply won't render
  }, []);
  return bands; // null until loaded
}

/* ═══════════════════════════════════════════════════════════════════════════
   SEGMENT SIZE FILTER UTILITY
   Returns a new profile where any run of consecutive non-diploid (CN≠2) bins
   shorter than minBins is replaced with CN=2.
   ═══════════════════════════════════════════════════════════════════════════ */
function filterSmallSegments(profile, minBins) {
  if (!minBins || minBins <= 1 || !profile) return profile;
  const result = profile.slice();
  let i = 0;
  while (i < result.length) {
    const cn = Math.round(result[i]);
    if (cn === 2) { i++; continue; }
    let j = i + 1;
    while (j < result.length && Math.round(result[j]) === cn) j++;
    if (j - i < minBins) {
      for (let k = i; k < j; k++) result[k] = 2;
    }
    i = j;
  }
  return result;
}

/* ═══════════════════════════════════════════════════════════════════════════
   GINI COEFFICIENT
   Standard discrete Gini coefficient over per-bin read counts — measures how
   unevenly reads are spread across the genome (0 = perfectly uniform coverage,
   →1 = all reads concentrated in a few bins). Also powers the Lorenz curve.
   ═══════════════════════════════════════════════════════════════════════════ */
function computeGini(values) {
  if (!values) return null;
  const v = values.filter(x => x != null && Number.isFinite(x) && x >= 0).sort((a, b) => a - b);
  const n = v.length;
  if (n === 0) return null;
  const sum = v.reduce((a, b) => a + b, 0);
  if (sum === 0) return 0;
  let weighted = 0;
  for (let i = 0; i < n; i++) weighted += (i + 1) * v[i];
  return (2 * weighted) / (n * sum) - (n + 1) / n;
}

const PUBLIC = process.env.PUBLIC_URL || "";
const LOGO_LIGHT = PUBLIC + "/logo_light.svg", LOGO_DARK = PUBLIC + "/logo_dark.svg";

/* Total copy-number colours come from the theme's --cn0..--cn6 tokens */
function cnPalette() {
  const rgb = [0, 1, 2, 3, 4, 5, 6].map((k) => hexRGB(css("--cn" + k))), na = hexRGB(css("--cnna"));
  return (v) => (v == null || isNaN(v) ? na : rgb[v <= 0 ? 0 : v >= 6 ? 6 : Math.round(v)]);
}
const asKey = (a, b) => {
  if (a == null || b == null) return null;
  const M = Math.max(a, b), m = Math.min(a, b);
  return M + m > 7 ? "7+" : AS_COLORS[`${M}|${m}`] ? `${M}|${m}` : "7+";
};
const swatch = (c) => `<span class="swatch" style="background:${c}"></span>`;

/** Bin index whose cumulative span holds gx (binary search). */
function binAt(bins, gx) {
  let lo = 0, hi = bins.end_cum.length - 1;
  while (lo < hi) { const mid = (lo + hi) >> 1; if (bins.end_cum[mid] < gx) lo = mid + 1; else hi = mid; }
  return lo;
}

/* ═══════════════════════════════════════════════════════════════════════════
   UPLOAD SCREEN
   ═══════════════════════════════════════════════════════════════════════════ */
async function readText(file) {
  if (/\.gz$/i.test(file.name) && typeof DecompressionStream !== "undefined")
    return new Response(file.stream().pipeThrough(new DecompressionStream("gzip"))).text();
  return file.text();
}

function UploadScreen({ onLoad, onToggleTheme }) {
  const [dragging, setDragging] = useState(false), [error, setError] = useState(null), [loading, setLoading] = useState(false);
  const fileRef = useRef();
  const handleFile = async (file) => {
    if (!file) return;
    if (/\.rds$/i.test(file.name)) {
      setError(`This is an .rds file — convert it to JSON first.\n\nIn R:\n  source("R/ascatsc_to_web.R")\n  rds_to_web("${file.name}", "ascat_data.json")\n\nThen open the resulting ascat_data.json here.`);
      return;
    }
    setLoading(true); setError(null);
    try {
      const d = JSON.parse(await readText(file));
      if (!d.bins || !d.profiles) throw new Error("Invalid JSON: missing 'bins' or 'profiles'.");
      if (!d.chr_info) throw new Error("Missing 'chr_info'.");
      onLoad(normalizeLoadedData(d), file.name.replace(/\.json(\.gz)?$/i, ""));
    } catch (err) { setError(err.message); setLoading(false); }
  };
  return (
    <div className="upload">
      <button className="iconbtn" style={{ position: "fixed", top: 16, right: 16 }} onClick={onToggleTheme} title="Light or dark (t)" aria-label="Switch between light and dark"><ThemeIcon /></button>
      <div className="inner">
        <div className="hero">
          <div className="logos"><img className="logo-light" src={LOGO_LIGHT} alt="ASCAT.sc" /><img className="logo-dark" src={LOGO_DARK} alt="ASCAT.sc" /></div>
          <h1>ASCAT.sc Explorer</h1>
          <p>Single-cell copy number, interactively</p>
        </div>
        <button className={"card drop" + (dragging ? " over" : "")}
          onDragOver={(e) => { e.preventDefault(); setDragging(true); }} onDragLeave={() => setDragging(false)}
          onDrop={(e) => { e.preventDefault(); setDragging(false); handleFile(e.dataTransfer.files[0]); }} onClick={() => fileRef.current?.click()}>
          <input ref={fileRef} type="file" accept=".json,.gz,.rds" hidden onChange={(e) => handleFile(e.target.files[0])} />
          {loading ? <><div className="spinner" /><div className="small">Opening…</div></> : <>
            <svg width="34" height="34" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="1.5"><path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4M7 10l5-5 5 5M12 5v12" /></svg>
            <div className="big">Drop <code>ascat_data.json</code> here</div>
            <div className="small">or click to choose a file · .json or .json.gz</div>
          </>}
        </button>
        {error && <div className="err">{error}</div>}
        <div className="actions"><button className="iconbtn" onClick={() => onLoad(generateDemoData(), "Demo dataset")}>Load the demo dataset · 200 cells, allele-specific</button></div>
        <div className="card">
          <div className="hd"><h2>Export from R</h2><span className="sub">from an ASCAT.sc result</span></div>
          <div className="bd"><pre className="mono">{`source("R/ascatsc_to_web.R")\nrds_to_web("your_results.rds", "ascat_data.json")`}</pre></div>
        </div>
      </div>
    </div>
  );
}

/* ═══════════════════════════════════════════════════════════════════════════
   DENDROGRAM CANVAS
   ═══════════════════════════════════════════════════════════════════════════ */
// Rows of the heatmap start below the chromosome labels and the ideogram
const HM_TOP = 40, HM_BOT = 4, LABEL_Y = 13, IDEO_Y = 22, IDEO_H = 12;

const DendrogramCanvas = memo(forwardRef(function DendrogramCanvas({ data, cellOrder, height, yRange, width = 92, theme, alleleMode = false }, ref) {
  const canvasRef = useRef();
  useImperativeHandle(ref, () => canvasRef.current, []);
  // In allele-specific mode, prefer the AS dendrogram so that LOH-driven
  // clustering differences are reflected in the dendrogram arms.
  const dendro = (alleleMode && data.dendrogram_as) ? data.dendrogram_as : data.dendrogram;
  const hasValidDendro = Array.isArray(dendro?.merge) && Array.isArray(dendro?.height) && dendro.merge.length > 0 && dendro.height.length === dendro.merge.length;
  const nCells = cellOrder.length;
  const plotH = height - HM_TOP - HM_BOT;
  const PAD = 6;

  // Build cellMapping and extract subtree when zoomed.
  // R merge uses -(profileCells index + 1) for leaves.
  const { segments, maxH } = useMemo(() => {
    if (!hasValidDendro || nCells < 2) return { segments: [], maxH: 0 };
    const profileCells = Object.keys(data.profiles || {});
    const posByCell = new Map(cellOrder.map((cell, i) => [cell, i]));
    const cellMapping = {};
    profileCells.forEach((cell, i) => {
      const pos = posByCell.get(cell);
      cellMapping[i] = Number.isFinite(pos) ? pos : null;
    });

    const yr = yRange || [0, nCells];
    const isZoomed = yr[0] > 0 || yr[1] < nCells;

    if (!isZoomed) return layoutDendrogram(dendro.merge, dendro.height, nCells, cellMapping);

    // Find which profileCell indices are visible in the zoomed row window.
    const visibleLeafSet = new Set();
    profileCells.forEach((cell, i) => {
      const pos = cellMapping[i];
      if (pos != null && pos >= yr[0] - 0.5 && pos < yr[1] - 0.5) visibleLeafSet.add(i);
    });
    if (visibleLeafSet.size < 2) return layoutDendrogram(dendro.merge, dendro.height, nCells, cellMapping);

    // For each merge node, count how many visible leaves it contains.
    const nodeVisCount = new Int32Array(dendro.merge.length);
    for (let i = 0; i < dendro.merge.length; i++) {
      const [l, r] = dendro.merge[i];
      const lc = l < 0 ? (visibleLeafSet.has(-l - 1) ? 1 : 0) : nodeVisCount[l - 1];
      const rc = r < 0 ? (visibleLeafSet.has(-r - 1) ? 1 : 0) : nodeVisCount[r - 1];
      nodeVisCount[i] = lc + rc;
    }

    // The subtree root: the lowest merge node holding every visible leaf (their MRCA).
    let subtreeRoot = dendro.merge.length - 1;
    for (let i = 0; i < dendro.merge.length; i++) {
      if (nodeVisCount[i] === visibleLeafSet.size) { subtreeRoot = i; break; }
    }

    // Collect only the merge rows that belong to the subtree rooted at subtreeRoot.
    const subtreeSteps = new Set();
    const stack = [subtreeRoot + 1]; // 1-based step code
    while (stack.length) {
      const code = stack.pop();
      if (code < 0) continue; // leaf
      const idx = code - 1;
      subtreeSteps.add(idx);
      const [l, r] = dendro.merge[idx];
      if (l > 0) stack.push(l);
      if (r > 0) stack.push(r);
    }

    // Filtered merge/height of the subtree, step codes remapped to new 1-based
    // indices, in topological order (children before parents, as hclust).
    const oldToNew = new Int32Array(dendro.merge.length).fill(-1);
    const newMerge = [];
    const newHeight = [];
    for (let i = 0; i < dendro.merge.length; i++) {
      if (!subtreeSteps.has(i)) continue;
      const newIdx = newMerge.length;
      oldToNew[i] = newIdx;
      const [l, r] = dendro.merge[i];
      const newL = l < 0 ? l : (oldToNew[l - 1] >= 0 ? oldToNew[l - 1] + 1 : l);
      const newR = r < 0 ? r : (oldToNew[r - 1] >= 0 ? oldToNew[r - 1] + 1 : r);
      newMerge.push([newL, newR]);
      newHeight.push(dendro.height[i]);
    }

    return layoutDendrogram(newMerge, newHeight, visibleLeafSet.size, cellMapping);
  }, [dendro, hasValidDendro, nCells, cellOrder, data.profiles, yRange]);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas || !hasValidDendro || nCells < 2 || maxH === 0) return;
    const ctx = fit(canvas, width, height);
    ctx.clearRect(0, 0, width, height);
    const yr = yRange || [0, nCells];
    const ySpan = yr[1] - yr[0] || 1;
    // sqrt-scaled heights, as the viewer: the many low merges stay readable
    const xScale = h => width - PAD - Math.sqrt(h / maxH) * (width - PAD * 2);
    const yScale = ci => HM_TOP + ((ci - yr[0]) / ySpan) * plotH;

    ctx.strokeStyle = css("--ink3");
    ctx.lineWidth = 0.8;
    ctx.beginPath();
    for (const s of segments) {
      const minY = Math.min(s.y1, s.y2);
      const maxY = Math.max(s.y1, s.y2);
      if (maxY < yr[0] || minY > yr[1]) continue;
      // Horizontal segments draw as-is; vertical segments clamp to viewport
      const drawY1 = s.y1 === s.y2 ? s.y1 : Math.max(yr[0], Math.min(yr[1], s.y1));
      const drawY2 = s.y1 === s.y2 ? s.y2 : Math.max(yr[0], Math.min(yr[1], s.y2));
      ctx.moveTo(xScale(s.x1), yScale(drawY1));
      ctx.lineTo(xScale(s.x2), yScale(drawY2));
    }
    ctx.stroke();
  }, [segments, maxH, hasValidDendro, nCells, width, height, yRange, plotH, theme]);

  if (!hasValidDendro || nCells < 2) return null;
  return <canvas ref={canvasRef} style={{ width, height, flexShrink: 0 }} />;
}));

/* ═══════════════════════════════════════════════════════════════════════════
   HEATMAP (2-layer canvas: data layer + overlay layer)
   ═══════════════════════════════════════════════════════════════════════════ */
const HeatmapPanel = memo(forwardRef(function HeatmapPanel({ data, cellOrder, chrInfo, selectedCell, onCellClick, height = 500, alleleMode = false, zoom, onZoomChange, showDendro = false, theme }, ref) {
  const DENDRO_W = showDendro && data.dendrogram ? 92 : 0;
  const nCells = cellOrder.length;
  const yRange = useMemo(() => zoom?.y || [0, nCells], [zoom, nCells]);
  const cytobands = useCytobands();
  const dendroRef = useRef();
  const heatmapRef = useRef();

  useImperativeHandle(ref, () => ({
    download(basename = "heatmap", format = "png") { saveCanvas([DENDRO_W ? dendroRef.current : null, heatmapRef.current], basename, format); },
  }), [DENDRO_W]);

  return (
    <div className="hmwrap" style={{ display: "flex", width: "100%" }}>
      {DENDRO_W > 0 && (
        <DendrogramCanvas ref={dendroRef} data={data} cellOrder={cellOrder} height={height} yRange={yRange} width={DENDRO_W} theme={theme} alleleMode={alleleMode} />
      )}
      <div style={{ flex: 1, minWidth: 0 }}>
        <HeatmapCanvas ref={heatmapRef} data={data} cellOrder={cellOrder} chrInfo={chrInfo}
          selectedCell={selectedCell} onCellClick={onCellClick} height={height}
          alleleMode={alleleMode} zoom={zoom} onZoomChange={onZoomChange} theme={theme} cytobands={cytobands} />
      </div>
    </div>
  );
}));

const HeatmapCanvas = memo(forwardRef(function HeatmapCanvas({ data, cellOrder, chrInfo, selectedCell, onCellClick, height, alleleMode, zoom, onZoomChange, theme, cytobands }, ref) {
  const dataCanvasRef = useRef();
  const overlayCanvasRef = useRef();
  useImperativeHandle(ref, () => dataCanvasRef.current, []);
  const containerRef = useRef();
  const cw = useWidthOf(containerRef);
  const brushRef = useRef(null);
  const isDragging = useRef(false);
  const prevOverlaySize = useRef({ w: 0, h: 0 });

  const nBins = data.bins.chr.length;
  const nCells = cellOrder.length;
  const plotW = cw;
  const plotH = height - HM_TOP - HM_BOT;
  const genomeMax = useMemo(() => nBins > 0 ? data.bins.end_cum[nBins - 1] : 1, [data, nBins]);

  const xRange = useMemo(() => zoom?.x || [0, genomeMax], [zoom, genomeMax]);
  const yRange = useMemo(() => zoom?.y || [0, nCells],    [zoom, nCells]);
  const xSpan = xRange[1] - xRange[0] || 1;
  const ySpan = yRange[1] - yRange[0] || 1;

  const MTT = HM_TOP;
  const toPxX = useCallback(gx => ((gx - xRange[0]) / xSpan) * plotW, [xRange, xSpan, plotW]);
  const toPxY = useCallback(ci => MTT + ((ci - yRange[0]) / ySpan) * plotH, [MTT, yRange, ySpan, plotH]);
  const toGx  = useCallback(px => xRange[0] + (px / plotW) * xSpan, [xRange, xSpan, plotW]);
  const toCi  = useCallback(py => yRange[0] + ((py - MTT) / plotH) * ySpan, [MTT, yRange, ySpan, plotH]);

  // ── DATA LAYER ──
  useEffect(() => {
    const canvas = dataCanvasRef.current;
    if (!canvas) return;
    const ctx = fit(canvas, cw, height);
    ctx.clearRect(0, 0, cw, height);
    if (nCells === 0) {
      ctx.fillStyle = css("--ink3"); ctx.font = `13px ${FONT}`; ctx.fillText("No cell passes the filters.", 12, MTT + 30);
      return;
    }
    const cnRGB = cnPalette(), asNA = hexRGB(css("--cnna"));
    const dark = document.documentElement.dataset.theme === "dark";

    // ── Ideogram (cytobands) ──
    if (cytobands && cytobands.length > 0) {
      const chrOffsets = {};
      chrInfo.forEach(ci => { chrOffsets[ci.chr] = ci.start_cum; });
      ctx.fillStyle = dark ? css("--field") : "#ffffff";
      ctx.fillRect(0, IDEO_Y, plotW, IDEO_H);
      for (const band of cytobands) {
        const offset = chrOffsets[band.chr];
        if (offset === undefined) continue;
        const cumS = offset + band.start;
        const cumE = offset + band.end;
        if (cumE <= xRange[0] || cumS >= xRange[1]) continue;
        const px0 = Math.max(0, Math.floor(toPxX(Math.max(cumS, xRange[0]))));
        const px1 = Math.min(plotW, Math.ceil(toPxX(Math.min(cumE, xRange[1]))));
        if (px1 <= px0) continue;
        ctx.fillStyle = band.stain === "gneg" && dark ? css("--field") : STAIN_COLORS[band.stain] || "#AAAAAA";
        ctx.fillRect(px0, IDEO_Y, px1 - px0, IDEO_H);
      }
      ctx.strokeStyle = css("--line2");
      ctx.lineWidth = 0.5;
      chrInfo.forEach(({ start_cum, end_cum }) => {
        if (end_cum <= xRange[0] || start_cum >= xRange[1]) return;
        const px0 = Math.max(0, Math.floor(toPxX(Math.max(start_cum, xRange[0]))));
        const px1 = Math.min(plotW, Math.ceil(toPxX(Math.min(end_cum, xRange[1]))));
        if (px1 > px0 + 1) ctx.strokeRect(px0 + 0.5, IDEO_Y + 0.5, px1 - px0 - 1, IDEO_H - 1);
      });
    }

    // ── Visible ranges ──
    const ciStart = Math.max(0, Math.floor(yRange[0]));
    const ciEnd   = Math.min(nCells, Math.ceil(yRange[1]));
    let biStart = 0, biEnd = nBins;
    for (let i = 0; i < nBins; i++)         { if (data.bins.end_cum[i]   >= xRange[0]) { biStart = i;   break; } }
    for (let i = nBins - 1; i >= 0; i--)   { if (data.bins.start_cum[i] <= xRange[1]) { biEnd   = i+1; break; } }

    // ── Heatmap tiles — ImageData pixel buffer (no per-bin draw calls) ──
    const tileW = Math.ceil(toPxX(data.bins.end_cum[biEnd - 1])) - Math.floor(toPxX(data.bins.start_cum[biStart]));
    const tileH = Math.ceil(toPxY(ciEnd)) - Math.floor(toPxY(ciStart));
    if (tileW > 0 && tileH > 0) {
      const tileX0 = Math.floor(toPxX(data.bins.start_cum[biStart]));
      const tileY0 = Math.floor(toPxY(ciStart));
      const offscreen = new OffscreenCanvas(tileW, tileH);
      const octx = offscreen.getContext("2d");
      const imgData = octx.createImageData(tileW, tileH);
      const buf = imgData.data;

      for (let ci = ciStart; ci < ciEnd; ci++) {
        const cell = cellOrder[ci];
        const rowY0 = Math.floor(toPxY(ci))     - tileY0;
        const rowY1 = Math.ceil(toPxY(ci + 1))  - tileY0;
        const rowH  = Math.max(0, Math.min(tileH, rowY1) - Math.max(0, rowY0));
        if (rowH < 1) continue;
        const startRow = Math.max(0, rowY0);

        for (let bi = biStart; bi < biEnd; bi++) {
          const colX0 = Math.floor(toPxX(data.bins.start_cum[bi])) - tileX0;
          const colX1 = Math.ceil(toPxX(data.bins.end_cum[bi]))    - tileX0;
          const colW  = Math.max(0, Math.min(tileW, colX1) - Math.max(0, colX0));
          if (colW < 1) continue;
          const startCol = Math.max(0, colX0);

          let r, g, b;
          if (alleleMode && data.nMajor?.[cell] && data.nMinor?.[cell]) {
            const k = asKey(data.nMajor[cell][bi], data.nMinor[cell][bi]);
            [r, g, b] = k ? AS_RGB[k] : asNA;
          } else {
            const prof = data.profiles[cell]; if (!prof) continue;
            [r, g, b] = cnRGB(prof[bi]);
          }

          // Fill first row of this bin
          const rowBase = (startRow * tileW + startCol) * 4;
          for (let dx = 0; dx < colW; dx++) {
            const idx = rowBase + dx * 4;
            buf[idx] = r; buf[idx + 1] = g; buf[idx + 2] = b; buf[idx + 3] = 255;
          }
          // Copy first row to remaining rows of this cell
          for (let dy = 1; dy < rowH; dy++) {
            const dstStart = ((startRow + dy) * tileW + startCol) * 4;
            buf.copyWithin(dstStart, rowBase, rowBase + colW * 4);
          }
        }
      }
      octx.putImageData(imgData, 0, 0);
      ctx.drawImage(offscreen, tileX0, tileY0);
    }

    // ── Chromosome separators (surface-coloured gaps) and labels ──
    ctx.fillStyle = css("--hm-gap");
    chrInfo.forEach(({ start_cum }) => {
      if (start_cum <= xRange[0] || start_cum > xRange[1]) return;
      ctx.fillRect(toPxX(start_cum) - 0.75, MTT, 1.5, plotH);
    });
    ctx.fillStyle = css("--ink2"); ctx.font = `11px ${FONT}`; ctx.textAlign = "center"; ctx.textBaseline = "alphabetic";
    let lastLabel = -1e9;
    chrInfo.forEach(({ chr, start_cum, end_cum }) => {
      if (end_cum < xRange[0] || start_cum > xRange[1]) return;
      const a = toPxX(Math.max(start_cum, xRange[0])), b = toPxX(Math.min(end_cum, xRange[1])), mid = (a + b) / 2;
      if (b - a > 14 && mid - lastLabel > 22) { ctx.fillText(chrLabel(chr), mid, LABEL_Y); lastLabel = mid; }
    });
  }, [data, cellOrder, cw, height, alleleMode, theme, xRange, yRange, nBins, nCells, plotW, plotH, chrInfo, toPxX, toPxY, MTT, cytobands]);

  // ── OVERLAY LAYER ──
  const drawOverlay = useCallback(() => {
    const canvas = overlayCanvasRef.current;
    if (!canvas) return;
    const dpr = Math.max(1, window.devicePixelRatio || 1);
    if (prevOverlaySize.current.w !== cw || prevOverlaySize.current.h !== height) {
      canvas.width = cw * dpr; canvas.height = height * dpr;
      canvas.style.width = cw + "px"; canvas.style.height = height + "px";
      prevOverlaySize.current = { w: cw, h: height };
    }
    const ctx = canvas.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, cw, height);

    if (selectedCell && nCells > 0) {
      const idx = cellOrder.indexOf(selectedCell);
      if (idx >= 0) {
        const y0 = toPxY(idx), y1 = toPxY(idx + 1);
        if (y1 > MTT && y0 < MTT + plotH) {
          const hh = Math.max(2, y1 - y0);
          ctx.strokeStyle = css("--accent"); ctx.lineWidth = 2;
          ctx.strokeRect(1, (y0 + y1) / 2 - hh / 2 - 1, plotW - 2, hh + 2);
        }
      }
    }
    const b = brushRef.current;
    if (b && isDragging.current) {
      const x0 = Math.min(b.sx, b.cx), y0 = Math.min(b.sy, b.cy);
      const w = Math.abs(b.cx - b.sx), h = Math.abs(b.cy - b.sy);
      ctx.fillStyle = css("--accent-soft");
      ctx.strokeStyle = css("--accent"); ctx.lineWidth = 1;
      ctx.fillRect(x0, y0, w, h); ctx.strokeRect(x0, y0, w, h);
    }
  }, [cw, height, selectedCell, cellOrder, nCells, toPxY, plotW, plotH, MTT]);

  useEffect(() => { drawOverlay(); }, [drawOverlay, theme]);

  const getXY = e => {
    const rect = overlayCanvasRef.current.getBoundingClientRect();
    return { x: e.clientX - rect.left, y: e.clientY - rect.top };
  };
  const hover = (e, p) => {
    const ci = Math.floor(toCi(p.y));
    if (p.y < MTT || ci < 0 || ci >= nCells) { tip(e, null); return; }
    const cell = cellOrder[ci], bi = binAt(data.bins, toGx(p.x));
    let val;
    if (alleleMode && data.nMajor?.[cell]) {
      const k = asKey(data.nMajor[cell][bi], data.nMinor[cell][bi]);
      val = k ? `${swatch(AS_COLORS[k])} ${k}` : "no value";
    } else {
      const v = data.profiles[cell]?.[bi];
      val = v == null ? "no value" : `${swatch(css("--cn" + Math.max(0, Math.min(6, Math.round(v)))))} ${v}`;
    }
    const ct = data.cell_types?.[cell];
    tip(e, `<b class="mono">${esc(cell)}</b>${ct ? ` <span class="k">· ${esc(ct)}</span>` : ""}<br><span class="k">${chrLabel(data.bins.chr[bi])}:${fmtBp(data.bins.start[bi])}–${fmtBp(data.bins.end[bi])}</span><br>${alleleMode ? "major|minor" : "copy number"} <b>${val}</b>`);
  };
  const handleMouseDown = e => {
    if (e.button !== 0) return;
    const p = getXY(e);
    brushRef.current = { sx: p.x, sy: p.y, cx: p.x, cy: p.y };
    isDragging.current = true;
    tip(e, null);
  };
  const handleMouseMove = e => {
    const p = getXY(e);
    if (!isDragging.current) { hover(e, p); return; }
    brushRef.current.cx = p.x; brushRef.current.cy = p.y;
    drawOverlay();
  };
  const handleMouseUp = e => {
    if (!isDragging.current) return;
    isDragging.current = false;
    const p = getXY(e);
    const b = brushRef.current;
    if (!b) return;
    const dx = Math.abs(b.sx - p.x), dy = Math.abs(b.sy - p.y);
    if (dx < 5 && dy < 5) {
      const ci = Math.floor(toCi(p.y));
      if (ci >= 0 && ci < nCells) onCellClick(cellOrder[ci]);
    } else if (dx > 8 || dy > 8) {
      const gx0 = toGx(Math.min(b.sx, p.x)), gx1 = toGx(Math.max(b.sx, p.x));
      const cy0 = toCi(Math.max(MTT, Math.min(b.sy, p.y))), cy1 = toCi(Math.max(b.sy, p.y));
      onZoomChange({
        x: [Math.max(0, gx0), Math.min(genomeMax, gx1)],
        y: [Math.max(0, cy0), Math.min(nCells, cy1)],
      });
    }
    brushRef.current = null;
    drawOverlay();
  };

  const handleWheel = useCallback(e => {
    const yZoomed = yRange[0] !== 0 || yRange[1] !== nCells;
    const xZoomed = xRange[0] !== 0 || xRange[1] !== genomeMax;
    if (!yZoomed && !xZoomed) return;
    e.preventDefault();

    let newX = xRange, newY = yRange;

    if (e.deltaY !== 0 && yZoomed) {
      const span = yRange[1] - yRange[0];
      const step = span * 0.1 * Math.sign(e.deltaY);
      const y0 = Math.max(0, Math.min(nCells - span, yRange[0] + step));
      newY = [y0, y0 + span];
    }

    if (e.deltaX !== 0 && xZoomed) {
      const span = xRange[1] - xRange[0];
      const step = span * 0.1 * Math.sign(e.deltaX);
      const x0 = Math.max(0, Math.min(genomeMax - span, xRange[0] + step));
      newX = [x0, x0 + span];
    }

    onZoomChange({ x: newX, y: newY });
  }, [yRange, xRange, nCells, genomeMax, onZoomChange]);

  // Attach wheel listener as non-passive so preventDefault works
  useEffect(() => {
    const canvas = overlayCanvasRef.current;
    if (!canvas) return;
    canvas.addEventListener("wheel", handleWheel, { passive: false });
    return () => canvas.removeEventListener("wheel", handleWheel);
  }, [handleWheel]);

  return (
    <div ref={containerRef} style={{ width: "100%", position: "relative" }}>
      <canvas ref={dataCanvasRef} style={{ width: cw, height }} />
      <canvas ref={overlayCanvasRef}
        style={{ width: cw, height, position: "absolute", top: 0, left: 0, cursor: "crosshair" }}
        onMouseDown={handleMouseDown} onMouseMove={handleMouseMove} onMouseUp={handleMouseUp}
        onMouseLeave={e => { tip(e, null); if (isDragging.current) { isDragging.current = false; brushRef.current = null; drawOverlay(); } }}
        onDoubleClick={() => onZoomChange(null)} />
    </div>
  );
}));

/* ═══════════════════════════════════════════════════════════════════════════
   PROFILE PLOT (canvas) — total CN & allele-specific, following the heatmap's
   genomic zoom. Read depth draws underneath on the same x scale.
   ═══════════════════════════════════════════════════════════════════════════ */
const PR = { l: 40, t: 12, r: 10, b: 26 };

function drawChromosomes(g, data, xS, gx0, gx1, top, h, labelY) {
  g.textAlign = "center"; g.textBaseline = "alphabetic"; g.font = `11px ${FONT}`;
  let last = -1e9;
  for (const { chr, start_cum, end_cum } of data.chr_info) {
    if (end_cum < gx0 || start_cum > gx1) continue;
    const a = xS(Math.max(start_cum, gx0)), b = xS(Math.min(end_cum, gx1));
    if (start_cum > gx0) { g.strokeStyle = css("--line2"); g.lineWidth = 1; g.beginPath(); g.moveTo(Math.round(a) + 0.5, top); g.lineTo(Math.round(a) + 0.5, top + h); g.stroke(); }
    const mid = (a + b) / 2;
    if (labelY != null && b - a > 12 && mid - last > 20) { g.fillStyle = css("--ink2"); g.fillText(chrLabel(chr), mid, labelY); last = mid; }
  }
}

/** Runs of one integer state within a chromosome: [first bin, last bin + 1, state]. */
function runsOf(data, vals) {
  const out = [], nb = vals.length, chr = data.bins.chr;
  for (let s = 0; s < nb;) {
    let e = s + 1; while (e < nb && vals[e] === vals[s] && chr[e] === chr[s]) e++;
    if (vals[s] != null) out.push([s, e, vals[s]]);
    s = e;
  }
  return out;
}

const ProfilePlot = memo(forwardRef(function ProfilePlot({ data, cellName, height = 300, alleleMode = false, showCi = true, zoomX, theme }, ref) {
  const containerRef = useRef(); const canvasRef = useRef();
  const width = useWidthOf(containerRef);
  useImperativeHandle(ref, () => canvasRef.current, []);
  const nBins = data.bins.chr.length, genomeMax = data.bins.end_cum[nBins - 1];
  const [gx0, gx1] = zoomX || [0, genomeMax];
  const isAS = !!(alleleMode && data.nMajor?.[cellName] && data.nMinor?.[cellName]);
  const raw = cellName ? data.raw?.[cellName] : null;

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas || !cellName || !data.profiles[cellName]) return;
    const g = fit(canvas, width, height); g.clearRect(0, 0, width, height);
    const plotW = width - PR.l - PR.r, plotH = height - PR.t - PR.b;
    const xS = v => PR.l + ((v - gx0) / (gx1 - gx0)) * plotW;
    const prof = data.profiles[cellName];
    const top = (arr) => { let m = 0; for (const v of arr) if (v != null && v > m) m = v; return m; };
    const ymax = Math.max(5, Math.min(12, (isAS ? top(data.nMajor[cellName]) : top(prof)) + 1));
    const Y = v => PR.t + plotH - (Math.max(-0.3, Math.min(ymax + 0.4, v)) / (ymax + 0.4)) * plotH;

    // grid: one line per copy number, solid at 2
    g.font = `11px ${FONT}`; g.strokeStyle = css("--grid"); g.lineWidth = 1; g.fillStyle = css("--ink3"); g.textAlign = "right"; g.textBaseline = "middle";
    for (let y = 0; y <= ymax; y++) {
      g.setLineDash(y === 2 ? [] : [3, 4]); g.beginPath(); g.moveTo(PR.l, Y(y)); g.lineTo(PR.l + plotW, Y(y)); g.stroke();
      g.fillText(y, PR.l - 8, Y(y));
    }
    g.setLineDash([]);
    drawChromosomes(g, data, xS, gx0, gx1, PR.t, plotH, PR.t + plotH + 18);
    g.save(); g.beginPath(); g.rect(PR.l, 0, plotW, height); g.clip();

    let b0 = 0, b1 = nBins;
    for (let i = 0; i < nBins; i++) if (data.bins.end_cum[i] >= gx0) { b0 = i; break; }
    for (let i = nBins - 1; i >= 0; i--) if (data.bins.start_cum[i] <= gx1) { b1 = i + 1; break; }
    const nz = b1 - b0;
    // bins on the copy-number scale
    if (raw) {
      g.fillStyle = css("--points"); g.globalAlpha = 0.55;
      const r = nz > 1500 ? 1.3 : nz > 400 ? 1.8 : 2.6;
      for (let i = b0; i < b1; i++) {
        const v = raw[i]; if (v == null || !isFinite(v)) continue;
        g.beginPath(); g.arc(xS((data.bins.start_cum[i] + data.bins.end_cum[i]) / 2), Y(v), r, 0, 7); g.fill();
      }
      g.globalAlpha = 1;
    }
    const runs = runsOf(data, prof);
    // integer calls (total) — or the two alleles
    if (!isAS) {
      g.strokeStyle = css("--call"); g.lineWidth = 3; g.lineCap = "round"; g.beginPath();
      for (const [s, e, cn] of runs) { g.moveTo(xS(data.bins.start_cum[s]) + 1, Y(Math.min(cn, ymax))); g.lineTo(xS(data.bins.end_cum[e - 1]) - 1, Y(Math.min(cn, ymax))); }
      g.stroke();
    } else {
      const ci = showCi && data.ci?.[cellName];
      for (const [arr, col, lo, hi, off] of [[data.nMajor[cellName], "--major", "nMajor_lower", "nMajor_upper", 0.07], [data.nMinor[cellName], "--minor", "nMinor_lower", "nMinor_upper", -0.07]]) {
        if (ci && ci[lo]) {
          g.fillStyle = css(col); g.globalAlpha = 0.2;
          for (let i = b0; i < b1; i++) {
            if (arr[i] == null || ci[lo][i] == null) continue;
            const yTop = Y(Math.max(arr[i], ci[hi][i])), yBot = Y(Math.min(arr[i], ci[lo][i]));
            g.fillRect(xS(data.bins.start_cum[i]), yTop, xS(data.bins.end_cum[i]) - xS(data.bins.start_cum[i]), yBot - yTop);
          }
          g.globalAlpha = 1;
        }
        g.strokeStyle = css(col); g.lineWidth = 3; g.lineCap = "round"; g.beginPath();
        for (const [s, e, v] of runsOf(data, arr)) { g.moveTo(xS(data.bins.start_cum[s]) + 1, Y(Math.min(v, ymax) + off)); g.lineTo(xS(data.bins.end_cum[e - 1]) - 1, Y(Math.min(v, ymax) + off)); }
        g.stroke();
      }
    }
    // segment means of the bins
    if (raw) {
      g.strokeStyle = css("--mean"); g.lineWidth = 1.2; g.lineCap = "butt"; g.beginPath();
      for (const [s, e] of runs) {
        let m = 0, k = 0; for (let i = s; i < e; i++) if (raw[i] != null && isFinite(raw[i])) { m += raw[i]; k++; }
        if (!k) continue; m /= k;
        g.moveTo(xS(data.bins.start_cum[s]), Y(m)); g.lineTo(xS(data.bins.end_cum[e - 1]), Y(m));
      }
      g.stroke();
    }
    g.restore();
    if (!raw) { g.fillStyle = css("--ink3"); g.textAlign = "left"; g.textBaseline = "alphabetic"; g.fillText("no per-bin values in this export; showing the copy-number calls", PR.l + 8, PR.t + 14); }
  }, [data, cellName, width, height, isAS, raw, showCi, gx0, gx1, nBins, theme]);

  const move = (e) => {
    if (!cellName) return;
    const r = canvasRef.current.getBoundingClientRect(), x = e.clientX - r.left, plotW = width - PR.l - PR.r;
    if (x < PR.l || x > PR.l + plotW) { tip(e, null); return; }
    const b = binAt(data.bins, gx0 + ((x - PR.l) / plotW) * (gx1 - gx0));
    const line = isAS ? `major|minor <b>${data.nMajor[cellName][b] ?? "–"}|${data.nMinor[cellName][b] ?? "–"}</b>` : `copy number <b>${data.profiles[cellName][b] ?? "–"}</b>`;
    const rd = data.reads?.[cellName]?.[b];
    tip(e, `<span class="k">${chrLabel(data.bins.chr[b])}:${fmtBp(data.bins.start[b])}–${fmtBp(data.bins.end[b])}</span><br>${line}`
      + (raw && raw[b] != null ? `<br><span class="k">bin on the CN scale</span> ${fmt(raw[b])}` : "")
      + (rd != null ? `<br><span class="k">reads</span> ${rd}` : ""));
  };

  const items = [];
  if (raw) items.push(["dot", "--points", "bin, on the copy-number scale"], ["line", "--mean", "segment mean"]);
  if (isAS) items.push(["thick", "--major", "major allele"], ["thick", "--minor", "minor allele"]);
  else items.push(["thick", "--call", "integer copy number"]);
  return (
    <div ref={containerRef} style={{ width: "100%" }}>
      {!cellName && <div className="empty">Click a row of the heatmap, or pick a cell on the left.</div>}
      <canvas ref={canvasRef} style={{ display: cellName ? "block" : "none" }} onMouseMove={move} onMouseLeave={(e) => tip(e, null)} />
      {cellName && <div className="legend lines">
        {items.map(([k, c, t]) => <span className="sw" key={t}><i style={{ width: k === "dot" ? 8 : 20, height: k === "thick" ? 4 : k === "dot" ? 8 : 2, borderRadius: k === "dot" ? 4 : 2, background: `var(${c})` }} />{t}</span>)}
        {isAS && showCi && data.ci?.[cellName] && <span className="sw"><i style={{ width: 20, height: 8, borderRadius: 2, background: "var(--major)", opacity: 0.25 }} />BAF 5–95% interval</span>}
      </div>}
    </div>
  );
}));

/* ═══════════════════════════════════════════════════════════════════════════
   READ DEPTH TRACK (canvas) — per-bin read counts on the profile's x scale,
   so it lines up bin for bin with the copy-number profile drawn above it.
   ═══════════════════════════════════════════════════════════════════════════ */
const ReadDepthPlot = memo(function ReadDepthPlot({ data, cellName, height = 120, zoomX, theme }) {
  const containerRef = useRef(); const canvasRef = useRef();
  const width = useWidthOf(containerRef);
  const reads = cellName ? data.reads?.[cellName] : null;
  const nBins = data.bins.chr.length, genomeMax = data.bins.end_cum[nBins - 1];
  const [gx0, gx1] = zoomX || [0, genomeMax];
  const T = 18, B = 6;

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas || !reads) return;
    const g = fit(canvas, width, height); g.clearRect(0, 0, width, height);
    const plotW = width - PR.l - PR.r, plotH = height - T - B;
    const xS = v => PR.l + ((v - gx0) / (gx1 - gx0)) * plotW;
    let b0 = 0, b1 = nBins;
    for (let i = 0; i < nBins; i++) if (data.bins.end_cum[i] >= gx0) { b0 = i; break; }
    for (let i = nBins - 1; i >= 0; i--) if (data.bins.start_cum[i] <= gx1) { b1 = i + 1; break; }
    let max = 1, sum = 0, n = 0;
    for (let i = b0; i < b1; i++) { const v = reads[i]; if (v == null) continue; if (v > max) max = v; sum += v; n++; }
    const Y = v => T + plotH - (Math.min(v, max) / max) * plotH;
    g.font = `11px ${FONT}`; g.strokeStyle = css("--grid"); g.lineWidth = 1; g.fillStyle = css("--ink3"); g.textAlign = "right"; g.textBaseline = "middle";
    for (const t of niceTicks(0, max, 3)) { g.beginPath(); g.moveTo(PR.l, Y(t)); g.lineTo(PR.l + plotW, Y(t)); g.stroke(); g.fillText(t, PR.l - 8, Y(t)); }
    drawChromosomes(g, data, xS, gx0, gx1, T, plotH, null);
    g.save(); g.beginPath(); g.rect(PR.l, 0, plotW, height); g.clip();
    g.fillStyle = css("--series1"); g.globalAlpha = 0.75;
    for (let i = b0; i < b1; i++) {
      const v = reads[i]; if (v == null) continue;
      const x0 = xS(data.bins.start_cum[i]), x1 = xS(data.bins.end_cum[i]), y = Y(v);
      g.fillRect(x0, y, Math.max(0.6, x1 - x0), T + plotH - y);
    }
    g.globalAlpha = 1;
    if (n) { g.strokeStyle = css("--ink3"); g.setLineDash([4, 4]); g.beginPath(); g.moveTo(PR.l, Y(sum / n)); g.lineTo(PR.l + plotW, Y(sum / n)); g.stroke(); g.setLineDash([]); }
    g.restore();
    g.fillStyle = css("--ink2"); g.textAlign = "left"; g.textBaseline = "alphabetic";
    g.fillText("Reads per bin", PR.l, 12);
    if (n) { g.fillStyle = css("--ink3"); g.textAlign = "right"; g.fillText(`dashed: mean ${fmt(sum / n, 1)}`, PR.l + plotW, 12); }
  }, [data, reads, width, height, gx0, gx1, nBins, theme]);

  return <div ref={containerRef} style={{ width: "100%", marginTop: reads ? 10 : 0 }}>{reads && <canvas ref={canvasRef} />}</div>;
});

/* ═══════════════════════════════════════════════════════════════════════════
   LORENZ CURVE (canvas) — cumulative share of reads vs. cumulative share of
   bins for the selected cell, sorted from lowest- to highest-covered bin.
   The gap between this curve and the diagonal line of equality is exactly
   what the Gini coefficient (annotated on-chart) quantifies.
   ═══════════════════════════════════════════════════════════════════════════ */
const LorenzCurvePlot = memo(function LorenzCurvePlot({ data, cellName, gini, height = 260, theme }) {
  const containerRef = useRef(); const canvasRef = useRef();
  const width = useWidthOf(containerRef, 320);
  const reads = cellName ? data.reads?.[cellName] : null;

  const points = useMemo(() => {
    if (!reads) return null;
    const v = reads.filter(x => x != null && x >= 0).sort((a, b) => a - b);
    const n = v.length;
    if (n === 0) return null;
    const total = v.reduce((a, b) => a + b, 0);
    if (total === 0) return null;
    const pts = [[0, 0]];
    let cum = 0;
    for (let i = 0; i < n; i++) { cum += v[i]; pts.push([(i + 1) / n, cum / total]); }
    return pts;
  }, [reads]);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas || !points) return;
    const g = fit(canvas, width, height); g.clearRect(0, 0, width, height);
    const box = { l: 52, t: 14, w: width - 66, h: height - 58 };
    const ticks = [0, 0.25, 0.5, 0.75, 1], pct = v => Math.round(v * 100) + "%";
    const { tx, ty } = axes(g, box, ticks, ticks, { x0: 0, x1: 1, y0: 0, y1: 1, xlab: "share of bins, fewest reads first", ylab: "share of reads", xfmt: pct, yfmt: pct });
    g.strokeStyle = css("--ink3"); g.setLineDash([4, 4]); g.beginPath(); g.moveTo(tx(0), ty(0)); g.lineTo(tx(1), ty(1)); g.stroke(); g.setLineDash([]);
    g.fillStyle = css("--accent-soft"); g.beginPath(); g.moveTo(tx(0), ty(0));
    for (const [x, y] of points) g.lineTo(tx(x), ty(y));
    g.lineTo(tx(1), ty(1)); g.closePath(); g.fill();
    g.strokeStyle = css("--accent"); g.lineWidth = 2; g.lineJoin = "round"; g.beginPath();
    points.forEach(([x, y], i) => (i ? g.lineTo(tx(x), ty(y)) : g.moveTo(tx(x), ty(y))));
    g.stroke();
    if (gini != null) { g.fillStyle = css("--ink"); g.font = `600 12px ${FONT}`; g.textAlign = "left"; g.textBaseline = "top"; g.fillText(`Gini ${gini.toFixed(3)}`, box.l + 10, box.t + 6); }
  }, [points, width, height, gini, theme]);

  return (
    <div ref={containerRef} style={{ width: "100%" }}>
      {points ? <canvas ref={canvasRef} /> : <div className="empty">No per-bin read counts in this export — re-export with the current ascatsc_to_web.R to include them.</div>}
    </div>
  );
});

/* ═══════════════════════════════════════════════════════════════════════════
   SMALL UI COMPONENTS
   ═══════════════════════════════════════════════════════════════════════════ */
function ThemeIcon() {
  return <svg width="14" height="14" viewBox="0 0 16 16" fill="none" stroke="currentColor" strokeWidth="1.5"><circle cx="8" cy="8" r="6" /><path d="M8 2a6 6 0 0 0 0 12z" fill="currentColor" /></svg>;
}

function Switch({ on, onChange, label }) {
  return <button className="switch" role="switch" aria-checked={on} aria-label={label} onClick={() => onChange(!on)} />;
}

function DownloadMenu({ onSave }) {
  const [open, setOpen] = useState(false);
  const menuRef = useRef();
  useEffect(() => {
    if (!open) return;
    const handler = (e) => { if (!menuRef.current?.contains(e.target)) setOpen(false); };
    document.addEventListener("mousedown", handler);
    return () => document.removeEventListener("mousedown", handler);
  }, [open]);
  return (
    <div className="menu" ref={menuRef}>
      <button className="iconbtn" onClick={() => setOpen(o => !o)} title="Save the figure" aria-haspopup="menu" aria-expanded={open}>
        <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2"><path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4M7 10l5 5 5-5M12 15V3" /></svg>Save
      </button>
      {open && <div className="pop" role="menu">{["PNG", "PDF"].map(f => <button key={f} role="menuitem" onClick={() => { onSave(f.toLowerCase()); setOpen(false); }}>{f}</button>)}</div>}
    </div>
  );
}

/** A slider whose far end means "off": cells beyond the value are flagged. */
function ThresholdSlider({ t, value, range, count, passOnly, onChange }) {
  const [lo, hi] = range, span = hi - lo || 1;
  const off = value == null;
  const pos = off ? (t.dir === "max" ? 1 : 0) : Math.max(0, Math.min(1, (value - lo) / span));
  return (
    <div className="slider">
      <div className="top"><span>{t.label}</span><span>{off ? "off" : t.f(value)}</span></div>
      <input type="range" min="0" max="1" step="0.001" value={pos} aria-label={t.label} style={{ "--p": (pos * 100).toFixed(1) + "%" }}
        onChange={(e) => {
          const p = +e.target.value, isOff = t.dir === "max" ? p >= 0.999 : p <= 0.001, v = lo + p * span;
          onChange(isOff ? null : t.int ? Math.round(v) : v);
        }} />
      <small className={"hint" + (count ? " on" : "")}>{count ? `${count.toLocaleString()} cell${count > 1 ? "s" : ""} ${passOnly ? "hidden" : "flagged"} by this threshold` : t.hint}</small>
    </div>
  );
}

const TAB_ICONS = {
  cn: <svg width="14" height="14" viewBox="0 0 16 16" fill="none" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true"><path d="M1.5 9.5h3v-4h3v6h3v-8h2.5" /><path d="M1.5 14h13" opacity=".45" /></svg>,
  qc: <svg width="14" height="14" viewBox="0 0 16 16" fill="none" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true"><circle cx="4" cy="11" r="1.3" /><circle cx="7.5" cy="6.5" r="1.3" /><circle cx="11" cy="9" r="1.3" /><circle cx="13" cy="4" r="1.3" /><path d="M1.5 14.5h13M1.5 1.5v13" opacity=".45" /></svg>,
  data: <svg width="14" height="14" viewBox="0 0 16 16" fill="none" stroke="currentColor" strokeWidth="1.5" strokeLinecap="round" strokeLinejoin="round" aria-hidden="true"><rect x="2.5" y="1.5" width="11" height="13" rx="2" /><path d="M5.5 5.5h5M5.5 8.5h5M5.5 11.5h3" /></svg>,
};

// Sidebar QC thresholds: a cell beyond any active one is flagged
const THRESH = [
  { key: "residual", metric: "median_residual", label: "Max median residual", dir: "max", reason: "high residual", tile: "High residual", hint: "distance of the bins from the integer fit", f: v => v.toFixed(2) },
  { key: "mapd", metric: "mapd", label: "Max MAPD", dir: "max", reason: "noisy", tile: "Noisy", hint: "noise of the cell's logR", f: v => v.toFixed(3) },
  { key: "coverage", metric: "coverage", label: "Min coverage", dir: "min", reason: "low coverage", tile: "Low coverage", hint: "ASCAT.sc's coverage of the cell", f: v => Math.round(v) + "×", int: true },
  { key: "bins_with_cna", metric: "bins_with_cna", label: "Max bins with CNA", dir: "max", reason: "many CNA bins", tile: "Many CNA bins", hint: "bins whose copy number is not 2", f: v => String(Math.round(v)), int: true },
  { key: "gini", metric: "gini", label: "Max Gini", dir: "max", reason: "uneven coverage", tile: "Uneven coverage", hint: "unevenness of the reads across bins", f: v => v.toFixed(3) },
];

function setFavicon(theme) {
  const el = document.getElementById("favicon");
  if (el) el.href = theme === "dark" ? LOGO_DARK : LOGO_LIGHT;
}

/* ═══════════════════════════════════════════════════════════════════════════
   MAIN APP
   ═══════════════════════════════════════════════════════════════════════════ */
export default function App() {
  const [data, setData] = useState(null);
  const [fileName, setFileName] = useState("");
  const [theme, setThemeState] = useState(() => document.documentElement.dataset.theme || "light");
  const [selectedCell, setSelectedCell] = useState(null);
  const [tab, setTab] = useState("cn");
  const heatmapPanelRef = useRef();
  const profilePlotRef = useRef();
  const searchRef = useRef();
  const [search, setSearch] = useState("");
  const [thresholds, setThresholds] = useState({ residual: null, mapd: null, coverage: null, bins_with_cna: null, gini: null });
  const [passOnly, setPassOnly] = useState(true);
  const [minSegmentMb, setMinSegmentMb] = useState(0);
  const [heatmapH, setHeatmapH] = useState(560);
  const [alleleMode, setAlleleMode] = useState(false);
  const [showCi, setShowCi] = useState(true);
  const [zoom, setZoom] = useState(null);
  const [showDendro, setShowDendro] = useState(true);
  const [cellTypeFilter, setCellTypeFilter] = useState("All");
  const [panel, setPanel] = useState(true);

  const setTheme = useCallback((t) => {
    document.documentElement.dataset.theme = t;
    try { localStorage.setItem("ascatsc-theme", t); } catch (e) { /* private mode */ }
    setFavicon(t); setThemeState(t);
  }, []);
  const toggleTheme = useCallback(() => setTheme(document.documentElement.dataset.theme === "dark" ? "light" : "dark"), [setTheme]);
  useEffect(() => { setFavicon(theme); }, [theme]);

  const hasAS = !!(data?.nMajor && Object.keys(data.nMajor).length > 0);
  const hasCi = !!(data?.ci && Object.keys(data.ci).length > 0);
  const hasDendro = !!(
    (data?.dendrogram && Array.isArray(data.dendrogram.merge) && Array.isArray(data.dendrogram.height) && data.dendrogram.merge.length > 0 && data.dendrogram.height.length === data.dendrogram.merge.length) ||
    (data?.dendrogram_as && Array.isArray(data.dendrogram_as.merge) && Array.isArray(data.dendrogram_as.height) && data.dendrogram_as.merge.length > 0 && data.dendrogram_as.height.length === data.dendrogram_as.merge.length)
  );
  const hasCellTypes = !!(data?.cell_types && Object.keys(data.cell_types).length > 0);

  // All unique cell type labels
  const cellTypeOptions = useMemo(() => {
    if (!hasCellTypes) return [];
    return [...new Set(Object.values(data.cell_types))].filter(Boolean).sort();
  }, [data, hasCellTypes]);

  const nBins = data?.bins?.chr?.length ?? 0;

  // Bin size in MB derived from the first bin's coordinates
  const binSizeMb = useMemo(() => {
    if (data?.bins?.start?.[0] != null && data?.bins?.end?.[0] != null)
      return (data.bins.end[0] - data.bins.start[0]) / 1e6;
    return 1;
  }, [data]);

  const minSegmentBins = minSegmentMb > 0 ? Math.max(1, Math.ceil(minSegmentMb / binSizeMb)) : 0;
  // Defer the expensive per-cell computation so the slider thumb moves instantly
  const deferredMinSegmentBins = useDeferredValue(minSegmentBins);

  // Single-pass: compute filtered profiles, allele arrays, quality, and heatmapData together
  const { heatmapData, effectiveQuality } = useMemo(() => {
    if (!data?.profiles) return { heatmapData: data, effectiveQuality: data?.quality ?? {} };
    if (!deferredMinSegmentBins) return { heatmapData: data, effectiveQuality: data.quality ?? {} };

    const profiles = {}, major = data.nMajor ? {} : null, minor = data.nMinor ? {} : null, quality = {};
    for (const cell of Object.keys(data.profiles)) {
      const orig = data.profiles[cell];
      const filtered = filterSmallSegments(orig, deferredMinSegmentBins);
      profiles[cell] = filtered;
      if (major && data.nMajor[cell]) {
        major[cell] = data.nMajor[cell].map((v, i) =>
          Math.round(orig[i] ?? 2) !== 2 && Math.round(filtered[i]) === 2 ? 1 : v);
        minor[cell] = data.nMinor[cell].map((v, i) =>
          Math.round(orig[i] ?? 2) !== 2 && Math.round(filtered[i]) === 2 ? 1 : v);
      }
      if (data.quality[cell]) {
        quality[cell] = { ...data.quality[cell], bins_with_cna: filtered.filter(v => Math.round(v) !== 2).length };
      }
    }
    return {
      heatmapData: { ...data, profiles, nMajor: major ?? data.nMajor, nMinor: minor ?? data.nMinor },
      effectiveQuality: quality,
    };
  }, [data, deferredMinSegmentBins]);

  // getQCs()-style metrics, computed once per dataset
  const qcm = useMemo(() => (data ? computeQCMetrics(data) : {}), [data]);
  const cells = useMemo(() => (data ? Object.keys(data.profiles) : []), [data]);
  const cellIndex = useMemo(() => new Map(cells.map((c, i) => [c, i])), [cells]);
  // every metric as an array aligned with `cells`; exported values win over computed ones
  const vals = useMemo(() => {
    const out = {};
    for (const [k] of METRIC_META) {
      const arr = cells.map(c => num(effectiveQuality[c]?.[k] ?? qcm[c]?.[k]));
      if (arr.some(v => !isNaN(v))) out[k] = arr;
    }
    return out;
  }, [cells, effectiveQuality, qcm]);
  const metAll = useMemo(() => metricSummaries(vals), [vals]);
  const ranges = useMemo(() => {
    const out = {};
    for (const t of THRESH) {
      const v = (vals[t.metric] || []).filter(x => !isNaN(x));
      if (v.length && Math.min(...v) < Math.max(...v)) out[t.key] = [Math.min(...v), Math.max(...v)];
    }
    return out;
  }, [vals]);

  // Default thresholds for a new dataset: median + 2 MAD for the residual and
  // MAPD, the 10th percentile of coverage, the 90th of CNA bins and Gini.
  useEffect(() => {
    if (!data) return;
    const all = Object.keys(data.profiles);
    const col = k => all.map(c => num(data.quality?.[c]?.[k] ?? qcm[c]?.[k])).filter(v => !isNaN(v)).sort((a, b) => a - b);
    const robust = k => { const v = col(k); if (!v.length) return null; const m = quantile(v, 0.5); const mad = 1.4826 * quantile(v.map(x => Math.abs(x - m)).sort((a, b) => a - b), 0.5); return mad > 0 ? Math.round((m + 2 * mad) * 1000) / 1000 : null; };
    const q = (k, p, f) => { const v = col(k); return v.length > 1 && v[0] < v[v.length - 1] ? f(quantile(v, p)) : null; };
    setThresholds({
      residual: robust("median_residual"),
      mapd: robust("mapd"),
      coverage: q("coverage", 0.1, Math.floor),
      bins_with_cna: q("bins_with_cna", 0.9, Math.ceil),
      gini: q("gini", 0.9, v => Math.round(v * 1000) / 1000),
    });
    setAlleleMode(false); setZoom(null); setShowDendro(true); setCellTypeFilter("All"); setMinSegmentMb(0); setPassOnly(true); setSearch("");
    const first = data.clustering_order?.[0] ?? all[0];
    setSelectedCell(first ?? null);
  }, [data, qcm]);

  // QC flags of every cell under the current thresholds
  const { passing, reasons, failCounts } = useMemo(() => {
    const p = new Array(cells.length).fill(true), r = new Array(cells.length).fill(""), counts = {};
    for (const t of THRESH) {
      const th = thresholds[t.key], v = vals[t.metric]; counts[t.key] = 0;
      if (th == null || !v || !ranges[t.key]) continue;
      for (let i = 0; i < cells.length; i++) {
        if (isNaN(v[i]) || !(t.dir === "max" ? v[i] > th : v[i] < th)) continue;
        p[i] = false; r[i] = r[i] ? r[i] + ", " + t.reason : t.reason; counts[t.key]++;
      }
    }
    return { passing: p, reasons: r, failCounts: counts };
  }, [cells, vals, thresholds, ranges]);
  const nPass = passing.filter(Boolean).length;

  const filteredCells = useMemo(() => {
    if (!data) return [];
    // In allele-specific mode use the AS-derived order (rbind maj/min clustering),
    // matching R's plot_allele_heatmap which recomputes hclust on the stacked matrix.
    const baseOrder = (alleleMode && data.clustering_order_as)
      ? data.clustering_order_as
      : (data.clustering_order || Object.keys(data.profiles));
    return baseOrder.filter(c => {
      const i = cellIndex.get(c);
      if (passOnly && i != null && !passing[i]) return false;
      if (cellTypeFilter !== "All" && data.cell_types?.[c] !== cellTypeFilter) return false;
      return true;
    });
  }, [data, alleleMode, passOnly, passing, cellIndex, cellTypeFilter]);
  const shownSet = useMemo(() => new Set(filteredCells), [filteredCells]);
  const listCells = useMemo(() => {
    const q = search.trim().toLowerCase();
    return q ? filteredCells.filter(c => c.toLowerCase().includes(q)) : filteredCells;
  }, [filteredCells, search]);

  useEffect(() => { setZoom(null); }, [filteredCells.length]);

  const navigateCell = useCallback(dir => {
    if (!filteredCells.length) return;
    const idx = filteredCells.indexOf(selectedCell);
    setSelectedCell(filteredCells[Math.max(0, Math.min(filteredCells.length - 1, idx < 0 ? 0 : idx + dir))]);
  }, [filteredCells, selectedCell]);

  useEffect(() => {
    const h = e => {
      if (e.target.tagName === "INPUT" || e.target.tagName === "SELECT") { if (e.key === "Escape") e.target.blur(); return; }
      if (e.metaKey || e.ctrlKey || e.altKey) return;
      if (e.key === "ArrowUp"   || e.key === "k") { e.preventDefault(); navigateCell(-1); }
      else if (e.key === "ArrowDown" || e.key === "j") { e.preventDefault(); navigateCell(1); }
      else if (e.key === "Escape") setZoom(null);
      else if (e.key === "/" && data) { e.preventDefault(); setPanel(true); setTimeout(() => searchRef.current?.focus(), 0); }
      else if (e.key === "t") toggleTheme();
    };
    window.addEventListener("keydown", h); return () => window.removeEventListener("keydown", h);
  }, [navigateCell, toggleTheme, data]);

  useEffect(() => {
    document.querySelector(".cell.on")?.scrollIntoView({ block: "nearest" });
  }, [selectedCell]);

  const loadData = useCallback((d, name) => { setData(d); setFileName(name); setTab("cn"); }, []);
  if (!data) return <UploadScreen onLoad={loadData} onToggleTheme={toggleTheme} />;

  const totalCells = cells.length;
  const selIdx = selectedCell != null ? cellIndex.get(selectedCell) : undefined;
  const metricOf = key => metAll.find(m => m.key === key);
  const valOf = key => (selIdx != null && vals[key] ? vals[key][selIdx] : NaN);
  const failsThreshold = metric => THRESH.some(t => t.metric === metric && thresholds[t.key] != null && ranges[t.key] && !isNaN(valOf(metric)) && (t.dir === "max" ? valOf(metric) > thresholds[t.key] : valOf(metric) < thresholds[t.key]));
  const activeThresholds = Object.fromEntries(THRESH.filter(t => thresholds[t.key] != null && ranges[t.key]).map(t => [t.metric, thresholds[t.key]]));
  const isZoomed = zoom !== null;
  const zoomLabel = () => {
    if (!zoom?.x) return "";
    const a = binAt(data.bins, zoom.x[0]), b = binAt(data.bins, zoom.x[1]);
    return `${chrLabel(data.bins.chr[a])}:${fmtBp(data.bins.start[a])} – ${chrLabel(data.bins.chr[b])}:${fmtBp(data.bins.end[b])}`;
  };
  const reset = () => { setData(null); setSelectedCell(null); setZoom(null); };

  const statusEl = selIdx == null ? null : passing[selIdx]
    ? <span className="status" style={{ color: "var(--good)" }}><svg width="12" height="12" viewBox="0 0 12 12"><circle cx="6" cy="6" r="6" fill="currentColor" /><path d="M3.4 6.2 5.2 8l3.4-3.6" stroke="#fff" strokeWidth="1.5" fill="none" /></svg>passes QC</span>
    : <span className="status" style={{ color: "var(--serious)" }}><svg width="12" height="12" viewBox="0 0 12 12"><path d="M6 .8 11.4 10.6H.6z" fill="currentColor" /><path d="M6 4.3v3M6 8.8v.2" stroke="#fff" strokeWidth="1.4" /></svg>flagged{reasons[selIdx] ? ": " + reasons[selIdx] : ""}</span>;
  const pills = [];
  if (selIdx != null) {
    const add = (k, label, v, title) => { if (!isNaN(valOf(k))) pills.push([k, label, v, title]); };
    add("ploidy", "ploidy", fmt(valOf("ploidy")));
    add("n_reads_total", "reads", fmtCount(valOf("n_reads_total")), `${fmtInt(valOf("n_reads_total"))} reads in this cell's bins`);
    add("reads_per_bin", "reads / bin", fmt(valOf("reads_per_bin"), 1));
    add("coverage", "coverage", fmtMetric("coverage", valOf("coverage")));
    add("mapd", "MAPD", fmt(valOf("mapd"), 3));
    add("median_residual", "residual", fmt(valOf("median_residual"), 3));
    add("gini", "Gini", fmt(valOf("gini"), 3));
    add("bins_with_cna", "CNA bins", fmtInt(valOf("bins_with_cna")), `bins whose copy number is not 2, of ${nBins}`);
    add("breakpoints", "breakpoints", fmtInt(valOf("breakpoints")));
  }

  const profileCard = (
    <div className="card" key="profile">
      <div className="hd">
        <h2>{selectedCell ? <span className="mono">{selectedCell}</span> : "Cell profile"}</h2>
        <span className="sub">{statusEl}{selectedCell && data.cell_types?.[selectedCell] && <span className="ok"> · {data.cell_types[selectedCell]}</span>}{alleleMode && hasAS && <span className="ok"> · allele-specific</span>}</span>
        <DownloadMenu onSave={fmtName => saveCanvas([profilePlotRef.current], `${fileName || "ascatsc"}_${selectedCell || "profile"}`, fmtName)} />
      </div>
      {pills.length > 0 && <div className="pills">{pills.map(([k, label, v, title]) => (
        <div className={"pill" + (failsThreshold(k) ? " bad" : "")} key={k} title={title || metricOf(k)?.help}><span className="k">{label}</span><span className="v">{v}</span></div>
      ))}</div>}
      <div className="bd">
        <ProfilePlot ref={profilePlotRef} data={data} cellName={selectedCell} height={300} alleleMode={alleleMode} showCi={showCi} zoomX={zoom?.x} theme={theme} />
        <ReadDepthPlot data={data} cellName={selectedCell} zoomX={zoom?.x} theme={theme} />
      </div>
    </div>
  );

  const qcTiles = [["Cells passing", nPass.toLocaleString(), `of ${totalCells.toLocaleString()} (${totalCells ? Math.round((100 * nPass) / totalCells) : 0}%)`]];
  for (const t of THRESH) if (thresholds[t.key] != null && ranges[t.key]) qcTiles.push([t.tile, failCounts[t.key].toLocaleString(), `${t.dir === "max" ? "above" : "below"} ${t.f(thresholds[t.key])}`]);

  const readsTotal = vals.n_reads_total ? vals.n_reads_total.filter(v => !isNaN(v)) : null;
  const datasetTiles = [["Cells", totalCells.toLocaleString(), hasCellTypes ? `${cellTypeOptions.length} cell types` : ""]];
  if (readsTotal?.length) datasetTiles.push(["Reads", fmtCount(readsTotal.reduce((a, v) => a + v, 0)), `median ${fmtCount(quantile([...readsTotal].sort((a, b) => a - b), 0.5))} per cell`]);
  datasetTiles.push(["Bins", nBins.toLocaleString(), binSizeMb ? fmtBp(binSizeMb * 1e6) + " each" : ""]);
  datasetTiles.push(["Chromosomes", data.chr_info.length, fmtBp(data.bins.end_cum[nBins - 1]) + " in bins"]);
  datasetTiles.push(["Allele-specific", hasAS ? "yes" : "no", hasAS ? (hasCi ? "with BAF intervals" : "major and minor copies") : "total copy number only"]);
  datasetTiles.push(["Clustering", hasDendro ? "Ward.D2" : "–", hasDendro ? "Manhattan distance of profiles" : "no dendrogram"]);

  return (
    <div className={"app" + (panel ? "" : " nopanel")}>
      <header>
        <div className="brand"><img className="logo logo-light" src={LOGO_LIGHT} alt="ASCAT.sc" /><img className="logo logo-dark" src={LOGO_DARK} alt="ASCAT.sc" />ASCAT.sc</div>
        <div className="project" title={fileName}>{fileName}</div>
        <div className="chips">
          <span className="chip">{totalCells.toLocaleString()} cells · {nBins.toLocaleString()} bins</span>
          {hasAS && <span className="chip as">allele-specific</span>}
        </div>
        <div className="spacer" />
        <button className="iconbtn" onClick={() => setPanel(v => !v)} title={panel ? "Hide the side panel" : "Show the side panel"} aria-pressed={panel}>
          <svg width="14" height="14" viewBox="0 0 16 16" fill="none" stroke="currentColor" strokeWidth="1.5"><rect x="1.5" y="2.5" width="13" height="11" rx="2" /><path d="M6 2.5v11" /></svg><span className="btnlabel">Panel</span>
        </button>
        <button className="iconbtn" onClick={reset} title="Open another file">
          <svg width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2"><path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4M7 10l5-5 5 5M12 5v12" /></svg><span className="btnlabel">Open</span>
        </button>
        <button className="iconbtn" onClick={toggleTheme} title="Light or dark (t)" aria-label="Switch between light and dark"><ThemeIcon /></button>
      </header>

      <aside>
        <div className="group">
          <h3>Copy number</h3>
          {hasAS && (
            <div className="seg">
              <button className={alleleMode ? "" : "on"} onClick={() => setAlleleMode(false)}>Total</button>
              <button className={alleleMode ? "on" : ""} onClick={() => setAlleleMode(true)}>Allele-specific</button>
            </div>
          )}
          {hasDendro && <div className="row"><div><label>Dendrogram</label><small>Ward clustering of profiles</small></div><Switch on={showDendro} onChange={setShowDendro} label="Dendrogram" /></div>}
          {alleleMode && hasCi && <div className="row"><div><label>BAF intervals</label><small>5–95% ribbons on the profile</small></div><Switch on={showCi} onChange={setShowCi} label="BAF intervals" /></div>}
          <div className="slider">
            <div className="top"><span>Min segment size</span><span>{minSegmentMb ? `${minSegmentMb} Mb` : "off"}</span></div>
            <input type="range" min="0" max={Math.max(1, Math.round(nBins * binSizeMb / 2))} step={Math.max(1, Math.ceil(binSizeMb))} value={minSegmentMb} aria-label="Minimum segment size"
              style={{ "--p": ((100 * minSegmentMb) / Math.max(1, Math.round(nBins * binSizeMb / 2))).toFixed(1) + "%" }}
              onChange={e => setMinSegmentMb(+e.target.value)} />
            <small className="hint">shorter non-diploid segments show as copy number 2</small>
          </div>
          <div className="slider">
            <div className="top"><span>Heatmap height</span><span>{heatmapH} px</span></div>
            <input type="range" min="250" max="900" step="10" value={heatmapH} aria-label="Heatmap height" style={{ "--p": ((100 * (heatmapH - 250)) / 650).toFixed(1) + "%" }} onChange={e => setHeatmapH(+e.target.value)} />
          </div>
        </div>

        <div className="group">
          <h3>Cells</h3>
          <div className="row"><div><label>Passing QC only</label><small>{nPass === totalCells ? "every cell passes the thresholds" : `hides ${(totalCells - nPass).toLocaleString()} flagged of ${totalCells.toLocaleString()} cells`}</small></div><Switch on={passOnly} onChange={setPassOnly} label="Passing QC only" /></div>
          {THRESH.filter(t => ranges[t.key]).map(t => (
            <ThresholdSlider key={t.key} t={t} value={thresholds[t.key]} range={ranges[t.key]} count={failCounts[t.key]} passOnly={passOnly}
              onChange={v => setThresholds(s => ({ ...s, [t.key]: v }))} />
          ))}
          {hasCellTypes && (
            <div className="row" style={{ display: "block" }}>
              <label>Cell type</label>
              <div className="seg wrap" style={{ marginTop: 8 }}>
                {["All", ...cellTypeOptions].map(ct => <button key={ct} className={cellTypeFilter === ct ? "on" : ""} onClick={() => setCellTypeFilter(ct)}>{ct}</button>)}
              </div>
            </div>
          )}
          <div className="row" style={{ background: "transparent", padding: "6px 2px 0" }}><span className="ok">{filteredCells.length.toLocaleString()} cells shown</span></div>
        </div>

        <div className="group">
          <h3>Find a cell</h3>
          <input ref={searchRef} className="search" value={search} onChange={e => setSearch(e.target.value)} placeholder="Barcode or name  ( / )" autoComplete="off" spellCheck="false" />
          <div className="cells" role="listbox" aria-label="Cells">
            {listCells.slice(0, 400).map(c => {
              const i = cellIndex.get(c), ok = i == null || passing[i];
              return (
                <button key={c} className={"cell" + (c === selectedCell ? " on" : "")} role="option" aria-selected={c === selectedCell} onClick={() => setSelectedCell(c)}>
                  <span className="dot" style={{ background: ok ? "var(--good)" : "var(--serious)" }} title={ok ? "passes QC" : "flagged: " + reasons[i]} />
                  <span className="nm">{c}</span>
                  {hasCellTypes && data.cell_types?.[c] && <span className="ct">{data.cell_types[c]}</span>}
                  <span className="pl">{i != null ? fmt(vals.ploidy?.[i]) : ""}</span>
                </button>
              );
            })}
            {listCells.length > 400 && <div className="more">and {(listCells.length - 400).toLocaleString()} more; refine the search</div>}
            {!listCells.length && <div className="more">No match.</div>}
          </div>
          <div className="navbtns"><button className="iconbtn" onClick={() => navigateCell(-1)}>↑ Previous</button><button className="iconbtn" onClick={() => navigateCell(1)}>Next ↓</button></div>
        </div>
        <div className="hint-keys"><kbd>↑</kbd> <kbd>↓</kbd> cells · <kbd>drag</kbd> zoom · <kbd>Esc</kbd> reset · <kbd>/</kbd> search · <kbd>t</kbd> theme</div>
      </aside>

      <main>
        <div className="tabs" role="tablist">
          {[["cn", "Copy number", "Heatmap and cell profiles"], ["qc", "Cell QC", "Per-cell quality metrics"], ["data", "Dataset", "What this export holds"]].map(([id, label, title]) => (
            <button key={id} className={"tab" + (tab === id ? " on" : "")} role="tab" aria-selected={tab === id} title={title} onClick={() => setTab(id)}>{TAB_ICONS[id]}{label}</button>
          ))}
        </div>

        {tab === "cn" && <section>
          <div className="card">
            <div className="hd">
              <h2>Genome-wide copy number</h2>
              <span className="sub">{filteredCells.length.toLocaleString()} of {totalCells.toLocaleString()} cells · {nBins.toLocaleString()} bins{isZoomed ? " (zoomed)" : ""}{alleleMode ? " · allele-specific" : ""}</span>
              {isZoomed && <button className="zoomchip" onClick={() => setZoom(null)}>{zoomLabel()}  ✕</button>}
              <DownloadMenu onSave={f => heatmapPanelRef.current?.download(`${fileName || "ascatsc"}_heatmap`, f)} />
            </div>
            <div className="bd">
              <HeatmapPanel ref={heatmapPanelRef} data={heatmapData} cellOrder={filteredCells} chrInfo={data.chr_info}
                selectedCell={selectedCell} onCellClick={setSelectedCell} height={heatmapH}
                alleleMode={alleleMode} zoom={zoom} onZoomChange={setZoom} showDendro={showDendro && hasDendro} theme={theme} />
              <div className="keyrow">
                <span className="keytitle">{alleleMode ? "Major | minor copies" : "Total copy number"}</span>
                {alleleMode
                  ? <div className="legend as" aria-label="allele-specific state colours">{Object.keys(AS_COLORS).filter(k => k !== "NA").map(k => <span className="sw" key={k}><i style={{ background: AS_COLORS[k] }} />{k}</span>)}</div>
                  : <div className="legend" aria-label="copy number colours">{[0, 1, 2, 3, 4, 5, 6].map(k => <span className="sw" key={k}><i style={{ background: `var(--cn${k})` }} />{k === 6 ? "6+" : k}</span>)}<span className="sw"><i style={{ background: "var(--cnna)" }} />NA</span></div>}
                <span className="ok" style={{ marginLeft: "auto" }}>drag to zoom · double-click to reset</span>
              </div>
            </div>
          </div>
          {profileCard}
          {selectedCell && <div className="grid2">
            <div className="card">
              <div className="hd"><h2>Coverage evenness</h2><span className="sub">Lorenz curve of the reads over bins</span></div>
              <div className="bd"><LorenzCurvePlot data={data} cellName={selectedCell} gini={valOf("gini")} theme={theme} /></div>
            </div>
            <div className="card">
              <div className="hd"><h2>This cell's QC</h2><span className="sub">against all cells · ▲ beyond median ± 3 MAD on the bad side</span></div>
              <div className="bd tablewrap"><table>
                <thead><tr><th>metric</th><th>cell</th><th>median</th><th /></tr></thead>
                <tbody>{metAll.map(m => {
                  const v = valOf(m.key), out = isOutlier(m, v), fail = failsThreshold(m.key);
                  return <tr key={m.key}><td title={m.help}>{m.label}</td><td style={{ fontWeight: 600, color: fail ? "var(--serious)" : undefined }}>{fmtMetric(m.key, v)}</td><td className="ok">{fmtMetric(m.key, m.median)}</td>
                    <td style={{ color: "var(--serious)", width: 1, whiteSpace: "nowrap" }}>{fail ? "threshold" : out ? "▲" : ""}</td></tr>;
                })}</tbody>
              </table></div>
            </div>
          </div>}
        </section>}

        {tab === "qc" && <section>
          <QCTab cells={cells} vals={vals} metAll={metAll} passing={passing} reasons={reasons} shown={shownSet} tiles={qcTiles} thresholds={activeThresholds}
            sel={selectedCell} onSelect={setSelectedCell} theme={theme}>{profileCard}</QCTab>
        </section>}

        {tab === "data" && <section>
          <div className="card"><div className="hd"><h2>Dataset</h2><span className="sub">{fileName}</span></div>
            <div className="bd" style={{ paddingBottom: 0 }}><div className="tiles">
              {datasetTiles.map(([k, v, s]) => <div className="tile" key={k}><div className="k">{k}</div><div className={"v" + (String(v).length > 9 ? " long" : "")} title={String(v)}>{v}</div><div className="s">{s}</div></div>)}
            </div></div>
          </div>
          <div className="card"><div className="hd"><h2>QC metrics</h2><span className="sub">{data.qc_metrics ? "from res$QC_metrics (getQCs), completed in the browser" : "computed in the browser from the export, as getQCs() does"}</span></div>
            <div className="bd"><span className="ok">{metAll.map(m => m.label).join(" · ")}</span></div>
          </div>
        </section>}
      </main>
    </div>
  );
}
