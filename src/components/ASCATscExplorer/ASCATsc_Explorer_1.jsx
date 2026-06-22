import { useState, useEffect, useRef, useMemo, useCallback, useDeferredValue, memo, forwardRef, useImperativeHandle } from "react";
import jsPDF from "jspdf";
import * as d3 from "d3";
import { cnColor, cnRGBDark, MAJOR_COLOR, MINOR_COLOR, CN_LEGEND_ENTRIES, AS_COLORS } from "../../constants/colors";
import { generateDemoData } from "./demoData";

/* Precomputed RGB lookup for allele-specific states (mirrors R's get_allele_specific_colors) */
function hexToRGBArr(h) {
  return [parseInt(h.slice(1,3),16), parseInt(h.slice(3,5),16), parseInt(h.slice(5,7),16)];
}
const AS_RGB = {};
for (const [state, hex] of Object.entries(AS_COLORS)) AS_RGB[state] = hexToRGBArr(hex);
const AS_RGB_NA = AS_RGB["NA"] || [190, 190, 190];

function asRGB(major, minor) {
  if (major == null || minor == null) return AS_RGB_NA;
  const maj = Math.max(major, minor), min = Math.min(major, minor);
  const total = maj + min;
  const key = total > 7 ? "7+" : `${maj}|${min}`;
  return AS_RGB[key] || AS_RGB_NA;
}

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
   UPLOAD SCREEN
   ═══════════════════════════════════════════════════════════════════════════ */
function UploadScreen({onLoad}){
  const[dragging,setDragging]=useState(false);const[error,setError]=useState(null);const[loading,setLoading]=useState(false);const fileRef=useRef();
  const handleFile=(file)=>{if(!file)return;if(file.name.endsWith(".rds")||file.name.endsWith(".RDS")){setError("This is an .rds file — convert to JSON first.\n\nIn R:\n  source(\"R/ascatsc_to_web.R\")\n  rds_to_web(\""+file.name+"\", \"ascat_data.json\")\n\nThen upload the resulting ascat_data.json here.\nSee the README at https://github.com/anagchar/ASCAT.sc-Explorer for details.");return;}setLoading(true);setError(null);const reader=new FileReader();reader.onload=e=>{try{const d=JSON.parse(e.target.result);if(!d.bins||!d.profiles)throw new Error("Invalid JSON: missing 'bins' or 'profiles'.");if(!d.chr_info)throw new Error("Missing 'chr_info'.");onLoad(normalizeLoadedData(d));}catch(err){setError(err.message);setLoading(false);}};reader.onerror=()=>{setError("Failed to read file");setLoading(false);};reader.readAsText(file);};
  return(
    <div className="min-h-screen flex items-center justify-center" style={{background:"linear-gradient(135deg,#1a1a1a 0%,#222222 50%,#1a1a1a 100%)"}}>
      <div className="w-full max-w-xl mx-4">
        <div className="text-center mb-10">
          <div className="flex items-center justify-center mb-3">
            <img src={process.env.PUBLIC_URL + "/ASCATsc_logo.svg"} alt="ASCAT.sc logo" className="h-16" />
          </div>
          <p className="text-gray-400 text-sm">Single-cell copy number visualization</p>
        </div>
        <div className={`border-2 border-dashed rounded-2xl p-12 text-center transition-all cursor-pointer ${dragging?"border-blue-400 bg-blue-500/10":"border-gray-600 hover:border-gray-500 bg-white/[0.02]"}`}
          onDragOver={e=>{e.preventDefault();setDragging(true);}} onDragLeave={()=>setDragging(false)}
          onDrop={e=>{e.preventDefault();setDragging(false);handleFile(e.dataTransfer.files[0]);}} onClick={()=>fileRef.current?.click()}>
          <input ref={fileRef} type="file" accept=".json,.gz,.rds" className="hidden" onChange={e=>handleFile(e.target.files[0])}/>
          {loading?(<div className="text-blue-400"><div className="inline-block w-8 h-8 border-2 border-blue-400 border-t-transparent rounded-full animate-spin mb-3"/><p>Loading...</p></div>):(
            <><svg className="mx-auto mb-4 text-gray-500" width="48" height="48" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="1.5"><path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4M7 10l5-5 5 5M12 5v12"/></svg>
              <p className="text-gray-300 text-lg mb-1">Drop <code className="text-blue-400 bg-blue-500/10 px-2 py-0.5 rounded text-sm">ascat_data.json</code> here</p></>)}
        </div>
        {error&&<div className="mt-4 p-4 bg-red-500/10 border border-red-500/30 rounded-xl text-red-400 text-sm whitespace-pre-wrap font-mono leading-relaxed">{error}</div>}
        <div className="mt-6 text-center">
          <button onClick={()=>onLoad(generateDemoData())} className="text-sm text-gray-500 hover:text-blue-400 transition-colors underline underline-offset-4 decoration-gray-700 hover:decoration-blue-400">Load demo dataset (200 cells, allele-specific)</button>
        </div>
        <div className="mt-8 rounded-xl p-4 text-left" style={{background:"rgba(255,255,255,0.03)",border:"1px solid #3a3a3a"}}>
          <div className="text-xs font-semibold text-gray-400 uppercase tracking-wider mb-3">How to export from R</div>
          <pre className="text-xs text-gray-500 font-mono leading-relaxed overflow-x-auto">{`source("R/ascatsc_to_web.R")
rds_to_web("your_results.rds", "ascat_data.json")`}</pre>
        </div>
      </div>
    </div>
  );
}

/* ═══════════════════════════════════════════════════════════════════════════
   DENDROGRAM CANVAS
   ═══════════════════════════════════════════════════════════════════════════ */
const DendrogramCanvas = memo(forwardRef(function DendrogramCanvas({ data, cellOrder, height, yRange, width = 100, lightMode, alleleMode = false }, ref) {
  const canvasRef = useRef();
  useImperativeHandle(ref, () => canvasRef.current, []);
  // In allele-specific mode, prefer the AS dendrogram so that LOH-driven
  // clustering differences are reflected in the dendrogram arms.
  const dendro = (alleleMode && data.dendrogram_as) ? data.dendrogram_as : data.dendrogram;
  const hasValidDendro = Array.isArray(dendro?.merge) && Array.isArray(dendro?.height) && dendro.merge.length > 0 && dendro.height.length === dendro.merge.length;
  const nCells = cellOrder.length;
  const MARGIN_TOP = 26;
  const MARGIN_BOT = 20;
  const plotH = height - MARGIN_TOP - MARGIN_BOT;
  const PAD = 8;

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
    // A node "spans" visible leaves if its count > 0.
    const nodeVisCount = new Int32Array(dendro.merge.length);
    for (let i = 0; i < dendro.merge.length; i++) {
      const [l, r] = dendro.merge[i];
      const lc = l < 0 ? (visibleLeafSet.has(-l - 1) ? 1 : 0) : nodeVisCount[l - 1];
      const rc = r < 0 ? (visibleLeafSet.has(-r - 1) ? 1 : 0) : nodeVisCount[r - 1];
      nodeVisCount[i] = lc + rc;
    }

    // Find the subtree root: the lowest merge node whose visible count equals
    // visibleLeafSet.size (i.e. the MRCA of all visible leaves).
    let subtreeRoot = dendro.merge.length - 1;
    for (let i = 0; i < dendro.merge.length; i++) {
      if (nodeVisCount[i] === visibleLeafSet.size) { subtreeRoot = i; break; }
    }

    // Collect only the merge rows that belong to the subtree rooted at subtreeRoot.
    // Walk the subtree and gather 1-based step indices.
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

    // Build a filtered merge/height containing only subtree steps,
    // remapping step codes to new 1-based indices.
    // Steps must be emitted in topological order (children before parents).
    const oldToNew = new Int32Array(dendro.merge.length).fill(-1);
    const newMerge = [];
    const newHeight = [];
    // Iterate original order (children always come before parents in hclust).
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
    const dpr = window.devicePixelRatio || 1;
    canvas.width = width * dpr; canvas.height = height * dpr;
    canvas.style.width = width + "px"; canvas.style.height = height + "px";
    const ctx = canvas.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, width, height);
    if (lightMode) { ctx.fillStyle = "#ffffff"; ctx.fillRect(0, 0, width, height); }

    const yr = yRange || [0, nCells];
    const ySpan = yr[1] - yr[0] || 1;
    const xScale = h => PAD + ((maxH - h) / maxH) * (width - PAD * 2);
    const yScale = ci => MARGIN_TOP + ((ci - yr[0]) / ySpan) * plotH;

    ctx.strokeStyle = lightMode ? "#000000" : "#64748b";
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
  }, [segments, maxH, hasValidDendro, nCells, width, height, yRange, plotH, lightMode]);

  if (!hasValidDendro || nCells < 2) return null;
  return <canvas ref={canvasRef} style={{ width, height, display: "block", flexShrink: 0 }} />;
}));

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
   HEATMAP (2-layer canvas: data layer + overlay layer)
   ═══════════════════════════════════════════════════════════════════════════ */
const HeatmapPanel = memo(forwardRef(function HeatmapPanel({ data, cellOrder, chrInfo, selectedCell, onCellClick, height = 500, alleleMode = false, zoom, onZoomChange, showDendro = false, lightMode = false }, ref) {
  const DENDRO_W = showDendro && data.dendrogram ? 100 : 0;
  const nCells = cellOrder.length;
  const yRange = useMemo(() => zoom?.y || [0, nCells], [zoom, nCells]);
  const cytobands = useCytobands();
  const dendroRef = useRef();
  const heatmapRef = useRef();

  useImperativeHandle(ref, () => ({
    download(basename = "heatmap", format = "png") {
      const hc = heatmapRef.current;
      const dc = dendroRef.current;
      if (!hc) return;
      const dW = dc ? dc.width : 0;
      const composite = document.createElement("canvas");
      composite.width = dW + hc.width;
      composite.height = hc.height;
      const ctx = composite.getContext("2d");
      ctx.fillStyle = lightMode ? "#ffffff" : "#1a1a1a";
      ctx.fillRect(0, 0, composite.width, composite.height);
      if (dc) ctx.drawImage(dc, 0, 0);
      ctx.drawImage(hc, dW, 0);
      if (format === "pdf") {
        const cssW = (dc ? dc.offsetWidth : 0) + hc.offsetWidth;
        const cssH = hc.offsetHeight;
        const dataUrl = composite.toDataURL("image/png");
        const pdf = new jsPDF({ orientation: cssW > cssH ? "landscape" : "portrait", unit: "px", format: [cssW, cssH], hotfixes: ["px_scaling"] });
        pdf.addImage(dataUrl, "PNG", 0, 0, cssW, cssH);
        pdf.save(basename + ".pdf");
      } else {
        const link = document.createElement("a");
        link.download = basename + ".png";
        link.href = composite.toDataURL("image/png");
        link.click();
      }
    }
  }), [lightMode]);

  return (
    <div className="flex w-full">
      {DENDRO_W > 0 && (
        <DendrogramCanvas ref={dendroRef} data={data} cellOrder={cellOrder} height={height} yRange={yRange} width={DENDRO_W} lightMode={lightMode} alleleMode={alleleMode} />
      )}
      <div className="flex-1 min-w-0">
        <HeatmapCanvas ref={heatmapRef} data={data} cellOrder={cellOrder} chrInfo={chrInfo}
          selectedCell={selectedCell} onCellClick={onCellClick} height={height}
          alleleMode={alleleMode} zoom={zoom} onZoomChange={onZoomChange} lightMode={lightMode} cytobands={cytobands} />
      </div>
    </div>
  );
}));

const HeatmapCanvas = memo(forwardRef(function HeatmapCanvas({ data, cellOrder, chrInfo, selectedCell, onCellClick, height, alleleMode, zoom, onZoomChange, lightMode, cytobands }, ref) {
  const dataCanvasRef = useRef();
  const overlayCanvasRef = useRef();
  useImperativeHandle(ref, () => dataCanvasRef.current, []);
  const containerRef = useRef();
  const [cw, setCw] = useState(800);
  const brushRef = useRef(null);
  const isDragging = useRef(false);
  const prevOverlaySize = useRef({ w: 0, h: 0 });

  const MARGIN = { top: 4, right: 75, bottom: 20, left: 0 };
  const IDEOGRAM_H = 16;
  const MARGIN_TOP_TOTAL = MARGIN.top + IDEOGRAM_H + 2;

  useEffect(() => {
    const ro = new ResizeObserver(entries => {
      const w = entries[0]?.contentRect.width;
      if (w && w > 0) setCw(Math.round(w));
    });
    if (containerRef.current) ro.observe(containerRef.current);
    return () => ro.disconnect();
  }, []);

  const nBins = data.bins.chr.length;
  const nCells = cellOrder.length;
  const plotW = cw - MARGIN.left - MARGIN.right;
  const plotH = height - MARGIN_TOP_TOTAL - MARGIN.bottom;
  const genomeMax = useMemo(() => nBins > 0 ? data.bins.end_cum[nBins - 1] : 1, [data, nBins]);

  const xRange = useMemo(() => zoom?.x || [0, genomeMax], [zoom, genomeMax]);
  const yRange = useMemo(() => zoom?.y || [0, nCells],    [zoom, nCells]);
  const xSpan = xRange[1] - xRange[0] || 1;
  const ySpan = yRange[1] - yRange[0] || 1;

  const ML = MARGIN.left, MTT = MARGIN_TOP_TOTAL;
  const toPxX = useCallback(gx => ML  + ((gx - xRange[0]) / xSpan) * plotW, [ML,  xRange, xSpan, plotW]);
  const toPxY = useCallback(ci => MTT + ((ci - yRange[0]) / ySpan) * plotH, [MTT, yRange, ySpan, plotH]);
  const toGx  = useCallback(px => xRange[0] + ((px - ML)  / plotW) * xSpan, [ML,  xRange, xSpan, plotW]);
  const toCi  = useCallback(py => yRange[0] + ((py - MTT) / plotH) * ySpan, [MTT, yRange, ySpan, plotH]);

  // ── DATA LAYER ──
  // alleleMode IS in the dep array — the fix for lag is to not debounce/memo it
  // separately; instead we ensure the draw is cheap by using pre-computed RGB fns.
  useEffect(() => {
    const canvas = dataCanvasRef.current;
    if (!canvas || nCells === 0) return;
    const dpr = window.devicePixelRatio || 1;
    canvas.width = cw * dpr; canvas.height = height * dpr;
    canvas.style.width = cw + "px"; canvas.style.height = height + "px";
    const ctx = canvas.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

    const bg    = lightMode ? "#ffffff" : "#1a1a1a";
    const chrTx = lightMode ? "#1e293b" : "#a0a0a0";
    const chrBd = lightMode ? "rgba(0,0,0,0.12)" : "rgba(255,255,255,0.08)";
    const cnRGB = cnRGBDark;

    ctx.fillStyle = bg; ctx.fillRect(0, 0, cw, height);

    // ── Ideogram (cytobands) ──
    if (cytobands && cytobands.length > 0) {
      const ideoY = MARGIN.top;
      // Always draw on a white background for correct G-banding appearance
      ctx.fillStyle = "#ffffff";
      ctx.fillRect(ML, ideoY, plotW, IDEOGRAM_H);
      const chrOffsets = {};
      chrInfo.forEach(ci => { chrOffsets[ci.chr] = ci.start_cum; });
      for (const band of cytobands) {
        const offset = chrOffsets[band.chr];
        if (offset === undefined) continue;
        const cumS = offset + band.start;
        const cumE = offset + band.end;
        if (cumE <= xRange[0] || cumS >= xRange[1]) continue;
        const px0 = Math.max(ML, Math.floor(toPxX(Math.max(cumS, xRange[0]))));
        const px1 = Math.min(ML + plotW, Math.ceil(toPxX(Math.min(cumE, xRange[1]))));
        if (px1 <= px0) continue;
        ctx.fillStyle = STAIN_COLORS[band.stain] || "#AAAAAA";
        ctx.fillRect(px0, ideoY, px1 - px0, IDEOGRAM_H);
      }
      // Chromosome outline boxes
      ctx.strokeStyle = "rgba(0,0,0,0.35)";
      ctx.lineWidth = 0.5;
      chrInfo.forEach(({ start_cum, end_cum }) => {
        if (end_cum <= xRange[0] || start_cum >= xRange[1]) return;
        const px0 = Math.max(ML, Math.floor(toPxX(Math.max(start_cum, xRange[0]))));
        const px1 = Math.min(ML + plotW, Math.ceil(toPxX(Math.min(end_cum, xRange[1]))));
        if (px1 > px0 + 1) ctx.strokeRect(px0 + 0.5, ideoY + 0.5, px1 - px0 - 1, IDEOGRAM_H - 1);
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
            [r, g, b] = asRGB(data.nMajor[cell][bi], data.nMinor[cell][bi]);
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
            const srcStart = rowBase;
            const dstStart = ((startRow + dy) * tileW + startCol) * 4;
            buf.copyWithin(dstStart, srcStart, srcStart + colW * 4);
          }
        }
      }
      octx.putImageData(imgData, 0, 0);
      ctx.drawImage(offscreen, tileX0, tileY0);
    }

    // ── Chromosome boundaries + labels ──
    ctx.strokeStyle = chrBd; ctx.lineWidth = 1;
    chrInfo.forEach(({ start_cum }) => {
      if (start_cum < xRange[0] || start_cum > xRange[1]) return;
      const px = Math.round(toPxX(start_cum)) + 0.5;
      ctx.beginPath(); ctx.moveTo(px, MTT); ctx.lineTo(px, MTT + plotH); ctx.stroke();
    });
    ctx.fillStyle = chrTx; ctx.font = "11px system-ui"; ctx.textAlign = "center";
    chrInfo.forEach(({ chr, mid_cum }) => {
      if (mid_cum < xRange[0] || mid_cum > xRange[1]) return;
      ctx.fillText(chr.replace("chr", ""), toPxX(mid_cum), MTT + plotH + 14);
    });

    // ── Right-side CN Legend ──
    const cnEntries = CN_LEGEND_ENTRIES;
    const legX = ML + plotW + 6;
    const boxW = 22, boxH = 15, gap = 2;
    if (!alleleMode) {
      ctx.font = "bold 8px system-ui"; ctx.fillStyle = chrTx; ctx.textAlign = "left";
      ctx.fillText("Copy no.", legX, MTT + 9);
      let cy = MTT + 14;
      ctx.font = "bold 9px system-ui";
      for (const [cn, label] of cnEntries) {
        ctx.fillStyle = cnColor(cn, false);
        ctx.fillRect(legX, cy, boxW, boxH);
        ctx.fillStyle = cn === 2 ? "#374151" : "#e2e8f0";
        ctx.textAlign = "center";
        ctx.fillText(label, legX + boxW / 2, cy + 10);
        cy += boxH + gap;
      }
      ctx.textAlign = "start";
    } else {
      // Allele-specific legend: state swatches grouped by total CN, matching R's build_allele_legend()
      const asGroups = [
        ["0|0"],
        ["1|0"],
        ["1|1","2|0"],
        ["2|1","3|0"],
        ["2|2","3|1","4|0"],
        ["3|2","4|1","5|0"],
        ["3|3","4|2","5|1","6|0"],
        ["4|3","5|2","6|1","7|0"],
        ["7+"],
      ];
      const swH = 10, swW = boxW, swGap = 1;
      ctx.font = "bold 8px system-ui"; ctx.fillStyle = chrTx; ctx.textAlign = "left";
      ctx.fillText("State", legX, MTT + 9);
      let cy = MTT + 13;
      for (const group of asGroups) {
        for (const state of group) {
          const [r, g, b] = AS_RGB[state] || AS_RGB_NA;
          ctx.fillStyle = `rgb(${r},${g},${b})`;
          ctx.fillRect(legX, cy, swW, swH);
          ctx.fillStyle = chrTx; ctx.font = "7px system-ui"; ctx.textAlign = "left";
          ctx.fillText(state, legX + swW + 2, cy + swH - 2);
          cy += swH + swGap;
        }
        cy += 2; // extra gap between groups
      }
      ctx.textAlign = "start";
    }
  }, [data, cellOrder, cw, height, alleleMode, lightMode, xRange, yRange, nBins, nCells, plotW, plotH, chrInfo, toPxX, toPxY, MTT, ML, cytobands]);

  // ── OVERLAY LAYER ──
  const drawOverlay = useCallback(() => {
    const canvas = overlayCanvasRef.current;
    if (!canvas) return;
    const dpr = window.devicePixelRatio || 1;
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
          ctx.strokeStyle = "#f97316"; ctx.lineWidth = 2;
          ctx.strokeRect(ML + 1, y0, plotW - 2, y1 - y0);
        }
      }
    }
    const b = brushRef.current;
    if (b && isDragging.current) {
      const x0 = Math.min(b.sx, b.cx), y0 = Math.min(b.sy, b.cy);
      const w = Math.abs(b.cx - b.sx), h = Math.abs(b.cy - b.sy);
      ctx.fillStyle = "rgba(59,130,246,0.12)";
      ctx.strokeStyle = "rgba(59,130,246,0.6)"; ctx.lineWidth = 1.5;
      ctx.fillRect(x0, y0, w, h); ctx.strokeRect(x0, y0, w, h);
    }
  }, [cw, height, selectedCell, cellOrder, nCells, toPxY, plotW, plotH, MTT, ML]);

  useEffect(() => { drawOverlay(); }, [drawOverlay]);

  const getXY = e => {
    const rect = overlayCanvasRef.current.getBoundingClientRect();
    return { x: e.clientX - rect.left, y: e.clientY - rect.top };
  };
  const handleMouseDown = e => {
    if (e.button !== 0) return;
    const p = getXY(e);
    brushRef.current = { sx: p.x, sy: p.y, cx: p.x, cy: p.y };
    isDragging.current = true;
  };
  const handleMouseMove = e => {
    if (!isDragging.current) return;
    const p = getXY(e);
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
      const cy0 = toCi(Math.min(b.sy, p.y)), cy1 = toCi(Math.max(b.sy, p.y));
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

  const canvasBg = lightMode ? "#f8fafc" : "#1a1a1a";
  return (
    <div ref={containerRef} className="w-full relative" style={{ background: canvasBg }}>
      <canvas ref={dataCanvasRef} style={{ width: cw, height, display: "block", borderRadius: "8px" }} />
      <canvas ref={overlayCanvasRef}
        style={{ width: cw, height, position: "absolute", top: 0, left: 0, cursor: "crosshair", borderRadius: "8px" }}
        onMouseDown={handleMouseDown} onMouseMove={handleMouseMove} onMouseUp={handleMouseUp}
        onMouseLeave={() => { if (isDragging.current) { isDragging.current = false; brushRef.current = null; drawOverlay(); } }}
        onDoubleClick={() => onZoomChange(null)} />
    </div>
  );
}));

/* ═══════════════════════════════════════════════════════════════════════════
   PROFILE PLOT (canvas) — total CN & allele-specific
   ═══════════════════════════════════════════════════════════════════════════ */
const ProfilePlot = memo(forwardRef(function ProfilePlot({ data, cellName, showRaw = true, height = 280, alleleMode = false, lightMode = false, showCi = true }, ref) {
  const containerRef = useRef(); const canvasRef = useRef(); const [width, setWidth] = useState(800);
  useImperativeHandle(ref, () => ({
    download(basename = "cell_profile", format = "png") {
      const canvas = canvasRef.current;
      if (!canvas) return;
      const out = document.createElement("canvas");
      out.width = canvas.width; out.height = canvas.height;
      const ctx = out.getContext("2d");
      ctx.fillStyle = lightMode ? "#ffffff" : "#1a1a1a";
      ctx.fillRect(0, 0, out.width, out.height);
      ctx.drawImage(canvas, 0, 0);
      if (format === "pdf") {
        const cssW = canvas.offsetWidth;
        const cssH = canvas.offsetHeight;
        const dataUrl = out.toDataURL("image/png");
        const pdf = new jsPDF({ orientation: cssW > cssH ? "landscape" : "portrait", unit: "px", format: [cssW, cssH], hotfixes: ["px_scaling"] });
        pdf.addImage(dataUrl, "PNG", 0, 0, cssW, cssH);
        pdf.save(basename + ".pdf");
      } else {
        const link = document.createElement("a");
        link.download = basename + ".png";
        link.href = out.toDataURL("image/png");
        link.click();
      }
    }
  }), [lightMode]);
  const ML = 48, MT = 24, MR = 20, MB = 36;
  useEffect(() => {
    const ro = new ResizeObserver(e => { const w = e[0]?.contentRect.width; if (w > 0) setWidth(Math.round(w)); });
    if (containerRef.current) ro.observe(containerRef.current);
    return () => ro.disconnect();
  }, []);

  useEffect(() => {
    const canvas = canvasRef.current;
    if (!canvas || !cellName || !data.profiles[cellName]) return;
    const dpr = window.devicePixelRatio || 1;
    canvas.width = width * dpr; canvas.height = height * dpr;
    canvas.style.width = width + "px"; canvas.style.height = height + "px";
    const ctx = canvas.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

    const plotW = width - ML - MR, plotH = height - MT - MB;
    const nBins = data.bins.chr.length, genomeMax = data.bins.end_cum[nBins - 1];
    const xS = v => ML + (v / genomeMax) * plotW;
    const isAS = alleleMode && data.nMajor?.[cellName] && data.nMinor?.[cellName];
    const prof = data.profiles[cellName];
    let maxCN = 5;
    if (isAS) maxCN = Math.max(5, (d3.max(data.nMajor[cellName].filter(v => v != null)) || 4) + 1);
    else      maxCN = Math.max(5, (d3.max(prof.filter(v => v != null)) || 4) + 1);
    const yS = v => MT + plotH - (v / maxCN) * plotH;

    const textCol = lightMode ? "#374151" : "#94a3b8";
    const axisCol = lightMode ? "#9ca3af" : "#374151";
    const refCol  = lightMode ? "#9ca3af" : "#4b5563";
    const bgEven  = lightMode ? "rgba(0,0,0,0.03)" : "rgba(255,255,255,0.02)";
    const bgOdd   = lightMode ? "rgba(0,0,0,0.06)" : "rgba(255,255,255,0.05)";
    const rawCol  = lightMode ? "rgba(156,163,175,0.4)" : "rgba(209,213,219,0.4)";

    // Background
    ctx.fillStyle = lightMode ? "#ffffff" : "transparent";
    ctx.clearRect(0, 0, width, height);
    if (lightMode) ctx.fillRect(0, 0, width, height);

    // Alternating chr backgrounds
    data.chr_info.forEach(({ start_cum, end_cum }, i) => {
      ctx.fillStyle = i % 2 === 0 ? bgEven : bgOdd;
      ctx.fillRect(xS(start_cum), MT, xS(end_cum) - xS(start_cum), plotH);
    });

    // Reference lines (dashed)
    ctx.strokeStyle = refCol; ctx.lineWidth = 1; ctx.setLineDash([4, 3]);
    (isAS ? [1, 2] : [2]).forEach(ref => {
      const py = yS(ref);
      ctx.beginPath(); ctx.moveTo(ML, py); ctx.lineTo(ML + plotW, py); ctx.stroke();
    });
    ctx.setLineDash([]);

    // Raw data dots
    if (!isAS && showRaw && data.raw?.[cellName]) {
      const rawVals = data.raw[cellName];
      ctx.fillStyle = rawCol;
      for (let i = 0; i < nBins; i += 2) {
        if (rawVals[i] == null) continue;
        const cx = xS((data.bins.start_cum[i] + data.bins.end_cum[i]) / 2);
        const cy = yS(rawVals[i]);
        ctx.beginPath(); ctx.arc(cx, cy, 2, 0, 6.2832); ctx.fill();
      }
    }

    // Segment lines
    const drawSegs = (vals, color, lw, yOff = 0) => {
      ctx.strokeStyle = color; ctx.lineWidth = lw; ctx.lineCap = "round";
      let ss = 0;
      ctx.beginPath();
      for (let i = 1; i <= nBins; i++) {
        if (i === nBins || vals[i] !== vals[ss] || data.bins.chr[i] !== data.bins.chr[ss]) {
          if (vals[ss] != null) {
            const py = yS(vals[ss]) + yOff;
            ctx.moveTo(xS(data.bins.start_cum[ss]), py);
            ctx.lineTo(xS(data.bins.end_cum[i - 1]), py);
          }
          ss = i;
        }
      }
      ctx.stroke();
    };

    if (isAS) {
      // CI ribbons — drawn before segment lines, same logic as R's plot_allele_profile
      const ciCell = showCi && data.ci?.[cellName];
      if (ciCell) {
        const { nMajor_lower, nMajor_upper, nMinor_lower, nMinor_upper } = ciCell;
        const nMaj = data.nMajor[cellName];
        const nMin = data.nMinor[cellName];
        ctx.globalAlpha = 0.2;
        // Major ribbon (red)
        ctx.fillStyle = MAJOR_COLOR;
        for (let i = 0; i < nBins; i++) {
          if (nMaj[i] == null || nMajor_lower[i] == null) continue;
          const yLo = yS(Math.max(nMaj[i], nMajor_upper[i]));
          const yHi = yS(Math.min(nMaj[i], nMajor_lower[i]));
          const x0 = xS(data.bins.start_cum[i]);
          const x1 = xS(data.bins.end_cum[i]);
          ctx.fillRect(x0, yLo, x1 - x0, yHi - yLo);
        }
        // Minor ribbon (blue)
        ctx.fillStyle = MINOR_COLOR;
        for (let i = 0; i < nBins; i++) {
          if (nMin[i] == null || nMinor_lower[i] == null) continue;
          const yLo = yS(Math.max(nMin[i], nMinor_upper[i]));
          const yHi = yS(Math.min(nMin[i], nMinor_lower[i]));
          const x0 = xS(data.bins.start_cum[i]);
          const x1 = xS(data.bins.end_cum[i]);
          ctx.fillRect(x0, yLo, x1 - x0, yHi - yLo);
        }
        ctx.globalAlpha = 1;
      }
      drawSegs(data.nMajor[cellName], MAJOR_COLOR, 2.5, -1);
      drawSegs(data.nMinor[cellName], MINOR_COLOR, 2.5, 1);
    } else {
      // Group segment lines by color so we batch strokes per CN value
      const segsByColor = new Map();
      let ss = 0;
      for (let i = 1; i <= nBins; i++) {
        if (i === nBins || prof[i] !== prof[ss] || data.bins.chr[i] !== data.bins.chr[ss]) {
          if (prof[ss] != null) {
            const col = cnColor(prof[ss], lightMode);
            if (!segsByColor.has(col)) segsByColor.set(col, []);
            segsByColor.get(col).push([data.bins.start_cum[ss], data.bins.end_cum[i - 1], prof[ss]]);
          }
          ss = i;
        }
      }
      ctx.lineWidth = 2.5; ctx.lineCap = "round";
      for (const [col, segs] of segsByColor) {
        ctx.strokeStyle = col; ctx.beginPath();
        for (const [s, e, cn] of segs) {
          const py = yS(cn);
          ctx.moveTo(xS(s), py); ctx.lineTo(xS(e), py);
        }
        ctx.stroke();
      }
    }

    // Y axis
    ctx.strokeStyle = axisCol; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(ML, MT); ctx.lineTo(ML, MT + plotH); ctx.stroke();
    ctx.fillStyle = textCol; ctx.font = "10px system-ui"; ctx.textAlign = "right";
    for (let cn = 0; cn <= maxCN; cn++) {
      const py = yS(cn);
      ctx.beginPath(); ctx.moveTo(ML - 4, py); ctx.lineTo(ML, py); ctx.stroke();
      ctx.fillText(cn, ML - 6, py + 3);
    }
    // Y label
    ctx.save(); ctx.translate(12, MT + plotH / 2); ctx.rotate(-Math.PI / 2);
    ctx.textAlign = "center"; ctx.font = "11px system-ui";
    ctx.fillText(isAS ? "Allele Copy Number" : "Copy Number", 0, 0);
    ctx.restore();

    // Chr labels
    ctx.fillStyle = textCol; ctx.font = "10px system-ui"; ctx.textAlign = "center";
    data.chr_info.forEach(({ chr, mid_cum }) => {
      ctx.fillText(chr.replace("chr", ""), xS(mid_cum), MT + plotH + 22);
    });

    // AS legend
    if (isAS) {
      const lgX = ML + 8, lgY = 14;
      ctx.lineWidth = 2.5; ctx.lineCap = "round";
      ctx.strokeStyle = MAJOR_COLOR; ctx.beginPath(); ctx.moveTo(lgX, lgY); ctx.lineTo(lgX + 18, lgY); ctx.stroke();
      ctx.fillStyle = textCol; ctx.textAlign = "left"; ctx.font = "10px system-ui";
      ctx.fillText("nMajor", lgX + 22, lgY + 4);
      ctx.strokeStyle = MINOR_COLOR; ctx.beginPath(); ctx.moveTo(lgX + 80, lgY); ctx.lineTo(lgX + 98, lgY); ctx.stroke();
      ctx.fillText("nMinor", lgX + 102, lgY + 4);
    }
  }, [data, cellName, width, height, showRaw, alleleMode, lightMode, showCi]);

  if (!cellName) return <div className="flex items-center justify-center h-40 text-gray-500 text-sm">Select a cell to view its profile</div>;
  return (
    <div ref={containerRef} className="w-full">
      <canvas ref={canvasRef} style={{ display: "block" }} />
    </div>
  );
}));

/* ═══════════════════════════════════════════════════════════════════════════
   QUALITY SCATTER (D3/SVG)
   ═══════════════════════════════════════════════════════════════════════════ */
const QualityScatter = memo(function QualityScatter({ data, thresholds, selectedCell, onCellClick, height = 320, lightMode = false }) {
  const containerRef = useRef(); const canvasRef = useRef(); const [width, setWidth] = useState(500);
  const ML = 56, MT = 16, MR = 20, MB = 44;
  useEffect(() => {
    const ro = new ResizeObserver(e => { const w = e[0]?.contentRect.width; if (w > 0) setWidth(Math.round(w)); });
    if (containerRef.current) ro.observe(containerRef.current);
    return () => ro.disconnect();
  }, []);

  // Store cell positions for click detection
  const cellPositionsRef = useRef([]);

  useEffect(() => {
    const canvas = canvasRef.current; if (!canvas) return;
    const dpr = window.devicePixelRatio || 1;
    canvas.width = width * dpr; canvas.height = height * dpr;
    canvas.style.width = width + "px"; canvas.style.height = height + "px";
    const ctx = canvas.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

    const plotW = width - ML - MR, plotH = height - MT - MB;
    const cells = Object.keys(data.quality);
    const vals = cells.map(c => data.quality[c]).filter(Boolean);
    const xMax = d3.max(vals, d => d.mapd) * 1.1 || 1;
    const yMax = d3.max(vals, d => d.median_residual) * 1.1 || 1;
    const xS = v => ML + (v / xMax) * plotW;
    const yS = v => MT + plotH - (v / yMax) * plotH;

    const textCol = lightMode ? "#374151" : "#94a3b8";
    const axisCol = lightMode ? "#9ca3af" : "#374151";

    ctx.clearRect(0, 0, width, height);
    if (lightMode) { ctx.fillStyle = "#ffffff"; ctx.fillRect(0, 0, width, height); }

    // Pass zone
    ctx.fillStyle = "rgba(16,185,129,0.06)";
    ctx.fillRect(ML, yS(thresholds.residual), xS(thresholds.mapd) - ML, plotH - (yS(thresholds.residual) - MT));

    // Threshold lines
    ctx.strokeStyle = "rgba(239,68,68,0.7)"; ctx.lineWidth = 1; ctx.setLineDash([5, 4]);
    ctx.beginPath(); ctx.moveTo(xS(thresholds.mapd), MT); ctx.lineTo(xS(thresholds.mapd), MT + plotH); ctx.stroke();
    ctx.beginPath(); ctx.moveTo(ML, yS(thresholds.residual)); ctx.lineTo(ML + plotW, yS(thresholds.residual)); ctx.stroke();
    ctx.setLineDash([]);

    // Circles — batch by color
    const positions = [];
    const groups = { pass: [], fail: [], sel: [] };
    cells.forEach(c => {
      const q = data.quality[c]; if (!q || q.mapd == null || q.median_residual == null) return;
      const pB = thresholds.bins_with_cna == null || q.bins_with_cna == null || q.bins_with_cna <= thresholds.bins_with_cna;
      const pass = q.mapd <= thresholds.mapd && q.median_residual <= thresholds.residual && pB;
      const sel = c === selectedCell;
      const cx = xS(q.mapd), cy = yS(q.median_residual);
      positions.push({ c, cx, cy });
      (sel ? groups.sel : pass ? groups.pass : groups.fail).push({ cx, cy, c });
    });
    cellPositionsRef.current = positions;

    const drawGroup = (items, color, r, alpha, stroke) => {
      ctx.globalAlpha = alpha; ctx.fillStyle = color;
      ctx.beginPath();
      for (const { cx, cy } of items) { ctx.moveTo(cx + r, cy); ctx.arc(cx, cy, r, 0, 6.2832); }
      ctx.fill();
      if (stroke) {
        ctx.globalAlpha = 1; ctx.strokeStyle = stroke; ctx.lineWidth = 2;
        ctx.beginPath();
        for (const { cx, cy } of items) { ctx.moveTo(cx + r, cy); ctx.arc(cx, cy, r, 0, 6.2832); }
        ctx.stroke();
      }
      ctx.globalAlpha = 1;
    };
    drawGroup(groups.fail, "#ef4444", 2.5, 0.5, null);
    drawGroup(groups.pass, "#3b82f6", 2.5, 0.5, null);
    drawGroup(groups.sel,  "#f97316", 5,   1,   "#ffffff");

    // Axes
    ctx.strokeStyle = axisCol; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(ML, MT); ctx.lineTo(ML, MT + plotH); ctx.lineTo(ML + plotW, MT + plotH); ctx.stroke();
    ctx.fillStyle = textCol; ctx.font = "10px system-ui";
    // X ticks
    ctx.textAlign = "center";
    const xTicks = d3.ticks(0, xMax, 6);
    xTicks.forEach(t => {
      const px = xS(t);
      ctx.beginPath(); ctx.moveTo(px, MT + plotH); ctx.lineTo(px, MT + plotH + 4); ctx.stroke();
      ctx.fillText(t.toFixed(2), px, MT + plotH + 14);
    });
    // Y ticks
    ctx.textAlign = "right";
    const yTicks = d3.ticks(0, yMax, 6);
    yTicks.forEach(t => {
      const py = yS(t);
      ctx.beginPath(); ctx.moveTo(ML - 4, py); ctx.lineTo(ML, py); ctx.stroke();
      ctx.fillText(t.toFixed(2), ML - 6, py + 3);
    });
    // Axis labels
    ctx.fillStyle = textCol; ctx.font = "11px system-ui"; ctx.textAlign = "center";
    ctx.fillText("MAPD", ML + plotW / 2, MT + plotH + 34);
    ctx.save(); ctx.translate(12, MT + plotH / 2); ctx.rotate(-Math.PI / 2);
    ctx.fillText("Median Residual", 0, 0); ctx.restore();
  }, [data, thresholds, selectedCell, width, height, lightMode]);

  const handleClick = useCallback(e => {
    const rect = canvasRef.current.getBoundingClientRect();
    const px = e.clientX - rect.left, py = e.clientY - rect.top;
    for (const { c, cx, cy } of cellPositionsRef.current) {
      if (Math.hypot(px - cx, py - cy) <= 6) { onCellClick(c); return; }
    }
  }, [onCellClick]);

  return (
    <div ref={containerRef} className="w-full">
      <canvas ref={canvasRef} style={{ display: "block", cursor: "pointer" }} onClick={handleClick} />
    </div>
  );
});

/* ═══════════════════════════════════════════════════════════════════════════
   SMALL UI COMPONENTS
   ═══════════════════════════════════════════════════════════════════════════ */
const CellMetrics = memo(function CellMetrics({ data, cellName, thresholds, lightMode, effectiveQuality }) {
  if (!cellName || !data.quality[cellName]) return null;
  const q = (effectiveQuality ?? data.quality)[cellName];
  const nBins = data.bins?.chr?.length ?? 0;
  const pR = q.median_residual <= thresholds.residual;
  const pM = q.mapd <= thresholds.mapd;
  const pC = thresholds.coverage == null || q.coverage == null || q.coverage >= thresholds.coverage;
  const pB = thresholds.bins_with_cna == null || q.bins_with_cna == null || q.bins_with_cna <= thresholds.bins_with_cna;
  const boxBg = lightMode ? "#f1f5f9" : "#1f1f1f";
  const accent = "#b5860d";
  const MetricBox = ({ label, val, pass, thresh, fmt = v => v?.toFixed(3) ?? "—", suffix = "" }) => (
    <div className="rounded-lg p-3 border" style={{ background: boxBg, borderColor: lightMode ? "#e2e8f0" : "#313131" }}>
      <div className="text-xs mb-1" style={{ color: lightMode ? "#6b7280" : accent }}>{label}</div>
      <div className="text-lg font-semibold font-mono" style={{ color: lightMode ? (pass ? "#34d399" : "#f87171") : accent }}>{fmt(val)}{suffix}</div>
      <div className="text-xs mt-1" style={{ color: lightMode ? (pass ? "#059669" : "#dc2626") : accent }}>
        {pass ? "PASS" : "FAIL"} ({thresh != null ? `≤ ${thresh.toFixed ? thresh.toFixed(2) : thresh}` : "—"})
      </div>
    </div>
  );
  return (
    <div className="grid grid-cols-2 gap-3 text-sm">
      <MetricBox label="Median Residual" val={q.median_residual} pass={pR} thresh={thresholds.residual} />
      <MetricBox label="MAPD" val={q.mapd} pass={pM} thresh={thresholds.mapd} />
      {q.coverage != null && (
        <MetricBox label="Coverage" val={q.coverage} pass={pC}
          thresh={thresholds.coverage} fmt={v => v != null ? Math.round(v).toString() : "—"} suffix="×" />
      )}
      {q.bins_with_cna != null && (
        <MetricBox label="Bins w/ CNA" val={q.bins_with_cna} pass={pB}
          thresh={thresholds.bins_with_cna}
          fmt={v => v != null ? `${v}${nBins ? ` / ${nBins}` : ""}` : "—"} />
      )}
    </div>
  );
});

/* ═══════════════════════════════════════════════════════════════════════════
   THEME TOGGLE BUTTON
   ═══════════════════════════════════════════════════════════════════════════ */
const ThemeToggle = memo(function ThemeToggle({ lightMode, onToggle }) {
  return (
    <button onClick={onToggle}
      title={lightMode ? "Switch to dark mode" : "Switch to light mode"}
      className="flex items-center gap-1.5 text-xs px-2.5 py-1.5 rounded-lg transition-colors"
      style={{ background: lightMode ? "#e2e8f0" : "#333333", color: lightMode ? "#374151" : "#a0a0a0" }}>
      {lightMode
        ? <svg width="13" height="13" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2"><circle cx="12" cy="12" r="5"/><path d="M12 1v2M12 21v2M4.22 4.22l1.42 1.42M18.36 18.36l1.42 1.42M1 12h2M21 12h2M4.22 19.78l1.42-1.42M18.36 5.64l1.42-1.42"/></svg>
        : <svg width="13" height="13" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2"><path d="M21 12.79A9 9 0 1111.21 3 7 7 0 0021 12.79z"/></svg>}
      {lightMode ? "Light" : "Dark"}
    </button>
  );
});

function DownloadMenu({ onDownload, style }) {
  const [open, setOpen] = useState(false);
  const menuRef = useRef();
  useEffect(() => {
    if (!open) return;
    const handler = (e) => { if (!menuRef.current?.contains(e.target)) setOpen(false); };
    document.addEventListener("mousedown", handler);
    return () => document.removeEventListener("mousedown", handler);
  }, [open]);
  return (
    <div ref={menuRef} className="relative">
      <button onClick={() => setOpen(o => !o)}
        className="flex items-center gap-1 text-xs px-2.5 py-1 rounded-md transition-colors"
        style={style}>
        <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2"><path d="M21 15v4a2 2 0 01-2 2H5a2 2 0 01-2-2v-4M7 10l5 5 5-5M12 15V3"/></svg>
        Download
        <svg width="10" height="10" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2"><path d="M6 9l6 6 6-6"/></svg>
      </button>
      {open && (
        <div className="absolute right-0 mt-1 rounded-md shadow-lg z-50 overflow-hidden" style={{ background: style.background, border: "1px solid rgba(255,255,255,0.1)", minWidth: "80px" }}>
          {["PNG", "PDF"].map(fmt => (
            <button key={fmt} onClick={() => { onDownload(fmt.toLowerCase()); setOpen(false); }}
              className="w-full text-left text-xs px-3 py-1.5 hover:bg-blue-600/20 transition-colors"
              style={{ color: style.color }}>
              {fmt}
            </button>
          ))}
        </div>
      )}
    </div>
  );
}

/* ═══════════════════════════════════════════════════════════════════════════
   MOBILE DETECTION HOOK
   ═══════════════════════════════════════════════════════════════════════════ */
function useIsMobile() {
  const [mobile, setMobile] = useState(() => window.innerWidth < 768);
  useEffect(() => {
    const mq = window.matchMedia("(max-width: 767px)");
    const handler = e => setMobile(e.matches);
    mq.addEventListener("change", handler);
    return () => mq.removeEventListener("change", handler);
  }, []);
  return mobile;
}

/* ═══════════════════════════════════════════════════════════════════════════
   MAIN APP
   ═══════════════════════════════════════════════════════════════════════════ */
export default function App() {
  const isMobile = useIsMobile();
  const [data, setData] = useState(null);
  const [selectedCell, setSelectedCell] = useState(null);
  const [tab, setTab] = useState("heatmap");
  const heatmapPanelRef = useRef();
  const profilePlotRef = useRef();
  const [search, setSearch] = useState("");
  const [thresholds, setThresholds] = useState({ residual: 1.5, mapd: 2.0, coverage: null, bins_with_cna: null });
  const [minSegmentMb, setMinSegmentMb] = useState(0);
  const [heatmapH, setHeatmapH] = useState(650);
  const [alleleMode, setAlleleMode] = useState(false);
  const [showCi, setShowCi] = useState(true);
  const [zoom, setZoom] = useState(null);
  const [showDendro, setShowDendro] = useState(true);
  const [lightMode, setLightMode] = useState(false);
  const [cellTypeFilter, setCellTypeFilter] = useState("All");
  const [sidebarVisible, setSidebarVisible] = useState(true);

  const hasAS = !!(data?.nMajor && Object.keys(data.nMajor).length > 0);
  const hasCi = !!(data?.ci && Object.keys(data.ci).length > 0);
  const hasDendro = !!(
    (data?.dendrogram && Array.isArray(data.dendrogram.merge) && Array.isArray(data.dendrogram.height) && data.dendrogram.merge.length > 0 && data.dendrogram.height.length === data.dendrogram.merge.length) ||
    (data?.dendrogram_as && Array.isArray(data.dendrogram_as.merge) && Array.isArray(data.dendrogram_as.height) && data.dendrogram_as.merge.length > 0 && data.dendrogram_as.height.length === data.dendrogram_as.merge.length)
  );
  const hasCellTypes = !!(data?.cell_types && Object.keys(data.cell_types).length > 0);
  const hasCoverage = !!(data && Object.values(data.quality).some(q => q?.coverage != null));
  const hasBinsWithCna = !!(data && Object.values(data.quality).some(q => q?.bins_with_cna != null));

  // All unique cell type labels in clustering order
  const cellTypeOptions = useMemo(() => {
    if (!hasCellTypes) return [];
    const seen = new Set();
    const order = data.clustering_order || Object.keys(data.profiles);
    const types = [];
    for (const c of order) {
      const ct = data.cell_types[c];
      if (ct && !seen.has(ct)) { seen.add(ct); types.push(ct); }
    }
    return types.sort();
  }, [data, hasCellTypes]);

  // Coverage stats for slider range
  const coverageStats = useMemo(() => {
    if (!hasCoverage) return null;
    const vals = Object.values(data.quality).map(q => q?.coverage).filter(v => v != null);
    return { min: Math.floor(d3.min(vals)), max: Math.ceil(d3.max(vals)) };
  }, [data, hasCoverage]);

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

  useEffect(() => {
    if (!data) return;
    const vals = Object.values(data.quality).filter(Boolean);
    const medR = d3.median(vals.map(v => v.median_residual)) || 0;
    const madR = d3.median(vals.map(v => Math.abs(v.median_residual - medR))) * 1.4826 || 0.5;
    const medM = d3.median(vals.map(v => v.mapd)) || 0;
    const madM = d3.median(vals.map(v => Math.abs(v.mapd - medM))) * 1.4826 || 0.5;
    // Default coverage threshold: bottom 10th percentile (keep high-coverage cells)
    let covThresh = null;
    if (hasCoverage) {
      const covVals = vals.map(v => v.coverage).filter(v => v != null).sort((a, b) => a - b);
      covThresh = Math.floor(d3.quantile(covVals, 0.10));
    }
    let cnaThresh = null;
    if (hasBinsWithCna) {
      const cnaVals = vals.map(v => v.bins_with_cna).filter(v => v != null).sort((a, b) => a - b);
      if (cnaVals.length) cnaThresh = Math.ceil(d3.quantile(cnaVals, 0.90));
    }
    setThresholds({
      residual: Math.round((medR + 2 * madR) * 100) / 100,
      mapd: Math.round((medM + 2 * madM) * 100) / 100,
      coverage: covThresh,
      bins_with_cna: cnaThresh,
    });
    setAlleleMode(false); setZoom(null); setShowDendro(true); setCellTypeFilter("All"); setMinSegmentMb(0);
  }, [data, hasCoverage, hasBinsWithCna]);

  const filteredCells = useMemo(() => {
    if (!data) return [];
    // In allele-specific mode use the AS-derived order (rbind maj/min clustering),
    // matching R's plot_allele_heatmap which recomputes hclust on the stacked matrix.
    const baseOrder = (alleleMode && data.clustering_order_as)
      ? data.clustering_order_as
      : (data.clustering_order || Object.keys(data.profiles));
    return baseOrder.filter(c => {
      const q = effectiveQuality[c]; if (!q) return true;
      if (q.median_residual != null && q.median_residual > thresholds.residual) return false;
      if (q.mapd != null && q.mapd > thresholds.mapd) return false;
      if (thresholds.coverage != null && q.coverage != null && q.coverage < thresholds.coverage) return false;
      if (thresholds.bins_with_cna != null && q.bins_with_cna != null && q.bins_with_cna > thresholds.bins_with_cna) return false;
      if (cellTypeFilter !== "All" && data.cell_types?.[c] !== cellTypeFilter) return false;
      if (search && !c.toLowerCase().includes(search.toLowerCase())) return false;
      return true;
    });
  }, [data, alleleMode, thresholds, search, cellTypeFilter, effectiveQuality]);

  useEffect(() => { setZoom(null); }, [filteredCells.length]);

  const totalCells = data ? ((alleleMode && data.clustering_order_as) ? data.clustering_order_as : (data.clustering_order || Object.keys(data.profiles))).length : 0;
  const passRate = totalCells > 0 ? Math.round((filteredCells.length / totalCells) * 100) : 0;

  const navigateCell = useCallback(dir => {
    if (!filteredCells.length) return;
    const idx = filteredCells.indexOf(selectedCell);
    setSelectedCell(filteredCells[Math.max(0, Math.min(filteredCells.length - 1, idx < 0 ? 0 : idx + dir))]);
  }, [filteredCells, selectedCell]);

  useEffect(() => {
    const h = e => {
      if (e.target.tagName === "INPUT") return;
      if (e.key === "ArrowUp"   || e.key === "k") { e.preventDefault(); navigateCell(-1); }
      if (e.key === "ArrowDown" || e.key === "j") { e.preventDefault(); navigateCell(1); }
      if (e.key === "Escape") setZoom(null);
    };
    window.addEventListener("keydown", h); return () => window.removeEventListener("keydown", h);
  }, [navigateCell]);

  if (!data) return <UploadScreen onLoad={setData} />;

  const TABS = [{ id: "heatmap", label: "Heatmap" }, { id: "profile", label: "Profile" }, { id: "quality", label: "QC" }];
  const isZoomed = zoom !== null;

  // Theme-aware style helpers
  const bg     = lightMode ? "#f8fafc" : "#1a1a1a";
  const bg2    = lightMode ? "#ffffff" : "#222222";
  const bgCard = lightMode ? "#ffffff" : "#2a2a2a";
  const bgSide = lightMode ? "#f1f5f9" : "#1a1a1a";
  const bgItem = lightMode ? "#e2e8f0" : "#333333";
  const border = lightMode ? "#e2e8f0" : "#3a3a3a";
  const text   = lightMode ? "#1e293b" : "#e2e8f0";
  const textMd = lightMode ? "#374151" : "#d1d5db";
  const textSm = lightMode ? "#6b7280" : "#a0a0a0";
  const textXs = lightMode ? "#9ca3af" : "#707070";

  /* ── MOBILE LAYOUT ─────────────────────────────────────────────────────── */
  if (isMobile) {
    return (
      <div style={{ minHeight: "100dvh", background: bg, color: text, fontFamily: "'Inter',system-ui,sans-serif", display: "flex", flexDirection: "column" }}>

        {/* Mobile header */}
        <header style={{ display: "flex", alignItems: "center", justifyContent: "space-between", padding: "10px 14px", borderBottom: `1px solid ${border}`, background: bg2, flexShrink: 0 }}>
          <div style={{ display: "flex", alignItems: "center", gap: 8 }}>
            <img src={process.env.PUBLIC_URL + "/ASCATsc_logo.svg"} alt="ASCAT.sc logo" style={{ width: 28, height: 28, borderRadius: 6 }} />
            <span style={{ fontWeight: 700, fontSize: 13, fontFamily: "'JetBrains Mono',monospace", color: text }}>ASCAT.sc</span>
            <span style={{ fontSize: 11, padding: "1px 7px", borderRadius: 999, background: "rgba(16,185,129,0.1)", color: "#34d399", border: "1px solid rgba(16,185,129,0.2)" }}>
              {totalCells} cells
            </span>
          </div>
          <div style={{ display: "flex", alignItems: "center", gap: 6 }}>
            <ThemeToggle lightMode={lightMode} onToggle={() => setLightMode(v => !v)} />
            <button onClick={() => setSidebarVisible(v => !v)}
              title="Settings"
              style={{ background: bgItem, border: "none", borderRadius: 8, padding: "6px 8px", cursor: "pointer", color: textSm, display: "flex", alignItems: "center", gap: 4 }}>
              <svg width="15" height="15" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                <circle cx="12" cy="12" r="3"/><path d="M19.4 15a1.65 1.65 0 00.33 1.82l.06.06a2 2 0 010 2.83 2 2 0 01-2.83 0l-.06-.06a1.65 1.65 0 00-1.82-.33 1.65 1.65 0 00-1 1.51V21a2 2 0 01-4 0v-.09A1.65 1.65 0 009 19.4a1.65 1.65 0 00-1.82.33l-.06.06a2 2 0 01-2.83-2.83l.06-.06A1.65 1.65 0 004.68 15a1.65 1.65 0 00-1.51-1H3a2 2 0 010-4h.09A1.65 1.65 0 004.6 9a1.65 1.65 0 00-.33-1.82l-.06-.06a2 2 0 012.83-2.83l.06.06A1.65 1.65 0 009 4.68a1.65 1.65 0 001-1.51V3a2 2 0 014 0v.09a1.65 1.65 0 001 1.51 1.65 1.65 0 001.82-.33l.06-.06a2 2 0 012.83 2.83l-.06.06A1.65 1.65 0 0019.4 9a1.65 1.65 0 001.51 1H21a2 2 0 010 4h-.09a1.65 1.65 0 00-1.51 1z"/>
              </svg>
            </button>
          </div>
        </header>

        {/* Settings drawer overlay */}
        {sidebarVisible && (
          <div style={{ position: "fixed", inset: 0, zIndex: 50, display: "flex", flexDirection: "column" }}>
            <div style={{ position: "absolute", inset: 0, background: "rgba(0,0,0,0.6)" }} onClick={() => setSidebarVisible(false)} />
            <div style={{ position: "relative", marginTop: "auto", background: bgSide, borderRadius: "18px 18px 0 0", maxHeight: "80dvh", overflowY: "auto", padding: 16, zIndex: 1 }}>
              <div style={{ width: 36, height: 4, borderRadius: 2, background: border, margin: "0 auto 16px" }} />
              <div style={{ display: "flex", justifyContent: "space-between", alignItems: "center", marginBottom: 14 }}>
                <span style={{ fontSize: 13, fontWeight: 600, color: text }}>Settings</span>
                <button onClick={() => setSidebarVisible(false)} style={{ background: "none", border: "none", cursor: "pointer", color: textSm, padding: 4 }}>
                  <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2"><path d="M18 6L6 18M6 6l12 12"/></svg>
                </button>
              </div>

              {/* Toggles */}
              {(hasAS || hasDendro) && (
                <div style={{ marginBottom: 16 }}>
                  <div style={{ fontSize: 10, fontWeight: 600, textTransform: "uppercase", letterSpacing: "0.08em", color: textSm, marginBottom: 8 }}>Display</div>
                  <div style={{ display: "flex", flexDirection: "column", gap: 8 }}>
                    {hasAS && (
                      <div style={{ display: "flex", alignItems: "center", justifyContent: "space-between", padding: "8px 10px", borderRadius: 8, background: bgItem }}>
                        <span style={{ fontSize: 13, color: textMd }}>Allele-specific</span>
                        <button onClick={() => setAlleleMode(!alleleMode)} style={{ position: "relative", width: 36, height: 20, borderRadius: 10, border: "none", cursor: "pointer", background: alleleMode ? "#b5860d" : "#4b5563", transition: "background 0.2s" }}>
                          <span style={{ position: "absolute", top: 2, left: alleleMode ? 18 : 2, width: 16, height: 16, borderRadius: "50%", background: "#fff", boxShadow: "0 1px 3px rgba(0,0,0,0.4)", transition: "left 0.2s ease-in-out" }} />
                        </button>
                      </div>
                    )}
                    {hasDendro && (
                      <div style={{ display: "flex", alignItems: "center", justifyContent: "space-between", padding: "8px 10px", borderRadius: 8, background: bgItem }}>
                        <span style={{ fontSize: 13, color: textMd }}>Dendrogram</span>
                        <button onClick={() => setShowDendro(!showDendro)} style={{ position: "relative", width: 36, height: 20, borderRadius: 10, border: "none", cursor: "pointer", background: showDendro ? "#b5860d" : "#4b5563", transition: "background 0.2s" }}>
                          <span style={{ position: "absolute", top: 2, left: showDendro ? 18 : 2, width: 16, height: 16, borderRadius: "50%", background: "#fff", boxShadow: "0 1px 3px rgba(0,0,0,0.4)", transition: "left 0.2s ease-in-out" }} />
                        </button>
                      </div>
                    )}
                  </div>
                </div>
              )}

              {/* Quality filters */}
              <div style={{ marginBottom: 16 }}>
                <div style={{ fontSize: 10, fontWeight: 600, textTransform: "uppercase", letterSpacing: "0.08em", color: textSm, marginBottom: 8 }}>Quality Filters</div>
                {[["Max Residual", "residual", 3], ["Max MAPD", "mapd", 4]].map(([label, key, max]) => (
                  <div key={key} style={{ marginBottom: 10 }}>
                    <div style={{ display: "flex", justifyContent: "space-between", fontSize: 12, marginBottom: 4 }}>
                      <span style={{ color: textSm }}>{label}</span>
                      <span style={{ fontFamily: "monospace", color: textMd }}>{thresholds[key].toFixed(2)}</span>
                    </div>
                    <input type="range" min="0" max={max} step="0.05" value={thresholds[key]}
                      onChange={e => setThresholds(t => ({ ...t, [key]: +e.target.value }))}
                      style={{ width: "100%", accentColor: "#3b82f6" }} />
                  </div>
                ))}
                <div style={{ fontSize: 12, color: textXs }}>{filteredCells.length} / {totalCells} cells ({passRate}%)</div>
              </div>

              {/* Cell search */}
              <div style={{ marginBottom: 16 }}>
                <div style={{ fontSize: 10, fontWeight: 600, textTransform: "uppercase", letterSpacing: "0.08em", color: textSm, marginBottom: 8 }}>Search Cell</div>
                <input type="text" placeholder="Barcode..." value={search} onChange={e => setSearch(e.target.value)}
                  style={{ width: "100%", padding: "8px 10px", borderRadius: 8, border: `1px solid ${border}`, background: bgItem, color: text, fontSize: 13, outline: "none", boxSizing: "border-box" }} />
              </div>

              {/* Cell list */}
              <div style={{ marginBottom: 16 }}>
                <div style={{ fontSize: 10, fontWeight: 600, textTransform: "uppercase", letterSpacing: "0.08em", color: textSm, marginBottom: 8 }}>Cell List</div>
                <div style={{ maxHeight: 160, overflowY: "auto", borderRadius: 8, border: `1px solid ${border}`, background: lightMode ? "#f8fafc" : "#222" }}>
                  {filteredCells.length === 0
                    ? <div style={{ padding: 12, textAlign: "center", fontSize: 12, color: textXs }}>No cells match</div>
                    : filteredCells.map(c => (
                      <button key={c} onClick={() => { setSelectedCell(c); setSidebarVisible(false); setTab("profile"); }}
                        style={{ width: "100%", textAlign: "left", padding: "8px 12px", fontSize: 12, fontFamily: "monospace", background: c === selectedCell ? "rgba(59,130,246,0.15)" : "transparent", color: c === selectedCell ? "#60a5fa" : textSm, border: "none", cursor: "pointer", display: "block" }}>
                        {c}
                      </button>
                    ))}
                </div>
              </div>

              {/* Nav buttons */}
              <div style={{ display: "flex", gap: 8, marginBottom: 16 }}>
                <button onClick={() => navigateCell(-1)} style={{ flex: 1, padding: "10px 0", borderRadius: 8, border: "none", background: lightMode ? "#e2e8f0" : "#b5860d", color: lightMode ? "#374151" : "#1a1a1a", fontWeight: 600, fontSize: 13, cursor: "pointer" }}>↑ Prev</button>
                <button onClick={() => navigateCell(1)} style={{ flex: 1, padding: "10px 0", borderRadius: 8, border: "none", background: lightMode ? "#e2e8f0" : "#b5860d", color: lightMode ? "#374151" : "#1a1a1a", fontWeight: 600, fontSize: 13, cursor: "pointer" }}>Next ↓</button>
              </div>

              <button onClick={() => { setData(null); setSelectedCell(null); setTab("heatmap"); setZoom(null); setSidebarVisible(false); }}
                style={{ width: "100%", padding: "10px 0", borderRadius: 8, border: `1px solid ${border}`, background: "transparent", color: textSm, fontSize: 13, cursor: "pointer" }}>
                Load new data
              </button>
            </div>
          </div>
        )}

        {/* Main content */}
        <div style={{ flex: 1, overflowY: "auto", padding: 12 }}>
          {tab === "heatmap" && (
            <div style={{ display: "flex", flexDirection: "column", gap: 12 }}>
              <div style={{ display: "flex", alignItems: "center", justifyContent: "space-between" }}>
                <span style={{ fontSize: 13, fontWeight: 600, color: textMd }}>{alleleMode ? "Allele-Specific Heatmap" : "CN Heatmap"}</span>
                <div style={{ display: "flex", alignItems: "center", gap: 8 }}>
                  {zoom && <button onClick={() => setZoom(null)} style={{ fontSize: 11, padding: "3px 8px", borderRadius: 6, background: "rgba(59,130,246,0.2)", color: "#93c5fd", border: "none", cursor: "pointer" }}>Reset zoom</button>}
                  <DownloadMenu onDownload={fmt => heatmapPanelRef.current?.download("heatmap", fmt)} style={{ background: bgItem, color: textSm }} />
                </div>
              </div>
              <div style={{ borderRadius: 10, overflow: "hidden", border: `1px solid ${border}`, background: bg2 }}>
                <HeatmapPanel ref={heatmapPanelRef} data={heatmapData} cellOrder={filteredCells} chrInfo={data.chr_info}
                  selectedCell={selectedCell} onCellClick={setSelectedCell} height={Math.max(280, Math.min(500, filteredCells.length * 2 + 60))}
                  alleleMode={alleleMode} zoom={zoom} onZoomChange={setZoom} showDendro={showDendro} lightMode={lightMode} />
              </div>
              {selectedCell && (
                <div style={{ borderRadius: 10, border: `1px solid ${border}`, padding: 12, background: bgCard }}>
                  <div style={{ fontSize: 11, fontWeight: 600, textTransform: "uppercase", letterSpacing: "0.07em", color: textSm, marginBottom: 6 }}>Profile — <span style={{ color: "#60a5fa", fontFamily: "monospace" }}>{selectedCell}</span></div>
                  <ProfilePlot data={data} cellName={selectedCell} height={180} alleleMode={alleleMode} lightMode={lightMode} showCi={showCi} />
                </div>
              )}
            </div>
          )}

          {tab === "profile" && (
            <div style={{ display: "flex", flexDirection: "column", gap: 12 }}>
              <div style={{ display: "flex", alignItems: "center", justifyContent: "space-between" }}>
                <span style={{ fontSize: 13, fontWeight: 600, color: textMd }}>Cell Profile</span>
                <DownloadMenu onDownload={fmt => profilePlotRef.current?.download(selectedCell || "cell_profile", fmt)} style={{ background: bgItem, color: textSm }} />
              </div>
              {selectedCell && <div style={{ fontSize: 12, fontFamily: "monospace", color: "#60a5fa", padding: "4px 10px", borderRadius: 6, background: bgItem, alignSelf: "flex-start" }}>{selectedCell}</div>}
              <div style={{ borderRadius: 10, border: `1px solid ${border}`, padding: 12, background: bgCard }}>
                <ProfilePlot ref={profilePlotRef} data={data} cellName={selectedCell} height={300} alleleMode={alleleMode} lightMode={lightMode} showCi={showCi} />
              </div>
              {selectedCell && (
                <div style={{ borderRadius: 10, border: `1px solid ${border}`, padding: 12, background: bgCard }}>
                  <div style={{ fontSize: 11, fontWeight: 600, textTransform: "uppercase", letterSpacing: "0.07em", color: textSm, marginBottom: 10 }}>Quality Metrics</div>
                  <CellMetrics data={data} cellName={selectedCell} thresholds={thresholds} lightMode={lightMode} effectiveQuality={effectiveQuality} />
                </div>
              )}
            </div>
          )}

          {tab === "quality" && (
            <div style={{ display: "flex", flexDirection: "column", gap: 12 }}>
              <span style={{ fontSize: 13, fontWeight: 600, color: textMd }}>Quality Overview</span>
              <div style={{ borderRadius: 10, border: `1px solid ${border}`, padding: 12, background: bgCard }}>
                <div style={{ fontSize: 11, color: textSm, marginBottom: 8 }}>MAPD vs Median Residual</div>
                <QualityScatter data={data} thresholds={thresholds} selectedCell={selectedCell} onCellClick={setSelectedCell} height={320} lightMode={lightMode} />
              </div>
              <div style={{ display: "grid", gridTemplateColumns: "1fr 1fr 1fr", gap: 8 }}>
                {[["Total", totalCells, "text-blue-400", "#60a5fa"], ["Pass", filteredCells.length, "text-emerald-400", "#34d399"], ["Rate", passRate + "%", "text-blue-400", "#60a5fa"]].map(([l, v, , c]) => (
                  <div key={l} style={{ borderRadius: 10, border: `1px solid ${border}`, padding: 10, background: bgCard }}>
                    <div style={{ fontSize: 10, color: textXs, marginBottom: 4 }}>{l}</div>
                    <div style={{ fontSize: 20, fontWeight: 700, fontFamily: "monospace", color: c }}>{v}</div>
                  </div>
                ))}
              </div>
            </div>
          )}
        </div>

        {/* Bottom tab bar */}
        <nav style={{ display: "flex", borderTop: `1px solid ${border}`, background: bg2, flexShrink: 0 }}>
          {TABS.map(t => (
            <button key={t.id} onClick={() => setTab(t.id)}
              style={{ flex: 1, padding: "12px 0 10px", fontSize: 12, fontWeight: 600, border: "none", background: "transparent", cursor: "pointer",
                color: tab === t.id ? "#3b82f6" : textSm,
                borderTop: tab === t.id ? "2px solid #3b82f6" : "2px solid transparent" }}>
              {t.label}
            </button>
          ))}
        </nav>
      </div>
    );
  }

  /* ── DESKTOP LAYOUT (unchanged) ─────────────────────────────────────────── */
  return (
    <div className="min-h-screen" style={{ background: bg, color: text, fontFamily: "'Inter',system-ui,sans-serif" }}>
      {/* Header */}
      <header className="flex items-center justify-between px-5 py-3 border-b" style={{ borderColor: border, background: bg2 }}>
        <div className="flex items-center gap-3 flex-wrap">
          <img src={process.env.PUBLIC_URL + "/ASCATsc_logo.svg"} alt="ASCAT.sc logo" className="w-8 h-8 rounded-lg object-contain" />
          <span className="font-bold tracking-tight text-sm" style={{ fontFamily: "'JetBrains Mono',monospace", color: text }}>ASCAT.sc Explorer</span>
          <span className="text-xs px-2 py-0.5 rounded-full bg-emerald-500/10 text-emerald-400 border border-emerald-500/20">
            {totalCells} cells &middot; {data.bins.chr.length} bins
          </span>
          {hasAS && <span className="text-xs px-2 py-0.5 rounded-full border" style={{ background: "rgba(181,134,13,0.15)", color: "#b5860d", borderColor: "rgba(181,134,13,0.3)" }}>Allele-specific</span>}
        </div>
        <div className="flex items-center gap-2">
          <button onClick={() => setSidebarVisible(v => !v)}
            title={sidebarVisible ? "Hide panel" : "Show panel"}
            className="flex items-center gap-1.5 text-xs px-2.5 py-1.5 rounded-lg transition-colors"
            style={{ background: lightMode ? "#e2e8f0" : "#333333", color: lightMode ? "#374151" : "#a0a0a0" }}>
            <svg width="13" height="13" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
              <rect x="3" y="3" width="18" height="18" rx="2"/><path d="M9 3v18"/>
            </svg>
            {sidebarVisible ? "Hide panel" : "Show panel"}
          </button>
          <ThemeToggle lightMode={lightMode} onToggle={() => setLightMode(v => !v)} />
          <button onClick={() => { setData(null); setSelectedCell(null); setTab("heatmap"); setZoom(null); }}
            className="text-xs px-3 py-1.5 rounded-lg hover:bg-white/5 transition-colors" style={{ color: textSm }}>Load new data</button>
        </div>
      </header>

      <div className="flex" style={{ height: "calc(100vh - 49px)" }}>
        {/* Sidebar */}
        {sidebarVisible && <aside className="flex-shrink-0 overflow-y-auto border-r" style={{ width: 280, borderColor: border, background: bgSide }}>
          <div className="p-4 space-y-4">
            {/* Tabs */}
            <div className="flex gap-1 p-1 rounded-lg" style={{ background: bgItem }}>
              {TABS.map(t => (
                <button key={t.id} onClick={() => setTab(t.id)}
                  className={`flex-1 text-xs py-2 rounded-md font-medium transition-all ${tab === t.id ? "bg-blue-600 text-white shadow" : ""}`}
                  style={tab !== t.id ? { color: textSm } : {}}>
                  {t.label}
                </button>
              ))}
            </div>

            {/* Toggles */}
            {(hasAS || hasDendro) && (
              <div className="space-y-2">
                {hasAS && (
                  <div className="flex items-center justify-between p-2.5 rounded-lg" style={{ background: bgItem }}>
                    <span className="text-xs" style={{ color: textMd }}>Allele-specific</span>
                    <button onClick={() => setAlleleMode(!alleleMode)} style={{ position:"relative", width:36, height:20, borderRadius:10, border:"none", cursor:"pointer", background: alleleMode ? "#b5860d" : "#4b5563", transition:"background 0.2s" }}>
                      <span style={{ position:"absolute", top:2, left: alleleMode ? 18 : 2, width:16, height:16, borderRadius:"50%", background:"#fff", boxShadow:"0 1px 3px rgba(0,0,0,0.4)", transition:"left 0.2s ease-in-out" }} />
                    </button>
                  </div>
                )}
                {hasDendro && (
                  <div className="flex items-center justify-between p-2.5 rounded-lg" style={{ background: bgItem }}>
                    <span className="text-xs" style={{ color: textMd }}>Dendrogram</span>
                    <button onClick={() => setShowDendro(!showDendro)} style={{ position:"relative", width:36, height:20, borderRadius:10, border:"none", cursor:"pointer", background: showDendro ? "#b5860d" : "#4b5563", transition:"background 0.2s" }}>
                      <span style={{ position:"absolute", top:2, left: showDendro ? 18 : 2, width:16, height:16, borderRadius:"50%", background:"#fff", boxShadow:"0 1px 3px rgba(0,0,0,0.4)", transition:"left 0.2s ease-in-out" }} />
                    </button>
                  </div>
                )}
                {alleleMode && hasCi && (
                  <div className="flex items-center justify-between p-2.5 rounded-lg" style={{ background: bgItem }}>
                    <span className="text-xs" style={{ color: textMd }}>CI ribbons</span>
                    <button onClick={() => setShowCi(!showCi)} style={{ position:"relative", width:36, height:20, borderRadius:10, border:"none", cursor:"pointer", background: showCi ? "#b5860d" : "#4b5563", transition:"background 0.2s" }}>
                      <span style={{ position:"absolute", top:2, left: showCi ? 18 : 2, width:16, height:16, borderRadius:"50%", background:"#fff", boxShadow:"0 1px 3px rgba(0,0,0,0.4)", transition:"left 0.2s ease-in-out" }} />
                    </button>
                  </div>
                )}
              </div>
            )}

            {/* Min Segment Size */}
            <div>
              <div className="text-xs font-semibold uppercase tracking-wider mb-2" style={{ color: textSm }}>Heatmap Display</div>
              <div className="mb-2">
                <div className="flex justify-between items-center text-xs mb-1">
                  <span style={{ color: textSm }}>Min Segment Size</span>
                  <div className="flex items-center gap-1 rounded-md px-2 py-0.5"
                    style={{ background: bgItem, border: `1px solid ${border}` }}>
                    <input
                      type="text" inputMode="numeric" pattern="[0-9]*"
                      value={minSegmentMb === 0 ? "" : String(minSegmentMb)}
                      placeholder="0"
                      onChange={e => {
                        const raw = e.target.value.replace(/[^0-9]/g, "");
                        const v = raw === "" ? 0 : Math.min(Math.round(nBins * binSizeMb / 2), parseInt(raw, 10));
                        setMinSegmentMb(v);
                      }}
                      className="bg-transparent border-none outline-none text-right font-mono"
                      style={{ width: 36, color: textMd, fontSize: 11 }} />
                    <span className="font-mono" style={{ color: textXs, fontSize: 10 }}>Mb</span>
                  </div>
                </div>
                <input type="range" min="0" max={Math.round(nBins * binSizeMb / 2)} step={Math.ceil(binSizeMb)}
                  value={minSegmentMb}
                  onChange={e => setMinSegmentMb(+e.target.value)}
                  className="w-full h-1 rounded-full appearance-none cursor-pointer"
                  style={{ accentColor: "#f59e0b", background: border }} />
                <div className="text-xs mt-1" style={{ color: textXs }}>
                  Segments shorter than this are shown as CN=2 (Mb)
                </div>
              </div>
            </div>

            {/* Quality Filters */}
            <div>
              <div className="text-xs font-semibold uppercase tracking-wider mb-2" style={{ color: textSm }}>Quality Filters</div>
              {[["Max Residual", "residual", 3], ["Max MAPD", "mapd", 4]].map(([label, key, max]) => (
                <div key={key} className="mb-2">
                  <div className="flex justify-between text-xs mb-1">
                    <span style={{ color: textSm }}>{label}</span>
                    <span className="font-mono" style={{ color: textMd }}>{thresholds[key].toFixed(2)}</span>
                  </div>
                  <input type="range" min="0" max={max} step="0.05" value={thresholds[key]}
                    onChange={e => setThresholds(t => ({ ...t, [key]: +e.target.value }))}
                    className="w-full h-1 rounded-full appearance-none cursor-pointer"
                    style={{ accentColor: "#3b82f6", background: border }} />
                </div>
              ))}
              {hasCoverage && coverageStats && (
                <div className="mb-2">
                  <div className="flex justify-between text-xs mb-1">
                    <span style={{ color: textSm }}>Min Coverage</span>
                    <span className="font-mono" style={{ color: textMd }}>
                      {thresholds.coverage != null ? `${thresholds.coverage}×` : "off"}
                    </span>
                  </div>
                  <input type="range" min={coverageStats.min} max={coverageStats.max} step="1"
                    value={thresholds.coverage ?? coverageStats.min}
                    onChange={e => setThresholds(t => ({ ...t, coverage: +e.target.value }))}
                    className="w-full h-1 rounded-full appearance-none cursor-pointer"
                    style={{ accentColor: "#10b981", background: border }} />
                  <button className="text-xs mt-1" style={{ color: textXs }}
                    onClick={() => setThresholds(t => ({ ...t, coverage: thresholds.coverage != null ? null : coverageStats.min }))}>
                    {thresholds.coverage != null ? "disable" : "enable"}
                  </button>
                </div>
              )}
              {hasBinsWithCna && nBins > 0 && (
                <div className="mb-2">
                  <div className="flex justify-between text-xs mb-1">
                    <span style={{ color: textSm }}>Max Bins w/ CNA</span>
                    <span className="font-mono" style={{ color: textMd }}>
                      {thresholds.bins_with_cna != null ? `${thresholds.bins_with_cna} / ${nBins}` : "off"}
                    </span>
                  </div>
                  <input type="range" min="0" max={nBins} step="1"
                    value={thresholds.bins_with_cna ?? nBins}
                    onChange={e => setThresholds(t => ({ ...t, bins_with_cna: +e.target.value }))}
                    className="w-full h-1 rounded-full appearance-none cursor-pointer"
                    style={{ accentColor: "#8b5cf6", background: border }} />
                  <button className="text-xs mt-1" style={{ color: textXs }}
                    onClick={() => setThresholds(t => ({ ...t, bins_with_cna: thresholds.bins_with_cna != null ? null : nBins }))}>
                    {thresholds.bins_with_cna != null ? "disable" : "enable"}
                  </button>
                </div>
              )}
              <div className="text-xs" style={{ color: textXs }}>{filteredCells.length} / {totalCells} cells ({passRate}%)</div>
            </div>

            {/* Cell Type Filter */}
            {hasCellTypes && (
              <div>
                <div className="text-xs font-semibold uppercase tracking-wider mb-2" style={{ color: textSm }}>Cell Type</div>
                <div className="flex flex-wrap gap-1">
                  {["All", ...cellTypeOptions].map(ct => (
                    <button key={ct} onClick={() => setCellTypeFilter(ct)}
                      className="text-xs px-2 py-1 rounded-md transition-colors"
                      style={{
                        background: cellTypeFilter === ct ? "#3b82f6" : bgItem,
                        color: cellTypeFilter === ct ? "#ffffff" : textSm,
                      }}>
                      {ct}
                    </button>
                  ))}
                </div>
              </div>
            )}

            {/* Search */}
            <div className="relative">
              <svg className="absolute left-2.5 top-2.5" style={{ color: textXs }} width="14" height="14" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2"><circle cx="11" cy="11" r="8" /><path d="M21 21l-4.35-4.35" /></svg>
              <input type="text" placeholder="Search barcode..." value={search} onChange={e => setSearch(e.target.value)}
                className="w-full pl-8 pr-3 py-2 rounded-lg text-xs font-mono border-0 outline-none focus:ring-1 focus:ring-blue-500"
                style={{ background: bgItem, color: text }} />
            </div>

            {/* Cell List */}
            <div>
              <div className="text-xs font-semibold uppercase tracking-wider mb-1" style={{ color: textSm }}>Cell List</div>
              <div className="rounded-lg overflow-hidden border" style={{ borderColor: border, maxHeight: 200, overflowY: "auto", background: lightMode ? "#f8fafc" : "#222222" }}>
                {filteredCells.length === 0
                  ? <div className="p-3 text-xs text-center" style={{ color: textXs }}>No cells match</div>
                  : filteredCells.map(c => (
                    <button key={c} onClick={() => setSelectedCell(c)}
                      className="w-full text-left px-3 py-1.5 text-xs font-mono truncate transition-colors"
                      style={{
                        background: c === selectedCell ? "rgba(59,130,246,0.15)" : "transparent",
                        color: c === selectedCell ? "#60a5fa" : textSm,
                      }}>{c}
                      {hasCellTypes && data.cell_types?.[c] && (
                        <span className="ml-1 opacity-50 text-xs">[{data.cell_types[c]}]</span>
                      )}
                    </button>
                  ))}
              </div>
            </div>

            {/* Nav */}
            <div className="flex gap-2">
              <button onClick={() => navigateCell(-1)} className="flex-1 py-1.5 text-xs rounded-lg transition-colors"
                style={{ background: lightMode ? "#e2e8f0" : "#b5860d", color: lightMode ? "#374151" : "#1a1a1a", fontWeight: 600 }}>&uarr; Prev</button>
              <button onClick={() => { if (selectedCell) setTab("profile"); }}
                className="flex-1 py-1.5 text-xs rounded-lg bg-blue-600 text-white hover:bg-blue-500">View</button>
              <button onClick={() => navigateCell(1)} className="flex-1 py-1.5 text-xs rounded-lg transition-colors"
                style={{ background: lightMode ? "#e2e8f0" : "#b5860d", color: lightMode ? "#374151" : "#1a1a1a", fontWeight: 600 }}>Next &darr;</button>
            </div>

            {/* Selected metrics */}
            {selectedCell && tab !== "profile" && (
              <div>
                <div className="text-xs font-semibold uppercase tracking-wider mb-1" style={{ color: textSm }}>Selected</div>
                <div className="text-xs font-mono text-blue-400 mb-1 truncate">{selectedCell}</div>
                {hasCellTypes && data.cell_types?.[selectedCell] && (
                  <div className="text-xs mb-2" style={{ color: textXs }}>Type: {data.cell_types[selectedCell]}</div>
                )}
                <CellMetrics data={data} cellName={selectedCell} thresholds={thresholds} lightMode={lightMode} effectiveQuality={effectiveQuality} />
              </div>
            )}

            {/* Hints */}
            <div className="text-xs space-y-0.5 pt-2 border-t" style={{ borderColor: border, color: textXs }}>
              <div><kbd className="px-1 py-0.5 rounded" style={{ background: bgItem, color: textXs }}>&uarr;&darr;</kbd> Navigate</div>
              <div><kbd className="px-1 py-0.5 rounded" style={{ background: bgItem, color: textXs }}>drag</kbd> Brush-zoom</div>
              <div><kbd className="px-1 py-0.5 rounded" style={{ background: bgItem, color: textXs }}>dblclick</kbd> / <kbd className="px-1 py-0.5 rounded" style={{ background: bgItem, color: textXs }}>Esc</kbd> Reset</div>
            </div>
          </div>
        </aside>}

        {/* Main */}
        <main className="flex-1 overflow-y-auto p-5" style={{ background: bg }}>
          {tab === "heatmap" && (
            <div className="space-y-4">
              <div className="flex items-center justify-between flex-wrap gap-2">
                <div className="flex items-center gap-3">
                  <h2 className="text-sm font-semibold" style={{ color: textMd }}>{alleleMode ? "Allele-Specific Heatmap" : "Copy Number Heatmap"}</h2>
                  {isZoomed && (
                    <button onClick={() => setZoom(null)} className="text-xs px-2.5 py-1 rounded-md bg-blue-600/20 text-blue-300 hover:bg-blue-600/30 flex items-center gap-1">
                      <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2"><path d="M15 3h6v6M9 21H3v-6M21 3l-7 7M3 21l7-7" /></svg>Reset Zoom
                    </button>
                  )}
                </div>
                <div className="flex items-center gap-3">
                  <span className="text-xs" style={{ color: textXs }}>{filteredCells.length} cells</span>
                  <input type="range" min="250" max="900" step="50" value={heatmapH} onChange={e => setHeatmapH(+e.target.value)}
                    className="w-20 h-1 rounded-full appearance-none cursor-pointer" style={{ accentColor: "#3b82f6", background: border }} />
                  <DownloadMenu onDownload={fmt => heatmapPanelRef.current?.download("heatmap", fmt)} style={{ background: bgItem, color: textSm }} />
                </div>
              </div>
              <div className="rounded-xl overflow-hidden border" style={{ borderColor: border, background: bg2 }}>
                <HeatmapPanel ref={heatmapPanelRef} data={heatmapData} cellOrder={filteredCells} chrInfo={data.chr_info}
                  selectedCell={selectedCell} onCellClick={setSelectedCell} height={heatmapH}
                  alleleMode={alleleMode} zoom={zoom} onZoomChange={setZoom} showDendro={showDendro} lightMode={lightMode} />
              </div>
              {isZoomed && (
                <div className="text-xs flex items-center gap-2" style={{ color: textXs }}>
                  <svg width="12" height="12" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2"><circle cx="11" cy="11" r="8" /><path d="M21 21l-4.35-4.35M11 8v6M8 11h6" /></svg>
                  Drag to zoom further &middot; Double-click or Esc to reset
                </div>
              )}
              {selectedCell && (
                <div className="rounded-xl border p-4" style={{ borderColor: border, background: bgCard }}>
                  <div className="flex items-center gap-2 mb-2">
                    <span className="text-xs font-semibold uppercase tracking-wider" style={{ color: textSm }}>Profile</span>
                    <span className="text-xs font-mono text-blue-400">{selectedCell}</span>
                    {alleleMode && <span className="text-xs" style={{ color: "#b5860d" }}>(allele-specific)</span>}
                  </div>
                  <ProfilePlot data={data} cellName={selectedCell} height={220} alleleMode={alleleMode} lightMode={lightMode} showCi={showCi} />
                </div>
              )}
            </div>
          )}

          {tab === "profile" && (
            <div className="space-y-4">
              <div className="flex items-center gap-3">
                <h2 className="text-sm font-semibold" style={{ color: textMd }}>Cell Profile</h2>
                {selectedCell && <span className="text-xs font-mono text-blue-400 px-2 py-0.5 rounded-md" style={{ background: bgItem }}>{selectedCell}</span>}
                {alleleMode && <span className="text-xs" style={{ color: "#b5860d" }}>(allele-specific)</span>}
                <DownloadMenu onDownload={fmt => profilePlotRef.current?.download(selectedCell || "cell_profile", fmt)} style={{ background: bgItem, color: textSm }} />
              </div>
              <div className="rounded-xl border p-4" style={{ borderColor: border, background: bgCard }}>
                <ProfilePlot ref={profilePlotRef} data={data} cellName={selectedCell} height={360} alleleMode={alleleMode} lightMode={lightMode} showCi={showCi} />
              </div>
              {selectedCell && (
                <div className="rounded-xl border p-4" style={{ borderColor: border, background: lightMode ? bgCard : "#252525" }}>
                  <div className="text-xs font-semibold uppercase tracking-wider mb-3" style={{ color: textSm }}>Quality Metrics</div>
                  <CellMetrics data={data} cellName={selectedCell} thresholds={thresholds} lightMode={lightMode} effectiveQuality={effectiveQuality} />
                </div>
              )}
            </div>
          )}

          {tab === "quality" && (
            <div className="space-y-4">
              <h2 className="text-sm font-semibold" style={{ color: textMd }}>Quality Overview</h2>
              <div className="rounded-xl border p-4" style={{ borderColor: border, background: bgCard }}>
                <div className="text-xs mb-2" style={{ color: textSm }}>MAPD vs Median Residual</div>
                <QualityScatter data={data} thresholds={thresholds} selectedCell={selectedCell} onCellClick={setSelectedCell} height={400} lightMode={lightMode} />
              </div>
              <div className="grid grid-cols-3 gap-3">
                {[["Total cells", totalCells, "text-blue-400"], ["Pass filter", filteredCells.length, "text-emerald-400"], ["Pass rate", passRate + "%", "text-blue-400"]].map(([l, v, c]) => (
                  <div key={l} className="rounded-xl border p-4" style={{ borderColor: border, background: bgCard }}>
                    <div className="text-xs mb-1" style={{ color: textXs }}>{l}</div>
                    <div className={`text-2xl font-bold font-mono ${c}`}>{v}</div>
                  </div>
                ))}
              </div>
            </div>
          )}
        </main>
      </div>
    </div>
  );
}
