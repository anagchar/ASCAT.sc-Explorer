// Per-cell QC metrics of ASCAT.sc's getQCs() (res$QC_metrics), as the
// ASCAT.sc-py viewer shows them (ascatsc/qcmetrics.py).
//
// When the JSON carries `qc_metrics` (exported from res$QC_metrics) those
// values are used as they are. Otherwise they are computed here from what the
// export has: the per-bin `raw` track, the integer `profiles` and the `reads`.
// Segments are the runs of constant copy number within a chromosome, and the
// solution is taken at purity 1 with the profile's mean copy number as ploidy
// (single cells), so these can differ slightly from getQCs() on the full fit.

import { quantile } from "./chart";

// [key, label, group, direction of a bad value, log scale, help]
export const METRIC_META = [
  ["n_reads_total", "Total reads", "Coverage", "lower", true, "Reads in the cell's bins. Few reads give noisy profiles."],
  ["reads_per_bin", "Reads per bin", "Coverage", "lower", true, "A cell's mean reads per bin."],
  ["coverage", "Coverage", "Coverage", "lower", false, "ASCAT.sc's coverage of the cell (allSolutions)."],
  ["gini", "Gini coefficient", "Coverage", "upper", false, "Unevenness of the reads across bins: 0 is uniform, towards 1 reads pile up in few bins."],
  ["mapd", "MAPD", "Noise", "upper", false, "Median absolute difference of neighbouring bins' GC-corrected logR, within chromosomes."],
  ["MAPD_gc_corrected", "MAPD (getQCs)", "Noise", "upper", false, "MAPD as getQCs() computes it."],
  ["spikiness", "Spikiness", "Noise", "upper", false, "Median second difference of logR: isolated outlier bins. Lower is smoother."],
  ["MBRSM_dispersion", "Segment dispersion", "Noise", "upper", false, "Scaled MAD of the bins' logR around their segment."],
  ["autocorr", "Autocorrelation", "Noise", "upper", false, "Lag-1 autocorrelation of segment residuals; high values suggest waves (GC, replication timing)."],
  ["median_residual", "Median residual", "Fit", "upper", false, "Median distance of the bins' copy number from the fitted integer profile."],
  ["distance_integer_logR", "LogR integer distance", "Fit", "upper", false, "Distance of segment copy numbers from the nearest integer. Lower is a better fit."],
  ["distance_integer_BAF", "BAF integer distance", "Fit", "upper", false, "Distance of segment BAFs from the nearest integer allelic state. Lower is cleaner."],
  ["MSRSI_non_integerness", "Non-integerness", "Fit", "upper", false, "Scaled MAD of segment logR from the logR its integer state predicts."],
  ["breakpoints", "Breakpoints", "Fit", "upper", false, "Copy-number state changes; very many usually means noise was segmented."],
  ["bins_with_cna", "Bins with CNA", "Solution", "upper", false, "Bins whose copy number is not 2 (after the minimum segment size)."],
  ["state_mode", "Modal CN state", "Solution", "both", false, "Most frequent copy-number state; usually 2."],
  ["ploidy", "Ploidy", "Solution", "both", false, "Mean copy number over the bins."],
  ["purity", "Purity", "Solution", "both", false, "Fitted purity; single cells are expected near 1."],
];
export const INT_METRICS = new Set(["breakpoints", "state_mode", "bins_with_cna", "n_reads_total"]);

const madScaled = (x) => {
  const v = x.filter((a) => !isNaN(a)).sort((a, b) => a - b);
  if (!v.length) return NaN;
  const m = quantile(v, 0.5);
  return 1.4826 * quantile(v.map((a) => Math.abs(a - m)).sort((a, b) => a - b), 0.5);
};
const medianOf = (v) => { const s = v.filter((a) => !isNaN(a)).sort((a, b) => a - b); return s.length ? quantile(s, 0.5) : NaN; };
const num = (v) => (v == null || !isFinite(v) ? NaN : +v);

function corr(a, b) {
  const n = a.length; let ma = 0, mb = 0;
  for (let i = 0; i < n; i++) { ma += a[i]; mb += b[i]; }
  ma /= n; mb /= n;
  let sab = 0, saa = 0, sbb = 0;
  for (let i = 0; i < n; i++) { const da = a[i] - ma, db = b[i] - mb; sab += da * db; saa += da * da; sbb += db * db; }
  return saa > 0 && sbb > 0 ? sab / Math.sqrt(saa * sbb) : NaN;
}

/** getQCs()-style metrics for every cell: { [cell]: { metric: value } }. */
export function computeQCMetrics(data) {
  const cells = Object.keys(data.profiles || {});
  const chr = data.bins.chr, nb = chr.length;
  const slices = [];
  for (let a = 0, i = 1; i <= nb; i++) if (i === nb || chr[i] !== chr[i - 1]) { slices.push([a, i]); a = i; }
  const widthMb = data.bins.start.map((s, i) => (data.bins.end[i] - s) / 1e6);
  // The R export stores raw as ploidy * 2^logR (copy-number scale, always > 0);
  // the Python one stores logR itself, which has negative values.
  const logStyle = !!data.raw && cells.slice(0, 20).some((c) => (data.raw[c] || []).some((v) => v != null && v < 0));

  const out = {};
  for (const cell of cells) {
    const r = {};
    const reads = data.reads?.[cell];
    if (reads) {
      let s = 0, n = 0; for (const v of reads) if (v != null) { s += v; n++; }
      r.n_reads_total = s; r.reads_per_bin = n ? s / n : NaN;
    }
    const prof = data.profiles[cell];
    let ps = 0, pn = 0; for (const v of prof) if (v != null && isFinite(v)) { ps += v; pn++; }
    const psi = pn ? ps / pn : NaN;
    r.ploidy = psi;
    // breakpoints and modal state, from the integer profile
    let bp = 0; const counts = new Map();
    for (const [a, b] of slices) for (let i = a; i < b; i++) {
      const v = prof[i]; if (v == null) continue;
      counts.set(v, (counts.get(v) || 0) + 1);
      if (i > a && prof[i - 1] != null && prof[i - 1] !== v) bp++;
    }
    r.breakpoints = bp;
    let mode = NaN, best = -1;
    for (const [v, k] of [...counts].sort((x, y) => x[0] - y[0])) if (k > best) { best = k; mode = v; }
    r.state_mode = mode;

    const raw = data.raw?.[cell];
    if (raw && raw.length === nb && isFinite(psi) && psi > 0) {
      const L = new Float64Array(nb);
      for (let i = 0; i < nb; i++) { const v = raw[i]; L[i] = v == null ? NaN : logStyle ? v : v > 0 ? Math.log2(v / psi) : NaN; }
      const d1 = [], d2 = [];
      for (const [a, b] of slices) {
        for (let i = a + 1; i < b; i++) { const d = Math.abs(L[i] - L[i - 1]); if (!isNaN(d)) d1.push(d); }
        for (let i = a + 2; i < b; i++) { const d = Math.abs(L[i] - 2 * L[i - 1] + L[i - 2]); if (!isNaN(d)) d2.push(d); }
      }
      r.mapd = medianOf(d1); r.spikiness = medianOf(d2);
      // segments: runs of one integer state within a chromosome
      const seg = new Float64Array(nb).fill(NaN);
      let dw = 0, dsum = 0;
      for (const [a, b] of slices) {
        for (let s = a; s < b;) {
          let e = s + 1; while (e < b && prof[e] === prof[s]) e++;
          let m = 0, k = 0, w = 0; for (let i = s; i < e; i++) { if (!isNaN(L[i])) { m += L[i]; k++; } w += widthMb[i]; }
          if (k && prof[s] != null) {
            m /= k; for (let i = s; i < e; i++) seg[i] = m;
            const cnU = psi * Math.pow(2, m); dsum += w * (cnU - Math.round(cnU)) ** 2; dw += w;
          }
          s = e;
        }
      }
      if (dw > 0) r.distance_integer_logR = dsum / dw;
      const resid = Array.from(L, (v, i) => v - seg[i]);
      r.MBRSM_dispersion = madScaled(resid);
      let acs = 0, aw = 0;
      for (const [a, b] of slices) {
        if (b - a < 3) continue;
        const x = [], y = [];
        for (let i = a; i < b - 1; i++) if (!isNaN(resid[i]) && !isNaN(resid[i + 1])) { x.push(resid[i]); y.push(resid[i + 1]); }
        if (x.length < 2) continue;
        const c = corr(x, y); if (isNaN(c)) continue;
        acs += c * (b - a); aw += b - a;
      }
      if (aw) r.autocorr = acs / aw;
      // distance of segment logR from the logR of its integer state, by (chromosome, state)
      const groups = new Map();
      for (let i = 0; i < nb; i++) {
        const s = prof[i]; if (s == null || !(s > 0) || isNaN(seg[i])) continue;
        const key = chr[i] + "|" + s; if (!groups.has(key)) groups.set(key, []);
        groups.get(key).push(seg[i] - Math.log2(s / psi));
      }
      if (groups.size) r.MSRSI_non_integerness = madScaled([...groups.values()].map(medianOf));
    }
    const pre = data.qc_metrics?.[cell];
    if (pre) for (const k in pre) if (pre[k] != null && isFinite(pre[k])) r[k] = +pre[k];
    out[cell] = r;
  }
  return out;
}

/** Robust summary of a metric over all cells, with median ± 3 MAD outliers on its bad side. */
export function summarise(meta, values) {
  const sc = values.map((v) => (meta.log ? (v > 0 ? Math.log10(v) : NaN) : v));
  const ok = sc.filter((v) => !isNaN(v)).sort((a, b) => a - b);
  const med = quantile(ok, 0.5), mad = 1.4826 * quantile(ok.map((v) => Math.abs(v - med)).sort((a, b) => a - b), 0.5);
  const lo = med - 3 * mad, hi = med + 3 * mad, back = (v) => (meta.log ? Math.pow(10, v) : v);
  let outliers = 0;
  sc.forEach((v) => { if (isNaN(v) || !(mad > 0)) return; if ((meta.dir !== "lower" && v > hi) || (meta.dir !== "upper" && v < lo)) outliers++; });
  const fin = values.filter((v) => !isNaN(v));
  let min = Infinity, max = -Infinity; for (const v of fin) { if (v < min) min = v; if (v > max) max = v; }
  return { n: fin.length, median: back(med), mad: meta.log ? null : mad, lo: back(lo), hi: back(hi), min, max, out: outliers };
}

/** The metrics present in `vals` ({ metric: number[] }), each with its summary. */
export function metricSummaries(vals) {
  return METRIC_META
    .filter(([k]) => vals[k] && vals[k].some((v) => !isNaN(v)))
    .map(([key, label, group, dir, log, help]) => {
      const m = { key, label, group, dir, log, help };
      return Object.assign(m, summarise(m, vals[key]));
    });
}

/** Whether v lies beyond median ± 3 MAD on the metric's bad side. */
export function isOutlier(m, v) {
  if (v == null || isNaN(v)) return false;
  if (m.log) { if (!(v > 0)) return m.dir !== "upper"; }
  return (m.dir !== "lower" && v > m.hi) || (m.dir !== "upper" && v < m.lo);
}

export function spearman(a, b) {
  const idx = []; for (let i = 0; i < a.length; i++) if (!isNaN(a[i]) && !isNaN(b[i])) idx.push(i);
  if (idx.length < 3) return NaN;
  const rk = (v) => {
    const o = idx.map((i, j) => [v[i], j]).sort((p, q) => p[0] - q[0]), r = new Array(o.length);
    for (let s = 0; s < o.length;) { let e = s; while (e + 1 < o.length && o[e + 1][0] === o[s][0]) e++; for (let k = s; k <= e; k++) r[o[k][1]] = (s + e) / 2; s = e + 1; }
    return r;
  };
  return corr(rk(a), rk(b));
}

export const fmtMetric = (key, v) => {
  if (v == null || !isFinite(v)) return "–";
  if (key === "n_reads_total") return v < 1e3 ? String(Math.round(v)) : v < 1e6 ? (v / 1e3).toFixed(v < 1e4 ? 2 : 1) + "k" : (v / 1e6).toFixed(v < 1e7 ? 2 : 1) + "M";
  if (INT_METRICS.has(key)) return String(Math.round(v * 10) / 10);
  const a = Math.abs(v);
  return a >= 100 ? v.toFixed(0) : a >= 10 ? v.toFixed(1) : a >= 1 ? v.toFixed(2) : v.toFixed(3);
};
export { num };
