// Canvas and formatting helpers shared by every plot — the same ones the
// ASCAT.sc-py viewer uses, so charts read the theme from the CSS tokens.

export const css = (v) => getComputedStyle(document.documentElement).getPropertyValue(v).trim();
export const hexRGB = (h) => [parseInt(h.slice(1, 3), 16), parseInt(h.slice(3, 5), 16), parseInt(h.slice(5, 7), 16)];
export const mixRGB = (a, b, t) => `rgb(${a.map((v, k) => Math.round(v + (b[k] - v) * t)).join(",")})`;
export const FONT = "-apple-system, BlinkMacSystemFont, 'SF Pro Text', sans-serif";

/** Size a canvas for the device pixel ratio; returns its 2D context in CSS pixels. */
export function fit(cv, w, h) {
  const dpr = Math.max(1, window.devicePixelRatio || 1);
  cv.width = Math.round(w * dpr); cv.height = Math.round(h * dpr);
  cv.style.width = w + "px"; cv.style.height = h + "px";
  const g = cv.getContext("2d"); g.setTransform(dpr, 0, 0, dpr, 0, 0); return g;
}

export const fmt = (v, d = 2) => (v == null || !isFinite(v) ? "–" : Number(v).toFixed(d));
export const fmtInt = (v) => (v == null || !isFinite(v) ? "–" : Math.round(v).toLocaleString());
// compact counts: 950, 12.0k, 3.62M (the exact number goes in a tooltip)
export const fmtCount = (v) => (v == null || !isFinite(v) ? "–" : v < 1e3 ? String(Math.round(v)) : v < 1e6 ? (v / 1e3).toFixed(v < 1e4 ? 2 : 1) + "k" : (v / 1e6).toFixed(v < 1e7 ? 2 : 1) + "M");
export const fmtBp = (b) => (b >= 1e6 ? (b / 1e6).toFixed(b % 1e6 ? 2 : 0) + " Mb" : (b / 1e3).toFixed(0) + " kb");
export const chrLabel = (c) => String(c).replace(/^chr/, "");
export const esc = (s) => String(s).replace(/[&<>"']/g, (c) => ({ "&": "&amp;", "<": "&lt;", ">": "&gt;", '"': "&quot;", "'": "&#39;" }[c]));

export function quantile(sorted, p) {
  if (!sorted.length) return NaN;
  const i = (sorted.length - 1) * p, a = Math.floor(i), b = Math.ceil(i);
  return sorted[a] + (sorted[b] - sorted[a]) * (i - a);
}
export function niceTicks(a, b, n = 5) {
  const s = Math.pow(10, Math.floor(Math.log10((b - a) / n))), st = [1, 2, 5, 10].map((m) => m * s).find((v) => (b - a) / v <= n) || s * 10;
  const out = []; for (let v = Math.ceil(a / st) * st; v <= b + 1e-12; v += st) out.push(+v.toFixed(10)); return out;
}
export function logTicks(a, b) {
  const out = [];
  for (let e = Math.floor(Math.log10(a)); e <= Math.ceil(Math.log10(b)); e++) for (const m of [1, 2, 5]) { const v = m * Math.pow(10, e); if (v >= a && v <= b) out.push(v); }
  return out;
}
export function axisTicks(a, b) { const t = logTicks(a, b); return t.length >= 3 ? t : niceTicks(a, b, 5).filter((v) => v > 0); }

/** Gridlines, tick labels and axis titles; returns the value → pixel maps. */
export function axes(g, box, xs, ys, opt) {
  const { x0, x1, y0, y1, xlog, ylog } = opt;
  const tx = (v) => box.l + ((xlog ? Math.log10(v) - Math.log10(x0) : v - x0) / (xlog ? Math.log10(x1) - Math.log10(x0) : x1 - x0)) * box.w;
  const ty = (v) => box.t + box.h - ((ylog ? Math.log10(v) - Math.log10(y0) : v - y0) / (ylog ? Math.log10(y1) - Math.log10(y0) : y1 - y0)) * box.h;
  g.font = `11px ${FONT}`; g.lineWidth = 1;
  g.strokeStyle = css("--grid"); g.fillStyle = css("--ink3"); g.textAlign = "right"; g.textBaseline = "middle";
  for (const v of ys) { const y = ty(v); g.beginPath(); g.moveTo(box.l, y); g.lineTo(box.l + box.w, y); g.stroke(); g.fillText(opt.yfmt ? opt.yfmt(v) : v, box.l - 6, y); }
  g.textAlign = "center"; g.textBaseline = "top";
  for (const v of xs) g.fillText(opt.xfmt ? opt.xfmt(v) : v, tx(v), box.t + box.h + 6);
  g.fillStyle = css("--ink2");
  if (opt.xlab) g.fillText(opt.xlab, box.l + box.w / 2, box.t + box.h + 22);
  if (opt.ylab) { g.save(); g.translate(12, box.t + box.h / 2); g.rotate(-Math.PI / 2); g.textBaseline = "middle"; g.fillText(opt.ylab, 0, 0); g.restore(); }
  return { tx, ty };
}

/** A row of [colour token, label, shape] legend items drawn on the canvas. */
export function legendRow(g, x, y, items) {
  g.font = `11.5px ${FONT}`; g.textBaseline = "middle"; g.textAlign = "left";
  for (const [c, t, shape] of items) {
    g.fillStyle = css(c); g.beginPath();
    if (shape === "tri") { g.moveTo(x + 4, y - 4); g.lineTo(x + 8, y + 3); g.lineTo(x, y + 3); }
    else if (shape === "line") { g.fillRect(x, y - 1.5, 14, 3); }
    else g.arc(x + 4, y, 3.5, 0, 7);
    g.fill(); g.fillStyle = css("--ink2"); g.fillText(t, x + (shape === "line" ? 20 : 12), y);
    x += g.measureText(t).width + (shape === "line" ? 34 : 26);
  }
}

export function roundTop(g, x, y, w, h, r) {
  if (h <= 0) return; r = Math.min(r, w / 2, h);
  g.beginPath(); g.moveTo(x, y + h); g.lineTo(x, y + r); g.quadraticCurveTo(x, y, x + r, y); g.lineTo(x + w - r, y);
  g.quadraticCurveTo(x + w, y, x + w, y + r); g.lineTo(x + w, y + h); g.fill();
}

/** Message centred (or at a corner) of an empty chart. */
export function emptyMsg(g, W, H, text) {
  g.fillStyle = css("--ink3"); g.font = `13px ${FONT}`; g.textAlign = "center"; g.textBaseline = "middle"; g.fillText(text, W / 2, H / 2);
}

/* ---------- tooltip: one fixed element shared by every chart ---------- */
let tipEl = null;
export function tip(e, html) {
  if (!tipEl) { tipEl = document.createElement("div"); tipEl.className = "tip"; document.body.appendChild(tipEl); }
  if (!html) { tipEl.style.display = "none"; return; }
  tipEl.innerHTML = html; tipEl.style.display = "block";
  const r = tipEl.getBoundingClientRect(); let x = e.clientX + 14, y = e.clientY + 14;
  if (x + r.width > window.innerWidth - 8) x = e.clientX - r.width - 14;
  if (y + r.height > window.innerHeight - 8) y = e.clientY - r.height - 14;
  tipEl.style.left = x + "px"; tipEl.style.top = y + "px";
}

/** Save a canvas as PNG or PDF, on the card's background colour. */
export async function saveCanvas(canvases, basename, format = "png") {
  const list = canvases.filter(Boolean);
  if (!list.length) return;
  const W = list.reduce((s, c) => s + c.width, 0), H = Math.max(...list.map((c) => c.height));
  const out = document.createElement("canvas"); out.width = W; out.height = H;
  const g = out.getContext("2d"); g.fillStyle = css("--panel"); g.fillRect(0, 0, W, H);
  let x = 0; for (const c of list) { g.drawImage(c, x, 0); x += c.width; }
  if (format === "pdf") {
    const { default: jsPDF } = await import("jspdf");
    const cssW = list.reduce((s, c) => s + c.offsetWidth, 0), cssH = Math.max(...list.map((c) => c.offsetHeight));
    const pdf = new jsPDF({ orientation: cssW > cssH ? "landscape" : "portrait", unit: "px", format: [cssW, cssH], hotfixes: ["px_scaling"] });
    pdf.addImage(out.toDataURL("image/png"), "PNG", 0, 0, cssW, cssH);
    pdf.save(basename + ".pdf");
  } else {
    const a = document.createElement("a"); a.download = basename + ".png"; a.href = out.toDataURL("image/png"); a.click();
  }
}
