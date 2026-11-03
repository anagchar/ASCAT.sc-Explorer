import { useEffect, useMemo, useRef, useState } from "react";
import { axes, axisTicks, css, emptyMsg, esc, fit, FONT, hexRGB, legendRow, mixRGB, niceTicks, roundTop, tip } from "./chart";
import { fmtMetric, spearman } from "./qc";

export function useWidth(ref, fallback = 600) {
  const [w, setW] = useState(fallback);
  useEffect(() => {
    if (!ref.current) return;
    const ro = new ResizeObserver((e) => { const v = e[0]?.contentRect.width; if (v > 0) setW(Math.round(v)); });
    ro.observe(ref.current);
    return () => ro.disconnect();
  }, [ref]);
  return w;
}

function statusHTML(pass, reason) {
  return pass ? `<span class="status" style="color:var(--good)">● passes QC</span>`
    : `<span class="status" style="color:var(--serious)">▲ flagged${reason ? ": " + esc(reason) : ""}</span>`;
}

/* ---------- Explore: any metric against any other ---------- */
function Scatter({ cells, mx, my, vals, passing, shown, reasons, sel, onSelect, thresholds, theme }) {
  const wrap = useRef(), cv = useRef(), pts = useRef([]), hover = useRef(null);
  const W = useWidth(wrap);
  useEffect(() => {
    const H = 340, g = fit(cv.current, W, H); g.clearRect(0, 0, W, H);
    if (!mx || !my) return emptyMsg(g, W, H, "No metric varies between the cells.");
    const xs = vals[mx.key], ys = vals[my.key];
    const okx = xs.filter((v) => !isNaN(v) && (!mx.log || v > 0)), oky = ys.filter((v) => !isNaN(v) && (!my.log || v > 0));
    if (!okx.length || !oky.length) return emptyMsg(g, W, H, "No cell has both values.");
    const range = (v, log) => {
      let a = Math.min(...v), b = Math.max(...v); if (a === b) { a -= log ? a * 0.5 : 0.5; b += log ? b * 0.5 : 0.5; }
      if (log) return [a / 1.15, b * 1.15]; const pad = (b - a) * 0.06; return [a - pad, b + pad];
    };
    const [x0, x1] = range(okx, mx.log), [y0, y1] = range(oky, my.log);
    const box = { l: 52, t: 30, w: W - 66, h: H - 74 };
    const tk = (a, b, log) => (log ? axisTicks(a, b) : niceTicks(a, b, 5));
    const f = (m) => (v) => fmtMetric(m.key, v).replace(/\.0+(?=[kM]?$)|(\.\d*?)0+(?=[kM]?$)/, "$1");
    const { tx, ty } = axes(g, box, tk(x0, x1, mx.log), tk(y0, y1, my.log), { x0, x1, y0, y1, xlog: mx.log, ylog: my.log, xlab: mx.label, ylab: my.label, xfmt: f(mx), yfmt: f(my) });
    // median ± 3 MAD of each axis, on the bad side
    g.setLineDash([4, 4]); g.strokeStyle = css("--ink3"); g.lineWidth = 1;
    const line = (v, isX) => {
      const lim = isX ? [x0, x1] : [y0, y1]; if (!(isFinite(v) && v > lim[0] && v < lim[1])) return;
      g.beginPath(); if (isX) { g.moveTo(tx(v), box.t); g.lineTo(tx(v), box.t + box.h); } else { g.moveTo(box.l, ty(v)); g.lineTo(box.l + box.w, ty(v)); } g.stroke();
    };
    for (const [m, isX] of [[mx, true], [my, false]]) for (const [v, side] of [[m.lo, "lower"], [m.hi, "upper"]]) if (m.dir === side || m.dir === "both") line(v, isX);
    // the sidebar's QC thresholds
    g.setLineDash([]); g.strokeStyle = css("--serious"); g.globalAlpha = 0.8;
    for (const [m, isX] of [[mx, true], [my, false]]) if (thresholds[m.key] != null) line(thresholds[m.key], isX);
    g.globalAlpha = 1;
    pts.current = [];
    for (const pass of [false, true]) for (let i = 0; i < cells.length; i++) {
      if (passing[i] !== pass || isNaN(xs[i]) || isNaN(ys[i]) || (mx.log && !(xs[i] > 0)) || (my.log && !(ys[i] > 0))) continue;
      const x = tx(xs[i]), y = ty(ys[i]); pts.current.push([x, y, i]);
      g.globalAlpha = shown.has(cells[i]) ? 0.85 : 0.18; g.fillStyle = css(pass ? "--series1" : "--series2");
      g.beginPath(); if (pass) g.arc(x, y, 3.2, 0, 7); else { g.moveTo(x, y - 4); g.lineTo(x + 3.6, y + 2.6); g.lineTo(x - 3.6, y + 2.6); } g.fill();
    }
    g.globalAlpha = 1;
    const si = sel != null ? cells.indexOf(sel) : -1;
    if (si >= 0 && !isNaN(xs[si]) && !isNaN(ys[si])) {
      const x = tx(xs[si]), y = ty(ys[si]);
      g.strokeStyle = css("--panel"); g.lineWidth = 2; g.beginPath(); g.arc(x, y, 7, 0, 7); g.stroke();
      g.strokeStyle = css("--accent"); g.beginPath(); g.arc(x, y, 5.5, 0, 7); g.stroke();
    }
    const items = [["--series1", "passes QC", "dot"], ["--series2", "flagged", "tri"], ["--ink3", "median ± 3 MAD", "line"]];
    if (thresholds[mx.key] != null || thresholds[my.key] != null) items.push(["--serious", "threshold", "line"]);
    legendRow(g, box.l, 12, items);
  }, [W, cells, mx, my, vals, passing, shown, sel, thresholds, theme]);
  const move = (e) => {
    const r = cv.current.getBoundingClientRect(), x = e.clientX - r.left, y = e.clientY - r.top;
    let best = null, bd = 64; for (const [px, py, i] of pts.current) { const d = (px - x) ** 2 + (py - y) ** 2; if (d < bd) { bd = d; best = i; } }
    cv.current.style.cursor = best == null ? "default" : "pointer"; hover.current = best;
    tip(e, best == null ? null : `<b class="mono">${esc(cells[best])}</b><br>${statusHTML(passing[best], reasons[best])}<br><span class="k">${esc(mx.label)}</span> ${fmtMetric(mx.key, vals[mx.key][best])}<br><span class="k">${esc(my.label)}</span> ${fmtMetric(my.key, vals[my.key][best])}`);
  };
  return <div ref={wrap}><canvas ref={cv} onMouseMove={move} onMouseLeave={(e) => tip(e, null)} onClick={() => hover.current != null && onSelect(cells[hover.current])} /></div>;
}

/* ---------- Spearman correlation matrix ---------- */
function Correlations({ met, vals, onPick, theme }) {
  const wrap = useRef(), cv = useRef(), geo = useRef(null);
  const W = useWidth(wrap);
  const rho = useMemo(() => met.map((a) => met.map((b) => (a === b ? 1 : spearman(vals[a.key], vals[b.key])))), [met, vals]);
  useEffect(() => {
    const m = met.length, left = 136, top = 8;
    const cell = Math.max(14, Math.min(34, Math.floor((W - left - 12) / Math.max(1, m)))), H = top + cell * m + 132;
    const g = fit(cv.current, W, H); g.clearRect(0, 0, W, H);
    if (!m) return;
    const neg = hexRGB(css("--cn0")), pos = hexRGB(css("--cn4")), mid = hexRGB(css("--panel2"));
    geo.current = { left, top, cell, m };
    g.font = `11px ${FONT}`;
    for (let r = 0; r < m; r++) for (let c = 0; c < m; c++) {
      const v = rho[r][c], x = left + c * cell, y = top + r * cell;
      g.fillStyle = r === c ? css("--field") : isNaN(v) ? css("--cnna") : v >= 0 ? mixRGB(mid, pos, v) : mixRGB(mid, neg, -v);
      g.beginPath(); if (g.roundRect) g.roundRect(x + 1, y + 1, cell - 2, cell - 2, 3); else g.rect(x + 1, y + 1, cell - 2, cell - 2); g.fill();
      if (cell >= 28 && !isNaN(v) && r !== c) {
        g.fillStyle = Math.abs(v) > 0.6 ? "#fff" : css("--ink2"); g.textAlign = "center"; g.textBaseline = "middle"; g.font = `10px ${FONT}`;
        g.fillText(v.toFixed(1).replace("-0.0", "0.0"), x + cell / 2, y + cell / 2); g.font = `11px ${FONT}`;
      }
    }
    g.fillStyle = css("--ink2"); g.textAlign = "right"; g.textBaseline = "middle";
    met.forEach((mt, r) => g.fillText(mt.label, left - 8, top + r * cell + cell / 2));
    met.forEach((mt, c) => { g.save(); g.translate(left + c * cell + cell / 2, top + m * cell + 6); g.rotate(-Math.PI / 4); g.fillText(mt.label, 0, 0); g.restore(); });
    const kx = left, ky = H - 12, kw = Math.min(220, cell * m);
    for (let k = 0; k < kw; k++) { const t = (k / (kw - 1)) * 2 - 1; g.fillStyle = t >= 0 ? mixRGB(mid, pos, t) : mixRGB(mid, neg, -t); g.fillRect(kx + k, ky, 1, 8); }
    g.fillStyle = css("--ink3"); g.textBaseline = "bottom";
    g.textAlign = "left"; g.fillText("−1", kx, ky - 2); g.textAlign = "center"; g.fillText("0", kx + kw / 2, ky - 2); g.textAlign = "right"; g.fillText("+1", kx + kw, ky - 2);
  }, [W, met, rho, theme]);
  const hit = (e) => {
    const c = geo.current; if (!c) return null;
    const r = cv.current.getBoundingClientRect(), x = e.clientX - r.left - c.left, y = e.clientY - r.top - c.top;
    const i = Math.floor(y / c.cell), j = Math.floor(x / c.cell);
    return i >= 0 && j >= 0 && i < c.m && j < c.m && i !== j ? [i, j] : null;
  };
  const move = (e) => {
    const h = hit(e); cv.current.style.cursor = h ? "pointer" : "default";
    tip(e, h ? `<b>${esc(met[h[0]].label)}</b> × <b>${esc(met[h[1]].label)}</b><br><span class="k">Spearman ρ</span> ${rho[h[0]][h[1]].toFixed(2)}` : null);
  };
  return <div ref={wrap}><canvas ref={cv} onMouseMove={move} onMouseLeave={(e) => tip(e, null)} onClick={(e) => { const h = hit(e); if (h) onPick(met[h[1]].key, met[h[0]].key); }} /></div>;
}

/* ---------- one metric's histogram, passing and flagged stacked ---------- */
function Dist({ m, values, passing, on, onClick, theme }) {
  const cv = useRef(), wrap = useRef();
  const W = useWidth(wrap, 200);
  useEffect(() => {
    const H = 92, g = fit(cv.current, W, H); g.clearRect(0, 0, W, H);
    const sc = values.map((v) => (m.log ? (v > 0 ? Math.log10(v) : NaN) : v)), ok = sc.filter((v) => !isNaN(v));
    if (!ok.length) return;
    let lo = Math.min(...ok), hi = Math.max(...ok); if (lo === hi) { lo -= 0.5; hi += 0.5; }
    const nb = 24, P = new Array(nb).fill(0), F = new Array(nb).fill(0);
    sc.forEach((v, i) => { if (isNaN(v)) return; const k = Math.min(nb - 1, Math.floor(((v - lo) / (hi - lo)) * nb)); (passing[i] ? P : F)[k]++; });
    const ymax = Math.max(...P.map((p, k) => p + F[k])) || 1, box = { l: 2, t: 4, w: W - 4, h: H - 22 }, bw = box.w / nb;
    for (let k = 0; k < nb; k++) {
      let y = box.t + box.h; const x = box.l + k * bw + 1, w = Math.max(1, bw - 2);
      for (const [n, c] of [[P[k], "--series1"], [F[k], "--series2"]]) {
        if (!n) continue; const h = box.h * (n / ymax); g.fillStyle = css(c);
        roundTop(g, x, y - h, w, h - (y < box.t + box.h ? 2 : 0), 2); y -= h;
      }
    }
    const X = (v) => box.l + ((v - lo) / (hi - lo)) * box.w, tr = (v) => (m.log ? Math.log10(v) : v), back = (v) => (m.log ? Math.pow(10, v) : v);
    g.setLineDash([3, 3]); g.strokeStyle = css("--ink3"); g.lineWidth = 1;
    for (const [v, side] of [[m.lo, "lower"], [m.hi, "upper"]]) {
      if (!isFinite(v) || !(m.dir === side || m.dir === "both")) continue;
      const x = X(tr(v)); if (x < box.l || x > box.l + box.w) continue; g.beginPath(); g.moveTo(x, box.t); g.lineTo(x, box.t + box.h); g.stroke();
    }
    g.setLineDash([]); g.strokeStyle = css("--ink"); g.lineWidth = 1.5;
    const xm = X(tr(m.median)); g.beginPath(); g.moveTo(xm, box.t); g.lineTo(xm, box.t + box.h); g.stroke();
    g.fillStyle = css("--ink3"); g.font = `10.5px ${FONT}`; g.textBaseline = "top";
    // the median's label, then the range's ends where they leave it room
    const ml = "median " + fmtMetric(m.key, m.median), mw = g.measureText(ml).width;
    const mc = Math.min(Math.max(xm, box.l + mw / 2), box.l + box.w - mw / 2), ly = box.t + box.h + 4;
    const lol = fmtMetric(m.key, back(lo)), hil = fmtMetric(m.key, back(hi));
    if (box.l + g.measureText(lol).width + 6 < mc - mw / 2) { g.textAlign = "left"; g.fillText(lol, box.l, ly); }
    if (box.l + box.w - g.measureText(hil).width - 6 > mc + mw / 2) { g.textAlign = "right"; g.fillText(hil, box.l + box.w, ly); }
    g.textAlign = "center"; g.fillStyle = css("--ink2"); g.fillText(ml, mc, ly);
  }, [W, m, values, passing, theme]);
  return (
    <button className={"dist" + (on ? " on" : "")} title={m.help + " Click to plot it on the x axis."} onClick={onClick}>
      <span className="dn">{m.label}</span><span className="dg">{m.group}</span>
      <span ref={wrap} style={{ display: "block" }}><canvas ref={cv} /></span>
    </button>
  );
}

/**
 * The Cell QC tab: flags, an explorer of any two metrics, their correlations,
 * the distribution of each and a table of robust statistics.
 *  - cells: cell names; vals: { metric: number[] } aligned with cells;
 *    metAll: metricSummaries(vals)
 *  - passing / reasons: per cell, from the sidebar's QC thresholds
 *  - shown: names of the cells the filters keep; thresholds: { metric: value }
 */
export default function QCTab({ cells, vals, metAll, passing, reasons, shown, tiles, thresholds, sel, onSelect, theme, children }) {
  const met = useMemo(() => metAll.filter((m) => new Set(vals[m.key].filter((v) => !isNaN(v))).size > 1), [metAll, vals]);
  const [ax, setAx] = useState({ x: null, y: null });
  const has = (k) => met.some((m) => m.key === k);
  const xKey = ax.x && has(ax.x) ? ax.x : has("reads_per_bin") ? "reads_per_bin" : has("mapd") ? "mapd" : met[0]?.key;
  const yKey = ax.y && has(ax.y) ? ax.y : has("mapd") && xKey !== "mapd" ? "mapd" : has("median_residual") ? "median_residual" : met[1]?.key ?? xKey;
  const mx = met.find((m) => m.key === xKey), my = met.find((m) => m.key === yKey);
  const groups = [...new Set(metAll.map((m) => m.group))];
  const dirTxt = { upper: "high", lower: "low", both: "either side" };

  if (!metAll.length) return <><div className="card"><div className="empty">This dataset has no per-cell QC metrics.</div></div>{children}</>;
  return (
    <>
      <div className="tiles">
        {tiles.map(([k, v, s]) => <div className="tile" key={k}><div className="k">{k}</div><div className="v">{v}</div><div className="s">{s}</div></div>)}
      </div>
      <div className="grid2">
        <div className="card">
          <div className="hd"><h2>Explore</h2><span className="sub">each dot is a cell · click to open it</span>
            <div className="axsel">
              <label>X <select aria-label="x axis metric" value={xKey || ""} onChange={(e) => setAx((a) => ({ ...a, x: e.target.value }))}>{met.map((m) => <option key={m.key} value={m.key}>{m.label}</option>)}</select></label>
              <label>Y <select aria-label="y axis metric" value={yKey || ""} onChange={(e) => setAx((a) => ({ ...a, y: e.target.value }))}>{met.map((m) => <option key={m.key} value={m.key}>{m.label}</option>)}</select></label>
            </div>
          </div>
          <div className="bd"><Scatter cells={cells} mx={mx} my={my} vals={vals} passing={passing} shown={shown} reasons={reasons} sel={sel} onSelect={onSelect} thresholds={thresholds} theme={theme} /></div>
        </div>
        <div className="card">
          <div className="hd"><h2>Metric correlations</h2><span className="sub">Spearman ρ · click a square to plot the pair</span></div>
          <div className="bd"><Correlations met={met} vals={vals} theme={theme} onPick={(x, y) => setAx({ x, y })} /></div>
        </div>
      </div>
      {children}
      <div className="card">
        <div className="hd"><h2>Distributions</h2><span className="sub">all cells, passing and flagged · solid: median · dashed: median ± 3 MAD on the bad side · click one to plot it</span></div>
        <div className="bd"><div className="dists">
          {met.map((m) => <Dist key={m.key} m={m} values={vals[m.key]} passing={passing} on={m.key === xKey} theme={theme} onClick={() => setAx((a) => ({ ...a, x: m.key }))} />)}
        </div></div>
      </div>
      <div className="card">
        <div className="hd"><h2>Statistics</h2><span className="sub">per metric, all cells · metrics of ASCAT.sc's getQCs()</span></div>
        <div className="bd tablewrap"><table>
          <thead><tr><th>metric</th><th>cells</th><th>median</th><th>MAD</th><th>min</th><th>max</th><th title="cells beyond median ± 3 MAD on the bad side (log scale for reads)">outliers</th></tr></thead>
          <tbody>
            {groups.map((gname) => [
              <tr className="grp" key={gname}><td colSpan={7}>{gname}</td></tr>,
              ...metAll.filter((m) => m.group === gname).map((m) => (
                <tr key={m.key}>
                  <td title={m.help}>{m.label} <span className="ok">· bad when {dirTxt[m.dir]}</span></td>
                  <td>{m.n.toLocaleString()}</td><td>{fmtMetric(m.key, m.median)}</td>
                  <td>{m.log ? <span className="ok">log</span> : fmtMetric(m.key, m.mad)}</td>
                  <td>{fmtMetric(m.key, m.min)}</td><td>{fmtMetric(m.key, m.max)}</td>
                  <td>{m.out ? <><span style={{ color: "var(--serious)", fontWeight: 600 }}>{m.out}</span> <span className="ok">({((100 * m.out) / m.n).toFixed(1)}%)</span></> : <span className="ok">0</span>}</td>
                </tr>
              )),
            ])}
          </tbody>
        </table></div>
      </div>
    </>
  );
}
