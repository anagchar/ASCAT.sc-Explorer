import * as d3 from "d3";

const DEMO_CELL_TYPES = ["Epithelial", "Stromal", "Immune", "Endothelial"];

// Manhattan distance + Ward.D2 hierarchical clustering, matching R's
// hclust(dist(..., method="manhattan"), method="ward.D2").
// Squared distances are stored throughout; heights are emitted as sqrt.
function buildDemoDendrogram(cells, profiles, nMajor = null, nMinor = null) {
  const n = cells.length;
  const nBins = profiles[cells[0]].length;
  const useAS = nMajor && nMinor;
  const maxNodes = 2 * n - 1;

  // Flat distance matrix (squared Manhattan)
  const dist = new Float64Array(maxNodes * maxNodes);
  const active = new Uint8Array(maxNodes);
  const nodeSize = new Int32Array(maxNodes);
  const nodeCode = new Int32Array(maxNodes);
  const setDist = (i, j, d) => { dist[i * maxNodes + j] = d; dist[j * maxNodes + i] = d; };
  const getDist = (i, j) => dist[i * maxNodes + j];

  for (let i = 0; i < n; i++) { active[i] = 1; nodeSize[i] = 1; nodeCode[i] = -(i + 1); }

  // Build per-cell vectors matching R: rbind(maj_mat, min_mat), NA -> -1
  // so the vector is [maj0..majN, min0..minN]
  const vectors = cells.map(c => {
    if (useAS && nMajor[c] && nMinor[c]) {
      const maj = nMajor[c], min = nMinor[c];
      const v = new Array(nBins * 2);
      for (let b = 0; b < nBins; b++) v[b]         = maj[b] ?? -1;
      for (let b = 0; b < nBins; b++) v[nBins + b]  = min[b] ?? -1;
      return v;
    }
    return profiles[c];
  });

  for (let i = 0; i < n; i++) {
    for (let j = i + 1; j < n; j++) {
      const a = vectors[i], b = vectors[j];
      let m = 0;
      for (let k = 0; k < a.length; k++) m += Math.abs((a[k] ?? -1) - (b[k] ?? -1));
      setDist(i, j, m * m); // store squared distance
    }
  }

  const merge = [], height = [];
  const children = new Map();
  let nextNode = n, activeCount = n;

  while (activeCount > 1) {
    let minD = Infinity, s = -1, t = -1;
    for (let i = 0; i < nextNode; i++) {
      if (!active[i]) continue;
      for (let j = i + 1; j < nextNode; j++) {
        if (!active[j]) continue;
        const d = getDist(i, j);
        if (d < minD) { minD = d; s = i; t = j; }
      }
    }
    if (s < 0 || !Number.isFinite(minD)) break;

    const stepCode = merge.length + 1;
    merge.push([nodeCode[s], nodeCode[t]]);
    height.push(Math.sqrt(minD)); // emit sqrt to match R convention
    children.set(stepCode, [nodeCode[s], nodeCode[t]]);

    const u = nextNode++;
    nodeCode[u] = stepCode;
    nodeSize[u] = nodeSize[s] + nodeSize[t];

    // Ward.D2 Lance-Williams update (on squared distances)
    const dST = getDist(s, t);
    for (let v = 0; v < u; v++) {
      if (!active[v] || v === s || v === t) continue;
      const nS = nodeSize[s], nT = nodeSize[t], nV = nodeSize[v], tot = nS + nT + nV;
      setDist(u, v, Math.max(0, ((nV + nS) / tot) * getDist(s, v) + ((nV + nT) / tot) * getDist(t, v) - (nV / tot) * dST));
    }

    active[s] = 0; active[t] = 0; active[u] = 1; activeCount--;
  }

  // Recover leaf order by traversing the merge tree
  const order = [];
  const stack = [merge.length]; // 1-based root step code
  while (stack.length) {
    const code = stack.pop();
    if (code < 0) { order.push(cells[-code - 1]); continue; }
    const kids = children.get(code);
    if (kids) { stack.push(kids[1]); stack.push(kids[0]); }
  }

  return { merge, height, order };
}

export function generateDemoData(nCells = 200, nBinsPerChr = 80) {
  const chrSizes=[249,243,198,191,181,171,159,145,138,134,135,133,114,107,102,90,83,80,59,64,47,51];
  const chrNames=chrSizes.map((_,i)=>`chr${i+1}`);
  const chrs=[],starts=[],ends=[],startCum=[],endCum=[],chrInfo=[];
  let cum=0;
  chrSizes.forEach((size,ci)=>{const bs=(size*1e6)/nBinsPerChr;const cs=cum;for(let b=0;b<nBinsPerChr;b++){chrs.push(chrNames[ci]);const s=Math.round(b*bs),e=Math.round((b+1)*bs);starts.push(s);ends.push(e);startCum.push(cum+s);endCum.push(cum+e);}cum+=size*1e6;chrInfo.push({chr:chrNames[ci],start_cum:cs,end_cum:cum,mid_cum:(cs+cum)/2});});
  const nBins=chrs.length;
  const cells=Array.from({length:nCells},(_,i)=>`CELL_${String(i+1).padStart(4,"0")}`);
  const profiles={},raw={},quality={},nMajor={},nMinor={},cell_types={};
  const pats=[
    ()=>({t:Array(nBins).fill(2),maj:Array(nBins).fill(1),min:Array(nBins).fill(1)}),
    ()=>{const t=Array(nBins).fill(2),m=Array(nBins).fill(1),n=Array(nBins).fill(1);for(let i=nBinsPerChr*3;i<nBinsPerChr*5;i++){t[i]=3;m[i]=2;}return{t,maj:m,min:n};},
    ()=>{const t=Array(nBins).fill(2),m=Array(nBins).fill(1),n=Array(nBins).fill(1);for(let i=0;i<nBinsPerChr;i++){t[i]=1;m[i]=1;n[i]=0;}for(let i=nBinsPerChr*7;i<nBinsPerChr*9;i++){t[i]=3;m[i]=2;}return{t,maj:m,min:n};},
    ()=>{const t=Array(nBins).fill(2),m=Array(nBins).fill(1),n=Array(nBins).fill(1);for(let i=nBinsPerChr*12;i<nBinsPerChr*14;i++){t[i]=1;m[i]=1;n[i]=0;}for(let i=nBinsPerChr*16;i<nBinsPerChr*18;i++){t[i]=4;m[i]=3;n[i]=1;}return{t,maj:m,min:n};},
    ()=>{const t=Array(nBins).fill(2),m=Array(nBins).fill(1),n=Array(nBins).fill(1);for(let i=nBinsPerChr*2;i<nBinsPerChr*4;i++){t[i]=0;m[i]=0;n[i]=0;}return{t,maj:m,min:n};},
  ];
  cells.forEach((c,ci)=>{
    const p=pats[Math.floor(Math.random()*pats.length)]();
    profiles[c]=p.t;nMajor[c]=p.maj;nMinor[c]=p.min;
    cell_types[c]=DEMO_CELL_TYPES[ci % DEMO_CELL_TYPES.length];
    const noise=Math.random()*0.3+0.1;
    raw[c]=p.t.map(v=>Math.max(0,v+(Math.random()-0.5)*noise*4));
    const res=p.t.map((v,i)=>Math.abs(raw[c][i]-v));
    const mr=d3.median(res);
    const vr=raw[c].filter(v=>v>0);
    const bl=d3.median(vr);
    let mapd=0;
    if(bl>0){
      // logR per bin, preserving position so chromosome boundaries are respected
      const logR=raw[c].map(v=>v>0?Math.log2(v/bl):NaN);
      // within-chromosome absolute diffs only — matching R's mapd(bins$logR, bins$chr)
      const diffs=[];
      for(let i=1;i<nBins;i++){if(chrs[i]===chrs[i-1]&&!isNaN(logR[i])&&!isNaN(logR[i-1]))diffs.push(Math.abs(logR[i]-logR[i-1]));}
      mapd=d3.median(diffs)||0;
    }
    const bins_with_cna = p.t.filter(v => Math.round(v) !== 2).length;
    quality[c]={median_residual:mr,mapd,bins_with_cna};
  });
  const order=[...cells].sort((a,b)=>{let da=0,db=0;for(let i=0;i<nBins;i+=10){da+=profiles[a][i];db+=profiles[b][i];}return da-db;});
  const dendro = buildDemoDendrogram(cells, profiles);
  const dendro_as = buildDemoDendrogram(cells, profiles, nMajor, nMinor);
  return { metadata:{n_cells:nCells,n_bins:nBins,is_allele_specific:true}, bins:{chr:chrs,start:starts,end:ends,start_cum:startCum,end_cum:endCum}, chr_info:chrInfo, profiles, raw, quality, clustering_order:order, nMajor, nMinor, dendrogram: dendro, dendrogram_as: dendro_as, cell_types };
}
