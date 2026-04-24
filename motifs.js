/* Graph / network motif generators.
   Deterministic per-seed so layouts stay stable across renders.
*/
(function () {
  function mulberry32(seed) {
    let t = seed >>> 0;
    return function () {
      t += 0x6D2B79F5;
      let r = Math.imul(t ^ (t >>> 15), 1 | t);
      r ^= r + Math.imul(r ^ (r >>> 7), 61 | r);
      return ((r ^ (r >>> 14)) >>> 0) / 4294967296;
    };
  }

  // Build a dotted-graph SVG. Nodes placed on perturbed lattice;
  // edges between nearest neighbors; a small subset highlighted blue.
  function graphSVG(opts) {
    const {
      width = 900, height = 600, cols = 9, rows = 7,
      jitter = 0.28, nodeR = 2.4, blueR = 3.4, seed = 7,
      connectK = 2, highlightFrac = 0.18, edgeOpacity = 0.25,
      cluster = null,        // {cx, cy, r} — pull nodes toward a center
    } = opts || {};

    const rnd = mulberry32(seed);
    const gx = width / (cols + 1);
    const gy = height / (rows + 1);
    const nodes = [];
    for (let i = 0; i < cols; i++) {
      for (let j = 0; j < rows; j++) {
        let x = gx * (i + 1) + (rnd() - 0.5) * gx * jitter * 2;
        let y = gy * (j + 1) + (rnd() - 0.5) * gy * jitter * 2;
        if (cluster) {
          const dx = cluster.cx - x, dy = cluster.cy - y;
          const d = Math.hypot(dx, dy);
          if (d > cluster.r) {
            const k = (d - cluster.r) / d * 0.45;
            x += dx * k; y += dy * k;
          }
        }
        nodes.push({ x, y, blue: false });
      }
    }
    // Highlight some near the center of mass
    const cxN = nodes.reduce((s,n)=>s+n.x,0)/nodes.length;
    const cyN = nodes.reduce((s,n)=>s+n.y,0)/nodes.length;
    const sorted = [...nodes].sort((a,b)=> (Math.hypot(a.x-cxN,a.y-cyN) - Math.hypot(b.x-cxN,b.y-cyN)));
    const nHi = Math.max(3, Math.floor(nodes.length * highlightFrac));
    for (let i = 0; i < nHi; i++) sorted[i].blue = true;

    // Edges: connect each node to its nearest k neighbors (dedup).
    const edges = new Set();
    const eArr = [];
    for (let i = 0; i < nodes.length; i++) {
      const dists = [];
      for (let j = 0; j < nodes.length; j++) {
        if (i === j) continue;
        dists.push({ j, d: Math.hypot(nodes[i].x-nodes[j].x, nodes[i].y-nodes[j].y) });
      }
      dists.sort((a,b)=>a.d-b.d);
      for (let k = 0; k < connectK; k++) {
        const j = dists[k].j;
        const key = i < j ? i+'-'+j : j+'-'+i;
        if (edges.has(key)) continue;
        edges.add(key);
        const blue = nodes[i].blue && nodes[j].blue;
        eArr.push({ a: i, b: j, blue });
      }
    }

    let svg = `<svg class="motif" viewBox="0 0 ${width} ${height}" preserveAspectRatio="xMidYMid meet" style="width:100%;height:100%;">`;
    // edges
    for (const e of eArr) {
      const a = nodes[e.a], b = nodes[e.b];
      svg += `<line class="${e.blue?'b':''}" x1="${a.x.toFixed(1)}" y1="${a.y.toFixed(1)}" x2="${b.x.toFixed(1)}" y2="${b.y.toFixed(1)}" ${e.blue?'':`opacity="${edgeOpacity}"`}/>`;
    }
    // nodes
    for (const n of nodes) {
      svg += `<circle class="${n.blue?'b':''}" cx="${n.x.toFixed(1)}" cy="${n.y.toFixed(1)}" r="${n.blue?blueR:nodeR}"/>`;
    }
    svg += `</svg>`;
    return svg;
  }

  // Sparse line-of-dots decoration — a thin band of connected nodes.
  function strandSVG(opts) {
    const { width = 900, height = 120, n = 14, seed = 3, nodeR = 2.2, blueEvery = 5 } = opts || {};
    const rnd = mulberry32(seed);
    const pts = [];
    for (let i = 0; i < n; i++) {
      const t = i / (n - 1);
      const x = 20 + t * (width - 40);
      const y = height/2 + (rnd() - 0.5) * height * 0.55;
      pts.push({ x, y, blue: (i % blueEvery) === 2 });
    }
    let svg = `<svg class="motif" viewBox="0 0 ${width} ${height}" preserveAspectRatio="none" style="width:100%;height:100%;">`;
    for (let i = 0; i < pts.length - 1; i++) {
      const a = pts[i], b = pts[i+1];
      const bl = a.blue && b.blue;
      svg += `<line class="${bl?'b':''}" x1="${a.x}" y1="${a.y}" x2="${b.x}" y2="${b.y}" ${bl?'':'opacity="0.3"'}/>`;
    }
    for (const p of pts) svg += `<circle class="${p.blue?'b':''}" cx="${p.x}" cy="${p.y}" r="${nodeR}"/>`;
    svg += `</svg>`;
    return svg;
  }

  // Concentric rings of hits. Graph edges connect each hit to nearby hits
  // on adjacent rings. A handful of radial "real track" chains are
  // highlighted blue — one hit per ring, curving slightly outward.
  function detectorSVG(opts) {
    const {
      size = 780,
      seed = 5,
      rings = 8,          // number of concentric rings
      innerR = 0.08,      // fraction of size
      outerR = 0.48,
      hitsMin = 18,       // hits on innermost ring
      hitsGrow = 4,       // extra hits per ring outward
      nTracks = 5,        // highlighted "real tracks"
      nodeR = 2.6,
      blueR = 4.4,
      edgeOpacity = 0.22,
      angularJitter = 0.55, // irregular spacing along the ring (fraction of step)
      radialJitter = 0.14,  // radial scatter per hit (fraction of gap)
    } = opts || {};

    const cx = size / 2, cy = size / 2;
    const rnd = mulberry32(seed);

    // Build rings of hits — irregular angular + radial spacing.
    const allHits = [];
    const perRing = [];
    const gap = (outerR - innerR) * size / (rings - 1);
    for (let r = 0; r < rings; r++) {
      const t = r / (rings - 1);
      const radius = (innerR + t * (outerR - innerR)) * size;
      const nHits = Math.round(hitsMin + r * hitsGrow);
      const offset = rnd() * Math.PI * 2;
      const step = (Math.PI * 2) / nHits;
      const ringArr = [];
      for (let i = 0; i < nHits; i++) {
        // Non-uniform angular placement
        const a = offset + i * step + (rnd() - 0.5) * step * angularJitter * 2;
        const rr = radius + (rnd() - 0.5) * gap * radialJitter * 2;
        const x = cx + Math.cos(a) * rr;
        const y = cy + Math.sin(a) * rr;
        const h = { x, y, ring: r, a, blue: false };
        ringArr.push(h);
        allHits.push(h);
      }
      perRing.push(ringArr);
    }

    // "Real track" chains: every track starts at the EXACT center
    // (interaction point). Then picks ~every other ring and always includes
    // the outermost ring so tracks reach the edge.
    const trackChains = [];
    const usedAngles = [];
    const ipPoint = { x: cx, y: cy, ring: -1, a: 0, blue: false, isIP: true };
    for (let t = 0; t < nTracks; t++) {
      let a0 = rnd() * Math.PI * 2;
      for (let tries = 0; tries < 20; tries++) {
        let minGap = Math.PI;
        for (const ua of usedAngles) {
          const d = Math.abs(Math.atan2(Math.sin(a0-ua), Math.cos(a0-ua)));
          if (d < minGap) minGap = d;
        }
        if (minGap > 0.55) break;
        a0 = rnd() * Math.PI * 2;
      }
      usedAngles.push(a0);
      const kappa = (rnd() - 0.5) * 0.30;
      const chain = [ipPoint]; // start at interaction point
      const ringStep = 2;
      const ringIdxs = [];
      for (let r = 0; r < rings; r += ringStep) ringIdxs.push(r);
      if (ringIdxs[ringIdxs.length - 1] !== rings - 1) ringIdxs.push(rings - 1);
      for (const r of ringIdxs) {
        const targetA = a0 + kappa * (r / (rings - 1)) * Math.PI;
        let best = null, bestD = Infinity;
        for (const h of perRing[r]) {
          const d = Math.abs(Math.atan2(Math.sin(h.a - targetA), Math.cos(h.a - targetA)));
          if (d < bestD) { bestD = d; best = h; }
        }
        if (best) { best.blue = true; chain.push(best); }
      }
      trackChains.push(chain);
    }

    // Neighbor edges between adjacent rings (combinatorial graph).
    const edges = [];
    for (let r = 0; r < rings - 1; r++) {
      for (const h of perRing[r]) {
        const cands = perRing[r+1].map(n => ({
          n, d: Math.hypot(n.x - h.x, n.y - h.y)
        })).sort((a,b)=>a.d-b.d).slice(0, 2);
        for (const c of cands) edges.push({ a: h, b: c.n });
      }
    }

    // Render — viewBox expanded to give space for chronometer bezel + lugs.
    const pad = size * 0.16; // bezel + lug padding
    const vbSize = size + pad * 2;
    const vbCx = vbSize / 2, vbCy = vbSize / 2;
    const offX = pad, offY = pad;
    let svg = `<svg class="motif detector" viewBox="0 0 ${vbSize} ${vbSize}" preserveAspectRatio="xMidYMid meet" style="width:100%;height:100%;overflow:visible;">`;

    // --- Chronometer bezel (analog, hairline) ---
    // Inline styles are used on bezel elements because CSS `.motif circle { fill: var(--ink) }`
    // overrides presentation attributes. Inline `style="..."` wins.
    const dialR = size * 0.50;                       // inner dial edge
    const bezelOuterR = size * 0.50 + pad * 0.38;    // outer bezel
    const bezelMidR = (dialR + bezelOuterR) / 2;
    const ink = 'var(--ink)';
    // Lugs — slimmer, with a narrower "neck" where they meet the bezel (analog watch style)
    const lugW = pad * 0.42, lugH = pad * 0.78;
    const lugXOff = bezelOuterR * 0.58;
    const lugs = [
      { cx: vbCx - lugXOff, cy: vbCy - bezelOuterR - lugH * 0.32 },
      { cx: vbCx + lugXOff, cy: vbCy - bezelOuterR - lugH * 0.32 },
      { cx: vbCx - lugXOff, cy: vbCy + bezelOuterR + lugH * 0.32 },
      { cx: vbCx + lugXOff, cy: vbCy + bezelOuterR + lugH * 0.32 },
    ];
    for (const lug of lugs) {
      // Trapezoidal lug with rounded top
      const x0 = lug.cx - lugW/2, x1 = lug.cx + lugW/2;
      const y0 = lug.cy - lugH/2, y1 = lug.cy + lugH/2;
      const narrow = lugW * 0.30;
      // Path: rounded outer edge + narrowing neck toward bezel
      const isTop = lug.cy < vbCy;
      if (isTop) {
        svg += `<path d="M ${(x0+narrow*0.3).toFixed(1)} ${y1.toFixed(1)} L ${x0.toFixed(1)} ${(y0+lugW*0.35).toFixed(1)} Q ${x0.toFixed(1)} ${y0.toFixed(1)} ${(x0+lugW*0.35).toFixed(1)} ${y0.toFixed(1)} L ${(x1-lugW*0.35).toFixed(1)} ${y0.toFixed(1)} Q ${x1.toFixed(1)} ${y0.toFixed(1)} ${x1.toFixed(1)} ${(y0+lugW*0.35).toFixed(1)} L ${(x1-narrow*0.3).toFixed(1)} ${y1.toFixed(1)}" style="fill:none;stroke:${ink};stroke-width:1.4;opacity:0.85;"/>`;
      } else {
        svg += `<path d="M ${(x0+narrow*0.3).toFixed(1)} ${y0.toFixed(1)} L ${x0.toFixed(1)} ${(y1-lugW*0.35).toFixed(1)} Q ${x0.toFixed(1)} ${y1.toFixed(1)} ${(x0+lugW*0.35).toFixed(1)} ${y1.toFixed(1)} L ${(x1-lugW*0.35).toFixed(1)} ${y1.toFixed(1)} Q ${x1.toFixed(1)} ${y1.toFixed(1)} ${x1.toFixed(1)} ${(y1-lugW*0.35).toFixed(1)} L ${(x1-narrow*0.3).toFixed(1)} ${y0.toFixed(1)}" style="fill:none;stroke:${ink};stroke-width:1.4;opacity:0.85;"/>`;
      }
    }
    // Crown at 3 o'clock — small fluted cylinder
    const crownW = pad * 0.22, crownH = pad * 0.38;
    const crownX = vbCx + bezelOuterR - 1;
    const crownY = vbCy;
    svg += `<rect x="${crownX.toFixed(1)}" y="${(crownY - crownH/2).toFixed(1)}" width="${crownW.toFixed(1)}" height="${crownH.toFixed(1)}" rx="1.5" style="fill:none;stroke:${ink};stroke-width:1.4;opacity:0.85;"/>`;
    // Crown stem (little nub protruding further)
    const stemW = crownW * 0.6, stemH = crownH * 0.4;
    svg += `<rect x="${(crownX + crownW - 0.5).toFixed(1)}" y="${(crownY - stemH/2).toFixed(1)}" width="${stemW.toFixed(1)}" height="${stemH.toFixed(1)}" rx="1" style="fill:none;stroke:${ink};stroke-width:1.2;opacity:0.85;"/>`;

    // Bezel rings — hairline, no fill. Inline style beats `.motif circle` CSS.
    svg += `<circle cx="${vbCx}" cy="${vbCy}" r="${bezelOuterR.toFixed(1)}" style="fill:none;stroke:${ink};stroke-width:1.8;opacity:0.95;"/>`;
    svg += `<circle cx="${vbCx}" cy="${vbCy}" r="${(bezelOuterR - pad*0.08).toFixed(1)}" style="fill:none;stroke:${ink};stroke-width:0.8;opacity:0.55;"/>`;
    svg += `<circle cx="${vbCx}" cy="${vbCy}" r="${dialR.toFixed(1)}" style="fill:none;stroke:${ink};stroke-width:0.8;opacity:0.55;"/>`;

    // Tick marks — 60 minute ticks, emphasize hour (every 5); blue at 12/3/6/9
    for (let i = 0; i < 60; i++) {
      const ang = (i / 60) * Math.PI * 2 - Math.PI / 2;
      const isHour = i % 5 === 0;
      const isQuarter = i % 15 === 0;
      const tIn  = bezelMidR - (isHour ? 7 : 2.5);
      const tOut = bezelMidR + (isHour ? 7 : 2.5);
      const sw = isHour ? (isQuarter ? 2.4 : 1.6) : 0.8;
      const col = isQuarter ? 'var(--accent)' : ink;
      const op = isHour ? 0.95 : 0.6;
      svg += `<line x1="${(vbCx + Math.cos(ang)*tIn).toFixed(1)}" y1="${(vbCy + Math.sin(ang)*tIn).toFixed(1)}" x2="${(vbCx + Math.cos(ang)*tOut).toFixed(1)}" y2="${(vbCy + Math.sin(ang)*tOut).toFixed(1)}" style="stroke:${col};stroke-width:${sw};stroke-linecap:round;opacity:${op};"/>`;
    }
    // Hour numerals 12/3/6/9
    const numerals = [
      { n: '12', ang: -Math.PI/2 },
      { n: '3',  ang: 0 },
      { n: '6',  ang: Math.PI/2 },
      { n: '9',  ang: Math.PI },
    ];
    const numR = dialR - 24;
    for (const nm of numerals) {
      const nx = vbCx + Math.cos(nm.ang) * numR;
      const ny = vbCy + Math.sin(nm.ang) * numR;
      svg += `<text x="${nx.toFixed(1)}" y="${ny.toFixed(1)}" text-anchor="middle" dominant-baseline="central" font-family="'JetBrains Mono', ui-monospace, monospace" font-size="15" font-weight="500" style="fill:${ink};opacity:0.6;letter-spacing:0.5px;">${nm.n}</text>`;
    }

    // --- Apply dial offset to all hit/track coordinates ---
    // Shift existing geometry into the padded viewBox.
    const shift = (x, y) => [x + offX, y + offY];

    // Dim combinatorial edges (all of them, unmarked)
    for (const e of edges) {
      const [x1, y1] = shift(e.a.x, e.a.y);
      const [x2, y2] = shift(e.b.x, e.b.y);
      svg += `<line x1="${x1.toFixed(1)}" y1="${y1.toFixed(1)}" x2="${x2.toFixed(1)}" y2="${y2.toFixed(1)}" stroke="var(--ink)" stroke-width="0.8" opacity="${edgeOpacity}"/>`;
    }

    // Thick blue track polylines — a single smooth arc per chain so the
    // curve has no inflection points (charged-particle trajectory).
    // Also collect extra hits sampled along each arc so tracks aren't just
    // 3-hit skeletons.
    const extraTrackHits = [];
    for (const chain of trackChains) {
      if (chain.length < 2) continue;
      const p0 = chain[0];
      const pN = chain[chain.length - 1];
      const pM = chain[Math.floor(chain.length / 2)];
      const ax = p0.x + offX, ay = p0.y + offY;
      const bx = pM.x + offX, by = pM.y + offY;
      const cx2 = pN.x + offX, cy2 = pN.y + offY;
      const d = 2 * (ax * (by - cy2) + bx * (cy2 - ay) + cx2 * (ay - by));
      let dStr;
      let arcInfo = null;
      // Extend endpoint slightly beyond outermost hit (outbound particle feel).
      const vxEnd = pN.x - cx, vyEnd = pN.y - cy;
      const extendFrac = 0.10;
      const pEx = { x: cx + vxEnd * (1 + extendFrac) + offX, y: cy + vyEnd * (1 + extendFrac) + offY };
      if (Math.abs(d) < 1e-3) {
        dStr = `M ${ax.toFixed(1)} ${ay.toFixed(1)} L ${pEx.x.toFixed(1)} ${pEx.y.toFixed(1)}`;
      } else {
        const ux = ((ax*ax + ay*ay) * (by - cy2) + (bx*bx + by*by) * (cy2 - ay) + (cx2*cx2 + cy2*cy2) * (ay - by)) / d;
        const uy = ((ax*ax + ay*ay) * (cx2 - bx) + (bx*bx + by*by) * (ax - cx2) + (cx2*cx2 + cy2*cy2) * (bx - ax)) / d;
        const r = Math.hypot(ax - ux, ay - uy);
        const cross = (bx - ax) * (cy2 - ay) - (by - ay) * (cx2 - ax);
        const sweep = cross > 0 ? 1 : 0;
        const angEx0 = Math.atan2(cy2 - uy, cx2 - ux);
        const angEx1 = Math.atan2(pEx.y - uy, pEx.x - ux);
        let dAng = angEx1 - angEx0;
        if (sweep === 1 && dAng < 0) dAng += Math.PI * 2;
        if (sweep === 0 && dAng > 0) dAng -= Math.PI * 2;
        const maxExt = 0.25;
        if (Math.abs(dAng) > maxExt) dAng = Math.sign(dAng) * maxExt;
        const angExt = angEx0 + dAng;
        const exX = ux + Math.cos(angExt) * r;
        const exY = uy + Math.sin(angExt) * r;
        dStr = `M ${ax.toFixed(1)} ${ay.toFixed(1)} A ${r.toFixed(1)} ${r.toFixed(1)} 0 0 ${sweep} ${exX.toFixed(1)} ${exY.toFixed(1)}`;
        arcInfo = { ux, uy, r, sweep };
      }
      svg += `<path d="${dStr}" fill="none" stroke="var(--accent)" stroke-width="3.2" stroke-linecap="round" stroke-linejoin="round" opacity="0.95"/>`;

      // Sample extra hits between consecutive chain hits along the arc.
      // Skip the interaction-point segment (don't draw a hit right at the center).
      for (let i = 0; i < chain.length - 1; i++) {
        const a = chain[i], b = chain[i+1];
        if (a.isIP) continue; // no extra hit in the IP→first-ring segment
        const t = 0.5;
        let x, y;
        const aShift = { x: a.x + offX, y: a.y + offY };
        const bShift = { x: b.x + offX, y: b.y + offY };
        if (arcInfo) {
          const ang0 = Math.atan2(aShift.y - arcInfo.uy, aShift.x - arcInfo.ux);
          const ang1 = Math.atan2(bShift.y - arcInfo.uy, bShift.x - arcInfo.ux);
          let delta = ang1 - ang0;
          if (arcInfo.sweep === 1 && delta < 0) delta += Math.PI * 2;
          if (arcInfo.sweep === 0 && delta > 0) delta -= Math.PI * 2;
          const ang = ang0 + delta * t;
          const jr = arcInfo.r + (rnd() - 0.5) * 3;
          x = arcInfo.ux + Math.cos(ang) * jr;
          y = arcInfo.uy + Math.sin(ang) * jr;
        } else {
          x = aShift.x + (bShift.x - aShift.x) * t + (rnd() - 0.5) * 3;
          y = aShift.y + (bShift.y - aShift.y) * t + (rnd() - 0.5) * 3;
        }
        extraTrackHits.push({ x, y });
      }
    }

    // All hits on top. Highlighted ones bigger/blue.
    for (const h of allHits) {
      if (h.blue) continue;
      svg += `<circle cx="${(h.x+offX).toFixed(1)}" cy="${(h.y+offY).toFixed(1)}" r="${nodeR}" fill="var(--ink)"/>`;
    }
    for (const h of allHits) {
      if (!h.blue) continue;
      svg += `<circle cx="${(h.x+offX).toFixed(1)}" cy="${(h.y+offY).toFixed(1)}" r="${blueR}" fill="var(--accent)"/>`;
    }
    // Extra sampled hits along fitted arcs
    for (const h of extraTrackHits) {
      svg += `<circle cx="${h.x.toFixed(1)}" cy="${h.y.toFixed(1)}" r="${(blueR*0.78).toFixed(1)}" fill="var(--accent)"/>`;
    }
    // Interaction point marker (small accent circle + ring)
    svg += `<circle cx="${vbCx}" cy="${vbCy}" r="${(blueR*1.5).toFixed(1)}" fill="none" stroke="var(--accent)" stroke-width="1.5" opacity="0.8"/>`;
    svg += `<circle cx="${vbCx}" cy="${vbCy}" r="${(blueR*0.65).toFixed(1)}" fill="var(--accent)"/>`;

    svg += `</svg>`;
    return svg;
  }

  window.Motifs = { graphSVG, strandSVG, detectorSVG };
})();
