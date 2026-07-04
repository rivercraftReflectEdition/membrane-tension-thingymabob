'use strict';
/* =========================================================================
 * App shell. The physics core above is solved once (dimensionless) in a
 * worker; every control change afterwards is a pure rescale:
 *   w = (σ·g·L³/F)·ŵ     slope = (σ·g·L²/F)·|∇ŵ|
 * ========================================================================= */
(function () {
const CFG = { nn: 41, nu: 0.34, maxIter: 30, floorFrac: 0.02 };
const P_SUN = 9.1e-6;                      // Pa, solar radiation pressure — the only static load in space
const SLEW_ACC = 0.65;                     // peak alpha ≈ 0.65·omega² for an overhead tracking pass
const DEFAULTS = { F: 250, g: 9.81, L: 15, massG: 350, tgt: 2.3, E_GPa: 4.4, rho: 1430,
                   phiDeg: 10, slew: 0.83 };
const state = { ...DEFAULTS };
/* shareable state via query string, e.g. ?g=1.62&F=400 (no storage APIs used) */
const QMAP = { F: 'F', g: 'g', L: 'L', m: 'massG', tgt: 'tgt', E: 'E_GPa', rho: 'rho',
               phi: 'phiDeg', slew: 'slew' };
const QLIM = { F: [0, 1000], g: [0, 9.81], L: [2, 20], massG: [100, 1000], tgt: [0.1, 20],
               E_GPa: [0.1, 500], rho: [500, 5000], phiDeg: [0, 45], slew: [0, 5] };
{
  const qs = new URLSearchParams(location.search);
  for (const [q, key] of Object.entries(QMAP)) {
    const v = parseFloat(qs.get(q));
    if (isFinite(v)) state[key] = Math.min(QLIM[key][1], Math.max(QLIM[key][0], v));
  }
}
const reduceMotion = matchMedia('(prefers-reduced-motion: reduce)').matches;

/* ------------------------------- theme -----------------------------------
 * Orange / white / black. Both palettes follow the same roles; canvases read
 * the active one at draw time. No storage — default follows the OS, the
 * header button flips it for the session. */
const THEMES = {
  light: {
    ink: '#0b0b0b', ink2: '#52514e', muted: '#898781', grid: '#e1e0d9',
    axis: '#c3c2b7', surface: '#fcfcfb', accent: '#e05500', space: '#0b0b0b',
    inkRGB: '11,11,11',
    ramp: ['#ffe9d6', '#ffcf9e', '#fca55c', '#e87613', '#b25000', '#6e3000'],
  },
  dark: {
    ink: '#f5f4f0', ink2: '#c3c2b7', muted: '#8f8d86', grid: '#2c2c2a',
    axis: '#4a4a46', surface: '#1a1a19', accent: '#d95926', space: '#f5f4f0',
    inkRGB: '245,244,240',
    ramp: ['#2b1c10', '#54300e', '#8a4a0a', '#c66a14', '#f68d33', '#ffbe8a'],
  },
};
for (const t of Object.values(THEMES))
  t.rampRGB = t.ramp.map(hx => [parseInt(hx.slice(1, 3), 16), parseInt(hx.slice(3, 5), 16),
                                parseInt(hx.slice(5, 7), 16)]);
let TH = THEMES.light;
const inkA = a => 'rgba(' + TH.inkRGB + ',' + a + ')';
function setTheme(name) {
  document.documentElement.dataset.theme = name;
  TH = THEMES[name] || THEMES.light;
  const btn = $('themeBtn');
  if (btn) {
    btn.textContent = name === 'dark' ? 'light' : 'dark';
    btn.setAttribute('aria-label', 'Switch to ' + (name === 'dark' ? 'light' : 'dark') + ' theme');
  }
  if (typeof drawAll === 'function') drawAll();
}

const $ = id => document.getElementById(id);
const els = {
  chips: { t: $('chipT'), sig: $('chipSig'), stress: $('chipStress') },
  g: {
    slope: $('stGslope'), slopeI: $('stGslopeI'), sag: $('stGsag'), sagI: $('stGsagI'),
    ten: $('stGten'), note: $('stGnote'), cv: $('cvG'), lbl: $('lblG'),
  },
  s: { slope: $('stSslope'), sag: $('stSsag'), ten: $('stSten'), note: $('stSnote'), cv: $('cvS') },
  reqLine: $('reqLine'), reqFlag: $('reqFlag'),
  cvSec: $('cvSec'), cvChart: $('cvChart'), cvSlew: $('cvSlew'), cvSlewAnim: $('cvSlewAnim'),
  phiStatus: $('phiStatus'), slewOut: $('slewOut'),
};

/* ------------------------------ formatting ------------------------------ */
const fmt = (v, d) => v.toLocaleString('en-US', { minimumFractionDigits: d, maximumFractionDigits: d });
function fmtAuto(v, d2 = 10, d1 = 100) {   // 2 dp, then 1, then 0 as magnitude grows
  if (!isFinite(v)) return '—';
  return fmt(v, v < d2 ? 2 : v < d1 ? 1 : 0);
}
function fmtSmart(v) {                     // like fmtAuto, but tiny values keep 2 sig figs
  if (!isFinite(v)) return '—';
  return v !== 0 && Math.abs(v) < 0.01 ? Number(v.toPrecision(2)).toString() : fmtAuto(v);
}
const unit = (s) => ` <small>${s}</small>`;

/* ------------------------------ physics glue ---------------------------- */
function derived() {
  const L = state.L, mass = state.massG / 1000, F = state.F;
  const g = state.g, E = state.E_GPa * 1e9, rho = state.rho;
  const sigma = analytic.sigma(mass, L);
  const t = analytic.thick(mass, rho, L);
  const N = analytic.tension(F, L);                    // reference convention F/(√2 L)
  const Feff = Math.max(F, 1);                         // slack guard for 1/F scales
  const p = analytic.pressure(sigma, g);
  return {
    L, mass, F, g, E, rho, sigma, t, N, p,
    sScale: sigma * g * L * L / Feff,                  // rad per unit slope coeff (gravity)
    wScale: sigma * g * L * L * L / Feff,              // m per unit ŵ (gravity)
    sScaleSun: P_SUN * L * L / Feff,                   // same, solar pressure only (space)
    wScaleSun: P_SUN * L * L * L / Feff,
    slack: F < 2,
  };
}

/* shapes solved once in the worker */
let shapes = null;                                     // {whatT, slopesT, whatR, slopesR, uhat, wstate, thetaW, coeffs, mesh:{nn,h,ne}}
const WC_ANALYTIC = Math.SQRT2 / (4 * Math.PI);        // ŵ_center of the analytic (circle-equiv) reference

/* ------------------------------ worker ----------------------------------
 * The dimensionless solve depends on ONE control: the catenary half-angle.
 * Solves are cached per whole degree; the worker stays alive and always
 * serves the latest requested angle (intermediate drag values are skipped). */
function packShapes(sh) {
  return {
    whatT: sh.truth.what, slopesT: sh.truth.slopes, whatS: sh.slew.what,
    whatR: sh.ref.what, slopesR: sh.ref.slopes,
    uhat: sh.inPlane.uhat, wstate: sh.inPlane.state, thetaW: sh.inPlane.thetaW,
    coeffs: {
      ST: sh.truth.rmsSlopeCoeff, PVT: sh.truth.pvCoeff, WcT: sh.truth.centerCoeff,
      SS: sh.slew.rmsSlopeCoeff, PVS: sh.slew.pvCoeff,
      SR: sh.ref.rmsSlopeCoeff, PVR: sh.ref.pvCoeff, WcR: sh.ref.centerCoeff,
    },
    mesh: { nn: sh.mesh.nn, h: sh.mesh.h, ne: sh.mesh.ne },
  };
}
const solver = { cache: new Map(), busy: false, queued: null, worker: null };
const phiKey = () => Math.round(state.phiDeg);
function requestShapes() {
  const k = phiKey();
  if (solver.cache.has(k)) {
    if (shapes !== solver.cache.get(k)) { applyShapes(solver.cache.get(k), false); }
    els.phiStatus.textContent = '';
    return;
  }
  els.phiStatus.textContent = '· solving…';
  solver.queued = k;
  pumpSolver();
}
function pumpSolver() {
  if (solver.busy || solver.queued === null) return;
  const k = solver.queued; solver.queued = null; solver.busy = true;
  const msg = { ...CFG, phi: k * Math.PI / 180, phiKey: k };
  if (solver.worker) solver.worker.postMessage(msg);
  else setTimeout(() => {                                // no-worker fallback
    const sh = buildShapes(msg.nn, msg.nu, msg);
    const p = packShapes(sh); p.phiKey = k;
    onSolved(p);
  }, 30);
}
function onSolved(p) {
  solver.busy = false;
  solver.cache.set(p.phiKey, p);
  if (p.phiKey === phiKey() && solver.queued === null) {
    els.phiStatus.textContent = '';
    applyShapes(p, !shapes);
  }
  pumpSolver();
}
function applyShapes(p, first) {
  shapes = p;
  /* per-element |grad| of the slew billow shape, for the maneuver preview */
  shapes.slopesS = elemSlopesOf(p.whatS, p.mesh);
  shapes.maxAbsS = p.whatS.reduce((a, v) => Math.max(a, Math.abs(v)), 0);
  if (!reduceMotion) { anim.gAmp.v = first ? 0 : Math.min(anim.gAmp.v, 0.5); }
  anim.gAmp.t = 1;
  update();
  scheduleTick();
  scheduleSlewAnim();
}
function elemSlopesOf(what, mesh) {
  const { nn, h, ne } = mesh, out = new Float64Array(ne);
  for (let e = 0; e < ne; e++) {
    const i = e % (nn - 1), j = (e / (nn - 1)) | 0, n0 = j * nn + i;
    const wxm = (what[n0 + 1] - what[n0] + what[n0 + nn + 1] - what[n0 + nn]) / (2 * h);
    const wym = (what[n0 + nn] - what[n0] + what[n0 + nn + 1] - what[n0 + 1]) / (2 * h);
    out[e] = Math.hypot(wxm, wym);
  }
  return out;
}
(function startSolver() {
  try {
    const src = document.getElementById('physics-core').textContent +
      ';self.onmessage=function(e){var sh=buildShapes(e.data.nn,e.data.nu,e.data);' +
      'var p=(' + packShapes.toString() + ')(sh);p.phiKey=e.data.phiKey;' +
      'postMessage(p,[p.whatT.buffer,p.slopesT.buffer,p.whatS.buffer,p.whatR.buffer,' +
      'p.slopesR.buffer,p.uhat.buffer,p.wstate.buffer,p.thetaW.buffer]);};';
    solver.worker = new Worker(URL.createObjectURL(new Blob([src], { type: 'text/javascript' })));
    solver.worker.onmessage = e => onSolved(e.data);
    /* worker dead or silent: fall back to solving on the main thread,
     * re-queuing the angle we were waiting for */
    const failover = () => {
      solver.worker = null; solver.busy = false;
      if (solver.queued === null) solver.queued = phiKey();
      pumpSolver();
    };
    solver.worker.onerror = failover;
    setTimeout(() => { if (!shapes && solver.busy) failover(); }, 6000);
  } catch (_) { solver.worker = null; }
  requestShapes();
})();

/* ------------------------------ animation ------------------------------- */
/* the membranes re-settle smoothly: sag amplitudes ease toward their target */
const anim = { gAmp: { v: 1, t: 1 }, raf: 0 };          // multiplier on ground sag
function tick() {
  anim.raf = 0;
  const a = anim.gAmp, k = 0.22;
  a.v += (a.t - a.v) * k;
  if (Math.abs(a.t - a.v) < 0.004) a.v = a.t; else scheduleTick();
  drawAll();
}
function scheduleTick() { if (!anim.raf) anim.raf = requestAnimationFrame(tick); }

/* ------------------------------ update ---------------------------------- */
function update() {
  const d = derived();

  /* chips */
  els.chips.t.textContent = fmt(d.t * 1e6, 2) + ' µm';
  els.chips.sig.textContent = fmt(d.sigma * 1e3, 2) + ' g/m²';
  els.chips.stress.textContent = d.t > 0 ? fmtAuto(d.N / d.t / 1e6) + ' MPa' : '—';

  /* per-environment readouts (ground = current g, space = 0) */
  els.g.lbl.textContent = 'g = ' + fmt(d.g, 2) + ' m/s²';
  const idealSlope = analytic.rmsSlope(d.p, d.L, Math.max(d.N, 1e-9)) * 1e3;
  const idealSag = analytic.sagC(d.p, d.L, Math.max(d.N, 1e-9)) * 1e3;
  const truthSlope = shapes ? shapes.coeffs.ST * d.sScale * 1e3 : NaN;
  const truthSag = shapes ? shapes.coeffs.PVT * d.wScale * 1e3 : NaN;

  const regime = truthSlope > 100 || idealSlope > 100;   // small-slope validity flag
  if (d.slack) {
    els.g.slope.innerHTML = 'slack'; els.g.sag.innerHTML = '—';
    els.g.slopeI.textContent = 'membrane untensioned'; els.g.sagI.textContent = ' ';
  } else {
    els.g.slope.innerHTML = (shapes ? fmtAuto(truthSlope) : '…') + unit('mrad') + (regime ? '†' : '');
    els.g.sag.innerHTML = (shapes ? fmtAuto(truthSag) : '…') + unit('mm');
    els.g.slopeI.textContent = 'uniform ideal ' + fmtAuto(idealSlope) + ' mrad';
    els.g.sagI.textContent = 'uniform ideal ' + fmtAuto(idealSag) + ' mm';
  }
  els.g.ten.innerHTML = fmtAuto(d.N) + unit('N/m');
  els.g.note.textContent = regime && !d.slack ? '† beyond small-slope regime' : ' ';

  /* space: solar radiation pressure is the only transverse load */
  const spSlope = shapes ? shapes.coeffs.ST * d.sScaleSun * 1e3 : NaN;
  const spSag = shapes ? shapes.coeffs.PVT * d.wScaleSun * 1e3 : NaN;
  if (d.slack) {
    els.s.slope.innerHTML = 'slack'; els.s.sag.innerHTML = '—';
  } else {
    els.s.slope.innerHTML = (shapes ? fmtSmart(spSlope) : '…') + unit('mrad');
    els.s.sag.innerHTML = (shapes ? fmtSmart(spSag) : '…') + unit('mm');
  }
  els.s.ten.innerHTML = fmtAuto(d.N) + unit('N/m');

  /* shared required-force line: same target, both loads */
  const tgt = state.tgt / 1e3;
  const FreqS = P_SUN * d.L * d.L / (2 * Math.sqrt(Math.PI) * tgt);   // sunlight only
  if (d.g < 0.005) {
    els.reqLine.textContent = '0 g — gravity gone: sunlight alone needs just ' +
      fmtSmart(FreqS) + ' N/corner for ' + fmt(state.tgt, 2) + ' mrad.';
    els.reqFlag.hidden = true;
  } else {
    const Freq = analytic.Freq(d.sigma, d.g, d.L, tgt);
    els.reqLine.textContent = 'To reach ' + fmt(state.tgt, 2) + ' mrad — ground ' +
      fmt(Freq, 0) + ' N/corner · space ' + fmtSmart(FreqS) +
      ' N (sunlight) · ×' + fmt(Freq / FreqS, 0);
    els.reqFlag.hidden = Freq <= 1000;
  }

  /* slew readout */
  if (shapes && !d.slack) {
    const alpha = SLEW_ACC * Math.pow(state.slew * Math.PI / 180, 2) * 180 / Math.PI;
    els.slewOut.textContent = 'peak α ' + fmtSmart(alpha) + ' °/s² · billow ' +
      fmtSmart(slewSlopeMrad(d, state.slew)) + ' mrad · with sunlight ' +
      fmtSmart(spaceTotalMrad(d, state.slew)) + ' mrad';
  } else els.slewOut.textContent = d.slack ? 'membrane slack' : '…';

  /* aria summaries */
  els.g.cv.setAttribute('aria-label', 'Ground membrane: RMS slope ' +
    (shapes ? fmtAuto(truthSlope) : '—') + ' milliradians, sag ' +
    (shapes ? fmtAuto(truthSag) : '—') + ' millimeters');
  els.s.cv.setAttribute('aria-label', 'Space membrane: essentially flat, RMS slope ' +
    (shapes ? fmtSmart(spSlope) : '—') + ' milliradians under sunlight alone');

  drawAll();
}

/* ------------------------------ canvases -------------------------------- */
function fitCanvas(cv) {
  const r = cv.getBoundingClientRect();
  const dpr = devicePixelRatio || 1;
  const w = Math.round(r.width * dpr), h = Math.round(r.height * dpr);
  if (cv.width !== w || cv.height !== h) { cv.width = w; cv.height = h; }
  const ctx = cv.getContext('2d');
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
  return { ctx, w: r.width, h: r.height };
}

/* sequential single-hue orange ramp: slope 0 → 2×target mrad (theme-aware) */
function rampColor(x) {                                  // x in [0,1]
  const R = TH.rampRGB;
  const t = Math.max(0, Math.min(0.9999, x)) * (R.length - 1);
  const i = Math.floor(t), f = t - i, a = R[i], b = R[i + 1];
  return 'rgb(' + Math.round(a[0] + (b[0] - a[0]) * f) + ',' +
    Math.round(a[1] + (b[1] - a[1]) * f) + ',' + Math.round(a[2] + (b[2] - a[2]) * f) + ')';
}

/* isometric projection constants — low-ish camera so the bow reads */
const KX = 0.866, KY = 0.34, KZ = 1.0, SAGFRAC = 0.24;
/* ONE fixed vertical exaggeration: the drawn sag varies continuously with the
 * controls (a taut mirror looks flat, a slack one droops), no rescale jumps.
 * Beyond ~SAGFRAC·L the drawn amplitude compresses smoothly (tanh), shape
 * preserved, so extreme slack states stay inside the panel. */
const EXZ = 50;
function sagCompression(Apeak, L) {                      // global multiplier <= 1
  const cap = SAGFRAC * L;
  if (!(Apeak > 0) || Apeak <= cap) return 1;
  return cap * Math.tanh(Apeak / cap) / Apeak;
}

function drawMembrane(cv, d, gEnv, amp) {
  const { ctx, w, h } = fitCanvas(cv);
  ctx.clearRect(0, 0, w, h);
  const L = d.L;
  /* fixed world box: width 2L·KX, height 2L·KY + allowance for sag */
  const pad = 10;
  const sagAllow = SAGFRAC * L * KZ * 1.1;
  const sc = Math.min((w - 2 * pad) / (2 * L * KX), (h - 2 * pad - 14) / (2 * L * KY + sagAllow));
  const cx = w / 2, cy = pad + L * KY * sc;
  const proj = (X, Y, z) => [cx + (X - Y) * KX * sc, cy + (X + Y) * KY * sc + z * KZ * sc];

  if (!shapes || d.slack) {                              // placeholder / slack flat sheet
    ctx.strokeStyle = inkA(.25); ctx.lineWidth = 1;
    if (d.slack) ctx.setLineDash([4, 4]);
    const c = [proj(-L/2,-L/2,0), proj(L/2,-L/2,0), proj(L/2,L/2,0), proj(-L/2,L/2,0)];
    ctx.beginPath(); c.forEach((p,i)=> i?ctx.lineTo(p[0],p[1]):ctx.moveTo(p[0],p[1])); ctx.closePath(); ctx.stroke();
    ctx.setLineDash([]);
    if (d.slack) {
      ctx.fillStyle = TH.muted; ctx.textAlign = 'center';
      ctx.font = '10px ' + getComputedStyle(document.body).getPropertyValue('--mono');
      ctx.fillText('untensioned — film hangs slack (out of model)', cx, cy + 4);
    }
    return;
  }

  const nn = shapes.mesh.nn, hh = shapes.mesh.h;
  const what = shapes.whatT, uhat = shapes.uhat;
  const uSc = d.t > 0 ? d.F * (1 - CFG.nu * CFG.nu) / (d.E * d.t) : 0;  // in-plane, TRUE scale
  /* ground: gravity sag; space: solar-pressure sag (µm — reads flat, honestly) */
  const wScaleEnv = gEnv ? d.wScale : d.wScaleSun;
  const sScaleEnv = gEnv ? d.sScale : d.sScaleSun;
  const compress = sagCompression(shapes.coeffs.PVT * wScaleEnv * amp * EXZ, L);
  const wSc = wScaleEnv * amp * EXZ * compress;          // meters of drawn sag per unit ŵ
  const slSc = sScaleEnv * 1e3;                          // mrad per unit slope coeff
  const slMax = 2 * state.tgt;

  /* node screen positions. In-plane displacement u = û·F(1−ν²)/(E·t) is
   * true scale. With catenary cords the edges are cut as arcs — the physics
   * runs on the square domain (few-% sagitta), but the outline is drawn
   * scalloped: transfinite warp, sagitta L·tan(φ/2)/2, zero at the corners. */
  const phiDraw = shapes.phiKey ?? state.phiDeg;         // draw the SOLVED geometry
  const sagit = phiDraw > 0 ? Math.tan(phiDraw * Math.PI / 360) / 2 * L : 0;
  const px = new Float64Array(nn * nn), py = new Float64Array(nn * nn);
  const wx = new Float64Array(nn * nn), wy = new Float64Array(nn * nn);
  for (let j = 0; j < nn; j++) for (let i = 0; i < nn; i++) {
    const n = j * nn + i;
    const xi = i * hh, eta = j * hh;
    let X = (i * hh - 0.5) * L + uhat[2 * n] * uSc;
    let Y = (j * hh - 0.5) * L + uhat[2 * n + 1] * uSc;
    if (sagit > 0) {
      Y += 4 * xi * (1 - xi) * sagit * (1 - 2 * eta);
      X += 4 * eta * (1 - eta) * sagit * (1 - 2 * xi);
    }
    wx[n] = X; wy[n] = Y;
    const p = proj(X, Y, what[n] * wSc);
    px[n] = p[0]; py[n] = p[1];
  }

  /* corner-plane frame (z = 0), dashed hairline behind the surface */
  ctx.save();
  ctx.strokeStyle = TH.axis; ctx.lineWidth = 1; ctx.setLineDash([3, 3]);
  const fr = [proj(-L/2,-L/2,0), proj(L/2,-L/2,0), proj(L/2,L/2,0), proj(-L/2,L/2,0)];
  ctx.beginPath(); fr.forEach((p,i)=> i?ctx.lineTo(p[0],p[1]):ctx.moveTo(p[0],p[1])); ctx.closePath(); ctx.stroke();
  ctx.restore();

  /* surface quads back-to-front, shaded by local slope */
  for (let s = 0; s <= 2 * (nn - 2); s++) {              // anti-diagonals: i+j = s
    for (let i = Math.max(0, s - (nn - 2)); i <= Math.min(nn - 2, s); i++) {
      const j = s - i;
      const e = j * (nn - 1) + i;
      const n0 = j * nn + i, n1 = n0 + 1, n2 = n0 + nn + 1, n3 = n0 + nn;
      const mr = shapes.slopesT[e] * slSc;
      ctx.fillStyle = rampColor(mr / slMax);
      ctx.beginPath();
      ctx.moveTo(px[n0], py[n0]); ctx.lineTo(px[n1], py[n1]);
      ctx.lineTo(px[n2], py[n2]); ctx.lineTo(px[n3], py[n3]);
      ctx.closePath(); ctx.fill();
      ctx.strokeStyle = ctx.fillStyle; ctx.lineWidth = 0.6; ctx.stroke(); // hide AA seams
    }
  }

  /* quiet wireframe every 5th grid line */
  ctx.strokeStyle = inkA(.10); ctx.lineWidth = 0.75;
  for (let j = 0; j < nn; j += 5) {
    ctx.beginPath();
    for (let i = 0; i < nn; i++) { const n = j * nn + i; i ? ctx.lineTo(px[n], py[n]) : ctx.moveTo(px[n], py[n]); }
    ctx.stroke();
    ctx.beginPath();
    for (let i = 0; i < nn; i++) { const n = i * nn + j; i ? ctx.lineTo(px[n], py[n]) : ctx.moveTo(px[n], py[n]); }
    ctx.stroke();
  }

  /* wrinkle hints: short strokes along the wrinkle direction (in-plane state —
   * present in both environments; it is a tension-field property, not gravity) */
  ctx.strokeStyle = inkA(.14); ctx.lineWidth = 0.75;
  const ne1 = nn - 1;
  for (let e = 0; e < shapes.mesh.ne; e++) {
    if (shapes.wstate[e] !== 1) continue;
    const i = e % ne1, j = (e / ne1) | 0;
    if (((i + 2 * j) % 3)) continue;                     // thin them out
    const n0 = j * nn + i, n2 = n0 + nn + 1;
    const mx = (px[n0] + px[n2]) / 2, my = (py[n0] + py[n2]) / 2;
    /* project the in-plane wrinkle direction */
    const th = shapes.thetaW[e], dx = Math.cos(th), dy = Math.sin(th);
    const sx = (dx - dy) * KX, sy = (dx + dy) * KY;
    const len = 0.55 * hh * L * sc / Math.max(Math.hypot(sx, sy), 1e-9);
    ctx.beginPath();
    ctx.moveTo(mx - sx * len / 2, my - sy * len / 2);
    ctx.lineTo(mx + sx * len / 2, my + sy * len / 2);
    ctx.stroke();
  }

  /* membrane edge outline — at φ>0 this is the catenary cord itself */
  ctx.strokeStyle = sagit > 0 ? inkA(.55) : inkA(.4);
  ctx.lineWidth = sagit > 0 ? 1.4 : 1;
  ctx.beginPath();
  const edge = [];
  for (let i = 0; i < nn; i++) edge.push(i);                       // bottom
  for (let j = 1; j < nn; j++) edge.push(j * nn + nn - 1);         // right
  for (let i = nn - 2; i >= 0; i--) edge.push((nn - 1) * nn + i);  // top
  for (let j = nn - 2; j >= 1; j--) edge.push(j * nn);             // left
  edge.forEach((n, k) => k ? ctx.lineTo(px[n], py[n]) : ctx.moveTo(px[n], py[n]));
  ctx.closePath(); ctx.stroke();

  /* corner-plane frame again in front, plus droop ticks at the mid-edges —
   * the gap between frame and film is what gravity costs */
  ctx.strokeStyle = inkA(.38); ctx.lineWidth = 1; ctx.setLineDash([3, 3]);
  ctx.beginPath(); fr.forEach((p, i) => i ? ctx.lineTo(p[0], p[1]) : ctx.moveTo(p[0], p[1]));
  ctx.closePath(); ctx.stroke();
  ctx.setLineDash([]);
  if (gEnv && wSc > 0) {
    ctx.strokeStyle = inkA(.30);
    const mid = (nn - 1) / 2;
    const midNodes = [mid, mid * nn, mid * nn + nn - 1, (nn - 1) * nn + mid];
    for (const n of midNodes) {                        // vertical drop, z = 0 → film
      const p0 = proj(wx[n], wy[n], 0);
      ctx.beginPath(); ctx.moveTo(p0[0], p0[1]);
      ctx.lineTo(px[n], py[n]); ctx.stroke();
    }
  }

  /* corner tabs + boom pull stubs */
  ctx.fillStyle = TH.ink; ctx.strokeStyle = TH.ink2; ctx.lineWidth = 1;
  const corners = [[0, 0, -1, -1], [nn - 1, 0, 1, -1], [0, nn - 1, -1, 1], [nn - 1, nn - 1, 1, 1]];
  for (const [ci, cj, sx, sy] of corners) {
    const n = cj * nn + ci;
    const bx = (sx - sy) * KX, by = (sx + sy) * KY;
    const bl = 0.05 * L * sc / Math.hypot(bx, by);
    ctx.beginPath(); ctx.moveTo(px[n], py[n]);
    ctx.lineTo(px[n] + bx * bl, py[n] + by * bl); ctx.stroke();
    ctx.beginPath(); ctx.arc(px[n] + bx * bl, py[n] + by * bl, 1.8, 0, 7); ctx.fill();
    ctx.fillRect(px[n] - 2, py[n] - 2, 4, 4);
  }

  /* annotation */
  ctx.fillStyle = TH.muted; ctx.font = '10px ' + getComputedStyle(document.body).getPropertyValue('--mono');
  ctx.textAlign = 'right';
  ctx.fillText((!gEnv ? 'sunlight only · ' : d.g < 1e-4 ? 'no gravity load · ' : '') +
    'vertical sag ×' + fmt(EXZ, 0) + (compress < 0.98 ? ' (compressed)' : '') +
    ' · in-plane true scale', w - 4, h - 4);
}

/* ------------------------- mid-span section chart ----------------------- */
const hover = { sec: null, chart: null, slew: null };
function drawSection() {
  const { ctx, w, h } = fitCanvas(els.cvSec);
  ctx.clearRect(0, 0, w, h);
  const d = derived();
  const mono = '10px ' + getComputedStyle(document.body).getPropertyValue('--mono');
  const m = { l: 34, r: 20, t: 16, b: 22 };
  const pw = w - m.l - m.r, ph = h - m.t - m.b;
  if (pw < 40 || ph < 30) return;

  /* series (mm, sag positive downward) along the mid row */
  let truth = null, ideal = null;
  if (shapes && !d.slack) {
    const nn = shapes.mesh.nn, mid = (nn - 1) / 2;
    truth = []; ideal = [];
    const refNorm = WC_ANALYTIC / Math.max(shapes.coeffs.WcR, 1e-12); // pin ref curve to the analytic w꜀ formula
    for (let i = 0; i < nn; i++) {
      truth.push(shapes.whatT[mid * nn + i] * d.wScale * 1e3 * anim.gAmp.v);
      ideal.push(shapes.whatR[mid * nn + i] * d.wScale * 1e3 * refNorm * anim.gAmp.v);
    }
  }
  const ymax = niceCeil(Math.max(1, truth ? Math.max(...truth, ...ideal) : 1));
  const X = i => m.l + (i / (CFG.nn - 1)) * pw;
  const Y = v => m.t + (v / ymax) * ph;

  /* frame + gridlines + ticks (sag in mm, downward) */
  ctx.strokeStyle = TH.grid; ctx.lineWidth = 1;
  ctx.font = mono; ctx.fillStyle = TH.muted;
  for (const f of [0, 0.5, 1]) {
    const y = m.t + f * ph;
    ctx.beginPath(); ctx.moveTo(m.l, y); ctx.lineTo(w - m.r, y); ctx.stroke();
    ctx.textAlign = 'right';
    ctx.fillText(fmtAuto(f * ymax), m.l - 5, y + 3.5);
  }
  ctx.textAlign = 'left'; ctx.fillText('sag · mm', m.l - 30, m.t - 6);
  for (const f of [0, 0.5, 1]) {
    ctx.textAlign = f === 0 ? 'left' : f === 1 ? 'right' : 'center';
    ctx.fillText(fmt(f * d.L, 1) + ' m', m.l + f * pw, h - 6);
  }

  if (!truth) { return; }

  /* space line at 0 (draw first, under others) */
  ctx.strokeStyle = TH.space; ctx.lineWidth = 2;
  ctx.beginPath(); ctx.moveTo(m.l, Y(0)); ctx.lineTo(w - m.r, Y(0)); ctx.stroke();

  /* ideal (dashed muted) */
  ctx.strokeStyle = TH.muted; ctx.lineWidth = 1.5; ctx.setLineDash([4, 3]);
  ctx.beginPath(); ideal.forEach((v, i) => i ? ctx.lineTo(X(i), Y(v)) : ctx.moveTo(X(i), Y(v))); ctx.stroke();
  ctx.setLineDash([]);

  /* truth (accent) */
  ctx.strokeStyle = TH.accent; ctx.lineWidth = 2;
  ctx.beginPath(); truth.forEach((v, i) => i ? ctx.lineTo(X(i), Y(v)) : ctx.moveTo(X(i), Y(v))); ctx.stroke();

  /* direct labels (identity never color-alone) */
  ctx.font = mono;
  ctx.textAlign = 'right';
  ctx.fillStyle = TH.space; ctx.fillText('space', w - m.r, Y(0) - 4);
  ctx.fillStyle = TH.accent; ctx.fillText('FEM', w - m.r, Y(truth[CFG.nn - 2]) + 12);
  ctx.textAlign = 'center';
  ctx.fillStyle = TH.muted;
  ctx.fillText('ideal', X((CFG.nn - 1) / 2), Y(ideal[(CFG.nn - 1) / 2]) + 12);

  /* hover readout */
  if (hover.sec) {
    const i = Math.round((hover.sec - m.l) / pw * (CFG.nn - 1));
    if (i >= 0 && i < CFG.nn) {
      const x = X(i);
      ctx.strokeStyle = TH.axis; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(x, m.t); ctx.lineTo(x, m.t + ph); ctx.stroke();
      ctx.fillStyle = TH.accent; ctx.beginPath(); ctx.arc(x, Y(truth[i]), 2.5, 0, 7); ctx.fill();
      ctx.fillStyle = TH.ink; ctx.textAlign = x > w / 2 ? 'right' : 'left';
      ctx.fillText(fmt(i / (CFG.nn - 1) * d.L, 1) + ' m · FEM ' + fmtAuto(truth[i]) +
        ' · ideal ' + fmtAuto(ideal[i]) + ' mm', x + (x > w / 2 ? -6 : 6), m.t + 10);
    }
  }
}

/* ------------------------- slope vs force chart -------------------------- */
function drawChart() {
  const { ctx, w, h } = fitCanvas(els.cvChart);
  ctx.clearRect(0, 0, w, h);
  const d = derived();
  const mono = '10px ' + getComputedStyle(document.body).getPropertyValue('--mono');
  const m = { l: 34, r: 18, t: 16, b: 30 };
  const pw = w - m.l - m.r, ph = h - m.t - m.b;
  if (pw < 40 || ph < 30) return;

  const Fmin = 20, Fmax = 1000;
  const kIdeal = d.sigma * d.g * d.L * d.L / (2 * Math.sqrt(Math.PI)) * 1e3; // mrad·N
  const kTruth = shapes ? shapes.coeffs.ST * d.sigma * d.g * d.L * d.L * 1e3 : NaN;
  const ymax = niceCeil(Math.min(40, Math.max(2.5 * state.tgt, kIdeal / Math.max(d.F, Fmin) * 1.3,
    shapes ? kTruth / Math.max(d.F, Fmin) * 1.15 : 0, 4)));
  const X = F => m.l + (F / Fmax) * pw;
  const Y = v => m.t + ph - (Math.min(v, ymax) / ymax) * ph;

  /* grid + ticks */
  ctx.font = mono; ctx.fillStyle = TH.muted; ctx.strokeStyle = TH.grid; ctx.lineWidth = 1;
  for (const f of [0, 0.25, 0.5, 0.75, 1]) {
    const y = m.t + f * ph;
    ctx.beginPath(); ctx.moveTo(m.l, y); ctx.lineTo(w - m.r, y); ctx.stroke();
    ctx.textAlign = 'right'; ctx.fillText(fmtAuto((1 - f) * ymax), m.l - 5, y + 3.5);
  }
  for (const F of [0, 250, 500, 750, 1000]) {
    ctx.textAlign = F === 0 ? 'left' : F === 1000 ? 'right' : 'center';
    ctx.fillText(fmt(F, 0), X(F), h - 18);
  }
  ctx.textAlign = 'left'; ctx.fillText('slope · mrad', m.l - 30, m.t - 6);
  ctx.textAlign = 'center'; ctx.fillText('corner force · N', m.l + pw / 2, h - 5);

  /* target hairline */
  ctx.strokeStyle = TH.ink2; ctx.lineWidth = 1; ctx.setLineDash([2, 3]);
  ctx.beginPath(); ctx.moveTo(m.l, Y(state.tgt)); ctx.lineTo(w - m.r, Y(state.tgt)); ctx.stroke();
  ctx.setLineDash([]);
  ctx.fillStyle = TH.ink2; ctx.textAlign = 'left';
  ctx.fillText('target ' + fmt(state.tgt, 2), m.l + 4, Y(state.tgt) - 4);

  const gzero = d.g < 0.005;

  /* space series: flat at ~0 */
  ctx.strokeStyle = TH.space; ctx.lineWidth = 2;
  ctx.beginPath(); ctx.moveTo(m.l, Y(0)); ctx.lineTo(w - m.r, Y(0)); ctx.stroke();
  ctx.fillStyle = TH.space; ctx.textAlign = 'right'; ctx.fillText('space', w - m.r - 2, Y(0) - 4);

  if (!gzero) {
    /* ideal (dashed muted) & truth (accent) hyperbolas — clip above ymax
     * (enter the plot where the curve crosses the top, don't smear along it) */
    const hyperbola = k => {
      ctx.beginPath();
      let started = false;
      const Fenter = Math.max(Fmin, k / ymax);           // slope(F) <= ymax from here on
      for (let F = Fenter; F <= Fmax; F += 5) {
        const yv = Y(k / F);
        started ? ctx.lineTo(X(F), yv) : ctx.moveTo(X(F), yv);
        started = true;
      }
      ctx.stroke();
    };
    ctx.strokeStyle = TH.muted; ctx.lineWidth = 1.5; ctx.setLineDash([4, 3]);
    hyperbola(kIdeal);
    ctx.setLineDash([]);
    if (shapes) { ctx.strokeStyle = TH.accent; ctx.lineWidth = 2; hyperbola(kTruth); }

    /* direct labels sit above their own curve mid-plot, clear of the corner;
     * nudge the ideal label along its curve if it lands on the target line */
    let Flab = 0.42 * Fmax;
    if (Math.abs(Y(Math.min(kIdeal / Flab, ymax)) - Y(state.tgt)) < 14) Flab = 0.62 * Fmax;
    ctx.textAlign = 'center';
    ctx.fillStyle = TH.muted; ctx.fillText('ideal', X(Flab), Y(Math.min(kIdeal / Flab, ymax)) - 6);
    if (shapes) { ctx.fillStyle = TH.accent; ctx.fillText('FEM', X(0.42 * Fmax), Y(kTruth / (0.42 * Fmax)) + 13); }

    /* required-force marker on the axis (the number lives in the shared line) */
    const Freq = analytic.Freq(d.sigma, d.g, d.L, state.tgt / 1e3);
    if (Freq <= Fmax) {
      ctx.strokeStyle = TH.ink2; ctx.lineWidth = 1.5;
      ctx.beginPath(); ctx.moveTo(X(Freq), m.t + ph); ctx.lineTo(X(Freq), m.t + ph - 7); ctx.stroke();
      ctx.lineWidth = 1;
    }

    /* current operating point on both ground curves */
    if (d.F >= Fmin) {
      if (shapes) {
        ctx.fillStyle = TH.accent;
        ctx.beginPath(); ctx.arc(X(d.F), Y(kTruth / d.F), 3.5, 0, 7); ctx.fill();
      }
      ctx.strokeStyle = TH.muted; ctx.fillStyle = TH.surface; ctx.lineWidth = 1.5;
      ctx.beginPath(); ctx.arc(X(d.F), Y(kIdeal / d.F), 3, 0, 7); ctx.fill(); ctx.stroke();
    }
  }

  /* hover crosshair + tooltip */
  if (hover.chart && !gzero) {
    const F = Math.round((hover.chart - m.l) / pw * Fmax / 5) * 5;
    if (F >= Fmin && F <= Fmax) {
      ctx.strokeStyle = TH.axis; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(X(F), m.t); ctx.lineTo(X(F), m.t + ph); ctx.stroke();
      ctx.fillStyle = TH.ink; ctx.font = mono;
      ctx.textAlign = X(F) > w / 2 ? 'right' : 'left';
      const parts = ['F ' + fmt(F, 0) + ' N'];
      if (shapes) parts.push('FEM ' + fmtAuto(kTruth / F));
      parts.push('ideal ' + fmtAuto(kIdeal / F) + ' mrad');
      ctx.fillText(parts.join(' · '), X(F) + (X(F) > w / 2 ? -6 : 6), m.t + 10);
    }
  }
}

function niceCeil(v) {
  if (!isFinite(v) || v <= 0) return 1;
  const e = Math.pow(10, Math.floor(Math.log10(v)));
  for (const m of [1, 1.5, 2, 2.5, 3, 4, 5, 6, 8, 10]) if (m * e >= v) return m * e;
  return 10 * e;
}

/* ------------------------- slew billow chart -----------------------------
 * Repointing in orbit: an overhead tracking pass at line-of-sight rate omega
 * peaks at angular acceleration alpha ≈ 0.65·omega². The tangential
 * acceleration alpha·x loads the film antisymmetrically about the slew axis;
 * its RMS slope adds to the (orthogonal) solar-pressure billow in quadrature. */
function slewSlopeMrad(d, omegaDegS) {
  if (!shapes) return NaN;
  const alpha = SLEW_ACC * Math.pow(omegaDegS * Math.PI / 180, 2);  // rad/s²
  return shapes.coeffs.SS * d.sigma * alpha * Math.pow(d.L, 3) / Math.max(d.F, 1) * 1e3;
}
function spaceTotalMrad(d, omegaDegS) {
  const srp = shapes ? shapes.coeffs.ST * d.sScaleSun * 1e3 : NaN;
  return Math.hypot(srp, slewSlopeMrad(d, omegaDegS));
}
function drawSlew() {
  const { ctx, w, h } = fitCanvas(els.cvSlew);
  ctx.clearRect(0, 0, w, h);
  const d = derived();
  const mono = '10px ' + getComputedStyle(document.body).getPropertyValue('--mono');
  const m = { l: 44, r: 18, t: 16, b: 30 };
  const pw = w - m.l - m.r, ph = h - m.t - m.b;
  if (pw < 40 || ph < 30 || !shapes || d.slack) return;

  const OMAX = 5;
  const srp = shapes.coeffs.ST * d.sScaleSun * 1e3;
  const ymax = niceCeil(Math.max(spaceTotalMrad(d, OMAX) * 1.15, srp * 2, 1e-6));
  const X = o => m.l + (o / OMAX) * pw;
  const Y = v => m.t + ph - (Math.min(v, ymax) / ymax) * ph;

  /* grid + ticks */
  ctx.font = mono; ctx.fillStyle = TH.muted; ctx.strokeStyle = TH.grid; ctx.lineWidth = 1;
  for (const f of [0, 0.5, 1]) {
    const y = m.t + f * ph;
    ctx.beginPath(); ctx.moveTo(m.l, y); ctx.lineTo(w - m.r, y); ctx.stroke();
    ctx.textAlign = 'right'; ctx.fillText(fmtSmart((1 - f) * ymax), m.l - 5, y + 3.5);
  }
  for (const o of [0, 1, 2, 3, 4, 5]) {
    ctx.textAlign = o === 0 ? 'left' : o === OMAX ? 'right' : 'center';
    ctx.fillText(fmt(o, 0), X(o), h - 18);
  }
  ctx.textAlign = 'left'; ctx.fillText('slope · mrad', m.l - 34, m.t - 6);
  ctx.textAlign = 'center'; ctx.fillText('slew rate · °/s', m.l + pw / 2, h - 5);

  /* sunlight-only floor */
  ctx.strokeStyle = TH.muted; ctx.setLineDash([4, 3]); ctx.lineWidth = 1.5;
  ctx.beginPath(); ctx.moveTo(m.l, Y(srp)); ctx.lineTo(w - m.r, Y(srp)); ctx.stroke();
  ctx.setLineDash([]);
  ctx.fillStyle = TH.muted; ctx.textAlign = 'left';
  ctx.fillText('sunlight only', m.l + 4, Y(srp) - 4);

  /* total curve */
  ctx.strokeStyle = TH.space; ctx.lineWidth = 2;
  ctx.beginPath();
  for (let o = 0; o <= OMAX + 1e-9; o += OMAX / 120) {
    const y = Y(spaceTotalMrad(d, o));
    o === 0 ? ctx.moveTo(X(o), y) : ctx.lineTo(X(o), y);
  }
  ctx.stroke();
  ctx.fillStyle = TH.space; ctx.textAlign = 'right';
  ctx.fillText('while slewing', w - m.r - 2, Y(spaceTotalMrad(d, OMAX)) - 6);

  /* target line only when it is on scale */
  if (state.tgt <= ymax) {
    ctx.strokeStyle = TH.ink2; ctx.setLineDash([2, 3]);
    ctx.beginPath(); ctx.moveTo(m.l, Y(state.tgt)); ctx.lineTo(w - m.r, Y(state.tgt)); ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = TH.ink2; ctx.textAlign = 'left';
    ctx.fillText('target ' + fmt(state.tgt, 2), m.l + 4, Y(state.tgt) - 4);
  }

  /* current operating point */
  ctx.fillStyle = TH.space;
  ctx.beginPath(); ctx.arc(X(state.slew), Y(spaceTotalMrad(d, state.slew)), 3.5, 0, 7); ctx.fill();

  /* hover */
  if (hover.slew != null) {
    const o = Math.max(0, Math.min(OMAX, (hover.slew - m.l) / pw * OMAX));
    ctx.strokeStyle = TH.axis; ctx.lineWidth = 1;
    ctx.beginPath(); ctx.moveTo(X(o), m.t); ctx.lineTo(X(o), m.t + ph); ctx.stroke();
    ctx.fillStyle = TH.ink; ctx.textAlign = X(o) > w / 2 ? 'right' : 'left';
    ctx.fillText(fmt(o, 2) + ' °/s · ' + fmtSmart(spaceTotalMrad(d, o)) + ' mrad',
      X(o) + (X(o) > w / 2 ? -6 : 6), m.t + 10);
  }
}

/* ------------------------- slew maneuver preview -------------------------
 * A repoint cycle (time-compressed): the mirror rocks between ±Θ; the film
 * lags the acceleration and S-bows, flattening through the coast. The billow
 * is drawn hugely exaggerated but its amplitude scales honestly with rate²,
 * so dragging the slider up makes the lag visible. */
const slewAnim = { raf: 0, t: 0, last: 0 };
const SLEW_T = 6;                                        // s per full cycle (compressed)
function slewAnimActive() {
  return !!shapes && !document.hidden && !reduceMotion && state.slew > 0.02;
}
function scheduleSlewAnim() {
  if (!slewAnim.raf && slewAnimActive()) slewAnim.raf = requestAnimationFrame(slewAnimTick);
}
function slewAnimTick(ts) {
  slewAnim.raf = 0;
  const dt = slewAnim.last ? Math.min(0.05, (ts - slewAnim.last) / 1000) : 0.016;
  slewAnim.last = ts;
  slewAnim.t += dt;
  drawSlewAnim();
  scheduleSlewAnim();
}
function drawSlewAnim() {
  const cv = els.cvSlewAnim;
  if (!cv) return;
  const { ctx, w, h } = fitCanvas(cv);
  ctx.clearRect(0, 0, w, h);
  const d = derived();
  const mono = '10px ' + getComputedStyle(document.body).getPropertyValue('--mono');
  if (!shapes || d.slack) {
    ctx.fillStyle = TH.muted; ctx.font = mono; ctx.textAlign = 'center';
    ctx.fillText(d.slack ? 'membrane slack' : '…', w / 2, h / 2);
    return;
  }
  const L = d.L;
  const pad = 8;
  const sc = Math.min((w - 2 * pad) / (2 * L * KX), (h - 2 * pad - 26) / (2 * L * KY + 0.5 * L));
  const cx = w / 2, cy = h / 2 - 2;
  const phase = (slewAnim.t % SLEW_T) / SLEW_T * 2 * Math.PI;
  const still = !slewAnimActive();
  const theta = still ? 0 : 0.20 * Math.sin(phase);      // pointing angle, rad (visual)
  const aNorm = still ? 1 : -Math.sin(phase);            // billow follows the acceleration
  /* drawn billow amplitude: honest ∝ rate², softly saturating in-frame */
  const cap = 0.22 * L;
  const Araw = 0.16 * L * Math.pow(state.slew / 5, 2) * 8;
  const Aeff = cap * Math.tanh(Araw / cap) * aNorm;
  const nn = shapes.mesh.nn, hh = shapes.mesh.h, STEP = 2;
  const norm = 1 / Math.max(shapes.maxAbsS, 1e-12);
  const proj = (X, Y, z) => {
    /* rotate about the in-plane centerline axis (parallel to Y edges) */
    const Xr = X * Math.cos(theta) - z * Math.sin(theta);
    const zr = X * Math.sin(theta) + z * Math.cos(theta);
    return [cx + (Xr - Y) * KX * sc, cy + (Xr + Y) * KY * sc + zr * KZ * sc];
  };
  const pts = new Map();
  const P = (i, j) => {
    const k = j * nn + i;
    let p = pts.get(k);
    if (!p) {
      const X = (i * hh - 0.5) * L, Y = (j * hh - 0.5) * L;
      p = proj(X, Y, shapes.whatS[k] * norm * Aeff);
      pts.set(k, p);
    }
    return p;
  };
  /* slew axis, dashed */
  ctx.strokeStyle = TH.axis; ctx.lineWidth = 1; ctx.setLineDash([3, 3]);
  const a0 = proj(0, -0.62 * L, 0), a1 = proj(0, 0.62 * L, 0);
  ctx.beginPath(); ctx.moveTo(a0[0], a0[1]); ctx.lineTo(a1[0], a1[1]); ctx.stroke();
  ctx.setLineDash([]);
  /* surface, back to front, shaded by local billow slope */
  const scale01 = Math.abs(Aeff) / cap;
  if (!shapes.maxSlopeS) shapes.maxSlopeS = shapes.slopesS.reduce((a, v) => Math.max(a, v), 1e-12);
  for (let s = 0; s <= 2 * (nn - 1 - STEP); s += STEP) {
    for (let i = Math.max(0, s - (nn - 1 - STEP)); i <= Math.min(nn - 1 - STEP, s); i += STEP) {
      const j = s - i;
      if (j > nn - 1 - STEP || j % STEP) continue;
      const e = j * (nn - 1) + i;
      ctx.fillStyle = rampColor(shapes.slopesS[e] / shapes.maxSlopeS * scale01);
      const p0 = P(i, j), p1 = P(i + STEP, j), p2 = P(i + STEP, j + STEP), p3 = P(i, j + STEP);
      ctx.beginPath();
      ctx.moveTo(p0[0], p0[1]); ctx.lineTo(p1[0], p1[1]);
      ctx.lineTo(p2[0], p2[1]); ctx.lineTo(p3[0], p3[1]);
      ctx.closePath(); ctx.fill();
      ctx.strokeStyle = ctx.fillStyle; ctx.lineWidth = 0.6; ctx.stroke();
    }
  }
  /* outline + corner marks */
  ctx.strokeStyle = inkA(.45); ctx.lineWidth = 1;
  ctx.beginPath();
  const last = nn - 1;
  const ring = [];
  for (let i = 0; i <= last; i += STEP) ring.push([i, 0]);
  for (let j = STEP; j <= last; j += STEP) ring.push([last, j]);
  for (let i = last - STEP; i >= 0; i -= STEP) ring.push([i, last]);
  for (let j = last - STEP; j >= STEP; j -= STEP) ring.push([0, j]);
  ring.forEach(([i, j], k) => { const p = P(i, j); k ? ctx.lineTo(p[0], p[1]) : ctx.moveTo(p[0], p[1]); });
  ctx.closePath(); ctx.stroke();
  ctx.fillStyle = TH.ink;
  for (const [i, j] of [[0, 0], [last, 0], [0, last], [last, last]]) {
    const p = P(i, j); ctx.fillRect(p[0] - 1.5, p[1] - 1.5, 3, 3);
  }
  /* annotation */
  ctx.fillStyle = TH.muted; ctx.font = mono; ctx.textAlign = 'right';
  ctx.fillText((still ? 'paused · ' : 'time compressed · ') +
    'billow exaggerated, grows ∝ rate²', w - 4, h - 4);
}

function drawAll() {
  const d = derived();
  drawMembrane(els.g.cv, d, true, anim.gAmp.v);
  drawMembrane(els.s.cv, d, false, 1);
  drawSection();
  drawChart();
  drawSlew();
  if (!slewAnimActive()) drawSlewAnim();                 // static frame when not animating
}

/* ------------------------------ controls -------------------------------- */
function bindPair(rngId, numId, key) {
  const rng = $(rngId), num = $(numId);
  const set = (v, live) => {
    state[key] = v;
    if (rng) rng.value = v;
    num.value = v;
    onControl(key, live);
  };
  if (rng) {
    rng.addEventListener('input', () => set(parseFloat(rng.value), true));  // no snap mid-drag
    rng.addEventListener('change', () => set(parseFloat(rng.value)));       // snap on release
  }
  num.addEventListener('change', () => {
    let v = parseFloat(num.value);
    if (!isFinite(v)) v = DEFAULTS[key];
    v = Math.min(parseFloat(num.max), Math.max(parseFloat(num.min), v));
    set(v);
  });
  return set;
}
const setters = {
  F: bindPair('rngF', 'numF', 'F'),
  phiDeg: bindPair('rngPhi', 'numPhi', 'phiDeg'),
  g: bindPair('rngG', 'numG', 'g'),
  L: bindPair('rngL', 'numL', 'L'),
  massG: bindPair('rngM', 'numM', 'massG'),
  tgt: bindPair(null, 'numTgt', 'tgt'),
  E_GPa: bindPair(null, 'numE', 'E_GPa'),
  rho: bindPair(null, 'numRho', 'rho'),
  slew: bindPair('rngSlew', 'numSlew', 'slew'),
};

const G_PRESETS = [0, 1.62, 3.72, 9.81];
const presetBtns = [...document.querySelectorAll('.presets button')];
function onControl(key, live) {
  if (key === 'g') {
    /* snap the star control to nearby presets — only on release, not mid-drag */
    if (!live)
      for (const p of G_PRESETS)
        if (Math.abs(state.g - p) < 0.22 && state.g !== p) {
          state.g = p; $('rngG').value = p; $('numG').value = p; break;
        }
    presetBtns.forEach(b => b.setAttribute('aria-pressed', String(+b.dataset.g === state.g)));
  }
  if (key === 'phiDeg') requestShapes();               // may need a fresh solve
  anim.gAmp.t = 1;
  update();
  scheduleTick();
  scheduleSlewAnim();
}
presetBtns.forEach(b => b.addEventListener('click', () => {
  /* discrete jump: let the membrane re-settle from a shallower state */
  if (!reduceMotion) anim.gAmp.v = 0.25;
  setters.g(+b.dataset.g);
}));
$('btnReset').addEventListener('click', () => {
  setters.F(DEFAULTS.F); setters.phiDeg(DEFAULTS.phiDeg); setters.g(DEFAULTS.g);
  setters.L(DEFAULTS.L); setters.massG(DEFAULTS.massG); setters.tgt(DEFAULTS.tgt);
  setters.E_GPa(DEFAULTS.E_GPa); setters.rho(DEFAULTS.rho); setters.slew(DEFAULTS.slew);
});

/* hover wiring */
function wireHover(cv, key, redraw) {
  cv.addEventListener('mousemove', e => { hover[key] = e.offsetX; redraw(); });
  cv.addEventListener('mouseleave', () => { hover[key] = null; redraw(); });
}
wireHover(els.cvSec, 'sec', drawSection);
wireHover(els.cvChart, 'chart', drawChart);
wireHover(els.cvSlew, 'slew', drawSlew);

/* resize */
new ResizeObserver(() => drawAll()).observe(document.querySelector('main'));

/* reflect (possibly query-overridden) state into the controls */
for (const [ids, key] of [[['rngF', 'numF'], 'F'], [['rngPhi', 'numPhi'], 'phiDeg'],
  [['rngG', 'numG'], 'g'], [['rngL', 'numL'], 'L'], [['rngM', 'numM'], 'massG'],
  [['numTgt'], 'tgt'], [['numE'], 'E_GPa'], [['numRho'], 'rho'],
  [['rngSlew', 'numSlew'], 'slew']])
  for (const id of ids) $(id).value = state[key];
presetBtns.forEach(b => b.setAttribute('aria-pressed', String(+b.dataset.g === state.g)));

/* theme: the pre-paint <head> script already chose (?theme= or OS); adopt it */
setTheme(document.documentElement.dataset.theme === 'dark' ? 'dark' : 'light');
matchMedia('(prefers-color-scheme: dark)')
  .addEventListener('change', e => setTheme(e.matches ? 'dark' : 'light'));
$('themeBtn').addEventListener('click', () =>
  setTheme(document.documentElement.dataset.theme === 'dark' ? 'light' : 'dark'));

document.addEventListener('visibilitychange', () => { slewAnim.last = 0; scheduleSlewAnim(); });

update();
})();
