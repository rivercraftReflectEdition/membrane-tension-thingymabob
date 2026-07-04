'use strict';
/* Acceptance tests for the physics core SHIPPED inside index.html.
 * The core is extracted from between the @physics-core sentinels and run
 * here, so what is tested is byte-for-byte what the page executes.
 * Run: node test/physics.test.js */

const fs = require('fs');
const path = require('path');

const html = fs.readFileSync(path.join(__dirname, '..', 'index.html'), 'utf8');
const m = html.match(/\/\* @physics-core-begin \*\/([\s\S]*?)\/\* @physics-core-end \*\//);
if (!m) { console.error('FAIL: physics-core sentinels not found in index.html'); process.exit(1); }

/* the built page must be fresh: its core must equal src/physics-core.js */
const srcCore = fs.readFileSync(path.join(__dirname, '..', 'src', 'physics-core.js'), 'utf8');
if (m[1].trim() !== srcCore.trim()) {
  console.error('FAIL: index.html is stale — run `node build.js`');
  process.exit(1);
}
console.log('PASS built index.html matches src/physics-core.js');

const module_ = { exports: {} };
(function (module) { eval(m[1]); })(module_);
const core = module_.exports;
const A = core.analytic;

let failures = 0;
function check(name, got, want, relTol) {
  const ok = Math.abs(got - want) <= Math.abs(want) * relTol;
  console.log('%s %s: got %s, want %s ±%d%%', ok ? 'PASS' : 'FAIL', name,
    Number(got.toPrecision(4)), Number(want.toPrecision(4)), relTol * 100);
  if (!ok) failures++;
}

/* -------- spec defaults: L=15 m, mass=0.35 kg, F=250 N, CP1 -------- */
const L = 15, mass = 0.35, F = 250, g = 9.81, rho = 1430;
const sigma = A.sigma(mass, L);
const t = A.thick(mass, rho, L);
const N = A.tension(F, L);
const p = A.pressure(sigma, g);

console.log('--- analytic (uniform-tension reference), spec sanity numbers ---');
check('areal density (g/m^2)', sigma * 1e3, 1.56, 0.02);
check('film thickness (um)', t * 1e6, 1.1, 0.05);
check('tension N (N/m)', N, 11.8, 0.02);
check('RMS slope (mrad)', A.rmsSlope(p, L, N) * 1e3, 3.9, 0.05);
check('central sag (mm)', A.sagC(p, L, N) * 1e3, 23, 0.05);
check('F for 1 mrad (N)', A.Freq(sigma, g, L, 1e-3), 970, 0.02);
check('F for 0.5 mrad ~ 2 kN', A.Freq(sigma, g, L, 0.5e-3), 1940, 0.02);

console.log('--- discrete FEM shapes (app startup solve, nn=41) ---');
const sh = core.buildShapes(41, 0.34, { maxIter: 30, floorFrac: 0.02 });
const sScale = sigma * g * L * L / F, wScale = sScale * L;

/* uniform-tension solver must reproduce the EXACT square drumhead */
check('ref center coeff (exact 0.10419)', sh.ref.centerCoeff, Math.SQRT2 * 0.073671, 0.01);
/* and sit within the equivalent-circle idealization of the closed forms */
check('ref center sag vs analytic formula (mm)', sh.ref.centerCoeff * wScale * 1e3,
  A.sagC(p, L, N) * 1e3, 0.10);
check('ref RMS slope vs analytic formula (mrad)', sh.ref.rmsSlopeCoeff * sScale * 1e3,
  A.rmsSlope(p, L, N) * 1e3, 0.10);

/* corner-loaded truth: converged coefficients (regression band from the
 * mesh-convergence study: S 0.144-0.150 over nn=31..51, PV ~0.101) */
check('truth RMS slope coeff', sh.truth.rmsSlopeCoeff, 0.146, 0.08);
check('truth P-V coeff', sh.truth.pvCoeff, 0.101, 0.06);
check('truth P-V sag at defaults (mm)', sh.truth.pvCoeff * wScale * 1e3, 21, 0.10);

/* static equilibrium: tension flux across the horizontal midline must carry
 * the two corner-force components, sqrt(2) in unit terms */
{
  const nn = 41, mesh = sh.mesh, ip = sh.inPlane;
  const jrow = (nn - 1) / 2 | 0;
  let flux = 0;
  for (let i = 0; i < nn - 1; i++) flux += ip.nhat[(jrow * (nn - 1) + i) * 3 + 1] * mesh.h;
  check('midline tension flux (= sqrt 2)', flux, Math.SQRT2, 0.01);
}

/* space: zero pressure -> exactly flat (pure scaling: slope = coeff * sigma*0*L^2/F) */
check('space RMS slope (mrad)', sh.truth.rmsSlopeCoeff * (sigma * 0 * L * L / F) * 1e3 + 1, 1, 1e-12);

/* payload fields the renderer depends on */
for (const [k, v] of Object.entries({
  slopesT: sh.truth.slopes.length, slopesR: sh.ref.slopes.length,
  thetaW: sh.inPlane.thetaW.length, state: sh.inPlane.state.length,
})) if (v !== 1600) { console.log('FAIL renderer field ' + k + ' length ' + v); failures++; }

console.log('--- catenary edge cords ---');
/* interior tension must match the closed form q = sqrt(2)·F·sin(phi)/(L·(cos+sin));
 * at phi = 45° the cords reproduce the conventional N = F/(sqrt(2)·L) exactly */
for (const deg of [10, 45]) {
  const phi = deg * Math.PI / 180;
  const shc = core.buildShapes(41, 0.34, { maxIter: 30, floorFrac: 0.02, phi });
  const mid = 20 * 40 + 20;
  const qExp = Math.SQRT2 * Math.sin(phi) / (Math.cos(phi) + Math.sin(phi));
  check('phi=' + deg + ' interior tension', shc.inPlane.nhat[mid * 3], qExp, 0.05);
  if (deg === 45)
    check('phi=45 RMS coeff -> analytic uniform', shc.truth.rmsSlopeCoeff, 0.2821, 0.15);
  if (deg === 10) {
    check('phi=10 RMS coeff (regression)', shc.truth.rmsSlopeCoeff, 0.887, 0.08);
    /* slew shape: antisymmetric billow about the centerline */
    const c = shc.slew.what[(41 * 41 - 1) / 2];
    if (Math.abs(c) > 1e-9) { console.log('FAIL slew centerline not ~0: ' + c); failures++; }
    else console.log('PASS slew shape antisymmetric (centerline ~0)');
    check('phi=10 slew RMS coeff (regression)', shc.slew.rmsSlopeCoeff, 0.148, 0.08);
  }
}

console.log(failures ? '\n' + failures + ' FAILURE(S)' : '\nall checks passed');
process.exit(failures ? 1 : 0);
