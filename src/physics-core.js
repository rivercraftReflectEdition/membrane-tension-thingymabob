/* ============================================================================
 * Membrane mirror physics core — dimensionless FEM on the unit square.
 *
 * Physics: a thin film carries load only in tension (no bending), so its
 * out-of-plane shape w obeys the drumhead equation
 *      div( N · grad w ) = -p          (N = 2x2 tension-resultant field, N/m)
 * with p = sigma·g (transverse gravity pressure, Pa) and w pinned only at
 * the 4 boom attachment corners. The in-plane tension field N comes from a
 * plane-stress solve with outward diagonal point forces F at the corners,
 * with a wrinkling model (a membrane cannot carry compression).
 *
 * Everything is solved ONCE in dimensionless form, then scaled:
 *   x = L·xi,   N = (F/L)·nhat,   w = (sigma·g·L^3 / F) · what
 * so  div(nhat grad what) = -1  on the unit square. The stress-field *shape*
 * depends only on Poisson's ratio (Michell), not on E or F, and the wrinkle
 * state is scale-invariant under pure corner loading — hence one solve at
 * startup covers every slider position. E enters only the in-plane strain
 * (edge scalloping displacement, u = uhat · F(1-nu^2)/(E t)) and rho only
 * the film thickness t = mass/(rho L^2).
 *
 * Small-slope (linear) regime throughout: results are in mrad, so the
 * geometric stiffening nonlinearity is negligible.
 *
 * v1 hooks (not modeled): thermal wrinkling, boom buckling, solar radiation
 * pressure (~9 uPa — would enter as an extra `p` term in space).
 * ==========================================================================*/

'use strict';

/* ---------- banded symmetric positive-definite solver (Cholesky LL^T) ----
 * Lower band storage: B[i*(bw+1)+d] = A[i][i-d], d = 0..bw.
 * Deterministic direct solve — an explicit time-stepper would need
 * femtosecond steps for a 1 um film's in-plane stiffness; we only want
 * static equilibrium, so we solve for it directly. */
function bandFactor(B, n, bw) {
  const w = bw + 1;
  for (let i = 0; i < n; i++) {
    const kmin = Math.max(0, i - bw);
    for (let j = kmin; j <= i; j++) {
      let s = B[i * w + (i - j)];
      const m = Math.max(kmin, j - bw);
      for (let k = m; k < j; k++) s -= B[i * w + (i - k)] * B[j * w + (j - k)];
      if (j === i) {
        if (s <= 0) throw new Error('matrix not positive definite at row ' + i);
        B[i * w] = Math.sqrt(s);
      } else {
        B[i * w + (i - j)] = s / B[j * w];
      }
    }
  }
}
function bandSolve(B, n, bw, rhs) {
  const w = bw + 1, x = Float64Array.from(rhs);
  for (let i = 0; i < n; i++) {                       // L y = b
    let s = x[i];
    const kmin = Math.max(0, i - bw);
    for (let k = kmin; k < i; k++) s -= B[i * w + (i - k)] * x[k];
    x[i] = s / B[i * w];
  }
  for (let i = n - 1; i >= 0; i--) {                  // L^T x = y
    let s = x[i];
    const kmax = Math.min(n - 1, i + bw);
    for (let k = i + 1; k <= kmax; k++) s -= B[k * w + (k - i)] * x[k];
    x[i] = s / B[i * w];
  }
  return x;
}

/* ---------- mesh: nn x nn nodes on [0,1]^2, bilinear quads ---------- */
function makeMesh(nn) {
  const h = 1 / (nn - 1);
  const nodes = nn * nn;
  const ne = (nn - 1) * (nn - 1);
  const elems = new Int32Array(ne * 4);               // ccw: n0 n1 n2 n3
  let e = 0;
  for (let j = 0; j < nn - 1; j++)
    for (let i = 0; i < nn - 1; i++) {
      const n0 = j * nn + i;
      elems[e * 4] = n0; elems[e * 4 + 1] = n0 + 1;
      elems[e * 4 + 2] = n0 + nn + 1; elems[e * 4 + 3] = n0 + nn;
      e++;
    }
  return { nn, h, nodes, ne, elems };
}

/* 2x2 Gauss points on the reference square [-1,1]^2 */
const GP = (() => {
  const a = 1 / Math.sqrt(3), pts = [];
  for (const xi of [-a, a]) for (const eta of [-a, a]) pts.push([xi, eta]);
  return pts;
})();

/* bilinear shape-function gradients in physical coords for square element of size h */
function shapeGrad(xi, eta, h) {
  const s = 2 / h;                                    // d(xi)/dx
  // N0=(1-xi)(1-eta)/4, N1=(1+xi)(1-eta)/4, N2=(1+xi)(1+eta)/4, N3=(1-xi)(1+eta)/4
  return {
    N:  [(1-xi)*(1-eta)/4, (1+xi)*(1-eta)/4, (1+xi)*(1+eta)/4, (1-xi)*(1+eta)/4],
    dx: [-(1-eta)/4*s, (1-eta)/4*s, (1+eta)/4*s, -(1+eta)/4*s],
    dy: [-(1-xi)/4*s, -(1+xi)/4*s, (1+xi)/4*s, (1-xi)/4*s],
  };
}

/* ============================================================================
 * IN-PLANE: plane stress with corner loads + wrinkling iteration.
 * Dimensionless: unit corner force, membrane modulus Et/(1-nu^2) = 1.
 * Returns per-element tension resultants nhat = [nxx, nyy, nxy] (units F/L)
 * and nodal displacements uhat (units F(1-nu^2)/(Et)).
 * Wrinkling: iterative membrane properties (Miller–Hedgepeth style) — where
 * the minor principal stress goes compressive the element "wrinkles": its
 * stiffness collapses to uniaxial along the major principal strain direction,
 * so the converged field carries (nearly) no compression. This is what makes
 * the free edges scallop and the corners carry the load.
 * ==========================================================================*/
function solveInPlane(mesh, nu, opts = {}) {
  const { nn, h, nodes, ne, elems } = mesh;
  const maxIter = opts.maxIter ?? 80;
  const kResid  = opts.kResid  ?? 1e-3;               // residual stiffness of wrinkled/slack dirs
  const rTab    = opts.rTab    ?? 0.03;               // corner tab radius (fraction of side) —
                                                      // real corners are reinforced tabs/cords,
                                                      // not points; also regularizes the FEM
  const phi     = opts.phi     ?? 0;                  // catenary cord arc half-angle (rad); 0 = no cords
  const ndof = nodes * 2;
  const bw = 2 * (nn + 1) + 1;

  const Dtaut = [[1, nu, 0], [nu, 1, 0], [0, 0, (1 - nu) / 2]];
  // per-element constitutive matrix (row-major 3x3), start taut
  const Dm = new Float64Array(ne * 9);
  for (let e = 0; e < ne; e++)
    for (let r = 0; r < 3; r++) for (let c = 0; c < 3; c++)
      Dm[e * 9 + r * 3 + c] = Dtaut[r][c];
  const state = new Int8Array(ne);                    // 0 taut, 1 wrinkled, 2 slack

  /* Load system. Two attachment schemes:
   *
   * phi = 0 — direct corner attach: unit force per corner, outward along the
   * diagonal, spread over the tab nodes (all nodes within rTab) with a linear
   * taper.
   *
   * phi > 0 — catenary edge cords: each edge is scalloped as a circular arc
   * of half-angle phi with a cord along it. Corner balance (two cord ends at
   * 45°∓phi from the boom diagonal) gives the cord tension
   *      T = F / (√2·(cos phi + sin phi))
   * and Young–Laplace turns cord tension into an outward normal line load on
   * the film edge, q = T/R with R = L/(2·sin phi):
   *      q = 2·T·sin phi / L
   * The boom force and the two cord-end pulls cancel at the tab, so the film
   * sees pure distributed edge traction — this is exactly what catenaries are
   * for. (At phi → 0 a straight cord routes the boom force corner-to-corner
   * without tensioning the film at all; the film goes slack. Real.)
   * Square-domain approximation: the traction acts on the straight edge; the
   * few-percent scallop depth is drawn but not meshed. */
  const f = new Float64Array(ndof);
  const c = Math.SQRT1_2;
  const cxy = [[0, 0], [1, 0], [0, 1], [1, 1]];
  const dir = [[-c, -c], [c, -c], [-c, c], [c, c]];
  const tabNodes = [];                                 // for out-of-plane pinning
  for (let k = 0; k < 4; k++) {
    const w = [], ids = [];
    for (let j = 0; j < nn; j++) for (let i = 0; i < nn; i++) {
      const d = Math.hypot(i * h - cxy[k][0], j * h - cxy[k][1]);
      if (d <= rTab + 1e-12) { ids.push(j * nn + i); w.push(1 - d / (rTab * 1.0001)); }
    }
    if (phi <= 0) {
      const wsum = w.reduce((a, b) => a + b, 0);
      for (let m = 0; m < ids.length; m++) {
        f[2 * ids[m]]     += dir[k][0] * w[m] / wsum;
        f[2 * ids[m] + 1] += dir[k][1] * w[m] / wsum;
      }
    }
    tabNodes.push(...ids);
  }
  const cordT = phi > 0 ? 1 / (Math.SQRT2 * (Math.cos(phi) + Math.sin(phi))) : 0;
  if (phi > 0) {
    const q = 2 * cordT * Math.sin(phi);               // per unit length, outward normal
    const qh = q * h / 2;                              // consistent nodal load per segment end
    for (let i = 0; i < nn - 1; i++) {
      for (const n of [i, i + 1])                f[2 * n + 1]                    -= qh; // bottom, (0,-1)
      for (const n of [i, i + 1])                f[2 * ((nn - 1) * nn + n) + 1]  += qh; // top, (0,+1)
      for (const j of [i, i + 1])                f[2 * (j * nn)]                 -= qh; // left, (-1,0)
      for (const j of [i, i + 1])                f[2 * (j * nn + nn - 1)]        += qh; // right, (+1,0)
    }
  }

  /* pin rigid-body modes at the center node (+ one rotation dof) —
   * loads are self-equilibrated so reactions are ~0 */
  const mid = (nn - 1) / 2;
  const fixed = new Set([
    2 * (mid * nn + mid), 2 * (mid * nn + mid) + 1,   // center: u = v = 0
    2 * (mid * nn + nn - 1) + 1,                      // mid-right edge: v = 0 (rotation)
  ]);

  const detJw = (h / 2) * (h / 2);                    // per Gauss point
  let u = new Float64Array(ndof);
  let strains = new Float64Array(ne * 3);
  let stresses = new Float64Array(ne * 3);
  let iters = 0, lastChanged = 0, maxDD = 0;

  for (let iter = 0; iter < maxIter; iter++) {
    maxDD = 0;
    /* assemble */
    const B = new Float64Array(ndof * (bw + 1));
    for (let e = 0; e < ne; e++) {
      const D = Dm.subarray(e * 9, e * 9 + 9);
      const Ke = new Float64Array(64);
      for (const [xi, eta] of GP) {
        const g = shapeGrad(xi, eta, h);
        // Bmat rows: [exx, eyy, gxy], cols: 8 dofs
        for (let a = 0; a < 4; a++) for (let b = 0; b < 4; b++) {
          const Ba = [[g.dx[a], 0], [0, g.dy[a]], [g.dy[a], g.dx[a]]];
          const Bb = [[g.dx[b], 0], [0, g.dy[b]], [g.dy[b], g.dx[b]]];
          for (let p = 0; p < 2; p++) for (let q = 0; q < 2; q++) {
            let s = 0;
            for (let r = 0; r < 3; r++) for (let cc = 0; cc < 3; cc++)
              s += Ba[r][p] * D[r * 3 + cc] * Bb[cc][q];
            Ke[(a * 2 + p) * 8 + (b * 2 + q)] += s * detJw;
          }
        }
      }
      for (let a = 0; a < 4; a++) for (let p = 0; p < 2; p++) {
        const gi = 2 * elems[e * 4 + a] + p;
        for (let b = 0; b < 4; b++) for (let q = 0; q < 2; q++) {
          const gj = 2 * elems[e * 4 + b] + q;
          if (gj <= gi) B[gi * (bw + 1) + (gi - gj)] += Ke[(a * 2 + p) * 8 + (b * 2 + q)];
        }
      }
    }
    /* apply pins (large-diagonal method keeps the band symmetric & SPD) */
    const rhs = Float64Array.from(f);
    for (const d of fixed) { B[d * (bw + 1)] += 1e8; rhs[d] = 0; }

    bandFactor(B, ndof, bw);
    u = bandSolve(B, ndof, bw, rhs);

    /* element-center strains & stresses, wrinkle-state update */
    let changed = 0;
    for (let e = 0; e < ne; e++) {
      const g = shapeGrad(0, 0, h);
      let exx = 0, eyy = 0, gxy = 0;
      for (let a = 0; a < 4; a++) {
        const ua = u[2 * elems[e * 4 + a]], va = u[2 * elems[e * 4 + a] + 1];
        exx += g.dx[a] * ua; eyy += g.dy[a] * va;
        gxy += g.dy[a] * ua + g.dx[a] * va;
      }
      strains[e * 3] = exx; strains[e * 3 + 1] = eyy; strains[e * 3 + 2] = gxy;

      /* trial stress with TAUT material decides the state (M–H criterion) */
      const sx = Dtaut[0][0] * exx + Dtaut[0][1] * eyy;
      const sy = Dtaut[1][0] * exx + Dtaut[1][1] * eyy;
      const txy = Dtaut[2][2] * gxy;
      const sm = (sx + sy) / 2, rr = Math.hypot((sx - sy) / 2, txy);
      const s2 = sm - rr;
      const em = (exx + eyy) / 2, er = Math.hypot((exx - eyy) / 2, gxy / 2);
      const e1 = em + er;

      let ns;                                          // new state
      if (s2 > 0) ns = 0;
      else if (e1 <= 0) ns = 2;
      else ns = 1;
      if (ns !== state[e]) changed++;
      state[e] = ns;

      /* target constitutive for next pass */
      let Dt;
      if (ns === 0) Dt = Dtaut;
      else if (ns === 2) Dt = [[kResid,0,0],[0,kResid,0],[0,0,kResid/2]];
      else {
        /* uniaxial along major principal STRAIN direction theta */
        const th = 0.5 * Math.atan2(gxy, exx - eyy);
        const cs = Math.cos(th), sn = Math.sin(th);
        // strain transform T (engineering shear), D_glob = T^T D* T with D* = diag(1, k, k/2·2)
        const T = [
          [cs*cs, sn*sn, cs*sn],
          [sn*sn, cs*cs, -cs*sn],
          [-2*cs*sn, 2*cs*sn, cs*cs - sn*sn],
        ];
        /* uniaxial modulus is Et, i.e. (1 - nu^2) in units of Et/(1-nu^2) */
        const Ds = [[1 - nu * nu, 0, 0], [0, kResid, 0], [0, 0, kResid]];
        Dt = [[0,0,0],[0,0,0],[0,0,0]];
        for (let r = 0; r < 3; r++) for (let cc = 0; cc < 3; cc++) {
          let s = 0;
          for (let k1 = 0; k1 < 3; k1++) for (let k2 = 0; k2 < 3; k2++)
            s += T[k1][r] * Ds[k1][k2] * T[k2][cc];
          Dt[r][cc] = s;
        }
      }
      /* stress resultants from the material that PRODUCED this solve — kept
       * in equilibrium with u even when the loop exits at maxIter (the damped
       * update below is for the NEXT pass only) */
      {
        const D = Dm.subarray(e * 9, e * 9 + 9);
        stresses[e * 3]     = D[0] * exx + D[1] * eyy + D[2] * gxy;
        stresses[e * 3 + 1] = D[3] * exx + D[4] * eyy + D[5] * gxy;
        stresses[e * 3 + 2] = D[6] * exx + D[7] * eyy + D[8] * gxy;
      }
      /* damped update to avoid state chatter: full steps first, then heavier
       * damping so the material field settles instead of oscillating */
      const alpha = iter < 3 ? 1.0 : (iter < 10 ? 0.5 : 0.25);
      for (let r = 0; r < 3; r++) for (let cc = 0; cc < 3; cc++) {
        const d = Dt[r][cc] - Dm[e * 9 + r * 3 + cc];
        Dm[e * 9 + r * 3 + cc] += alpha * d;
        if (Math.abs(d) > maxDD) maxDD = Math.abs(d);
      }
    }

    lastChanged = changed; iters = iter + 1;
    if (iter > 4 && changed === 0 && maxDD < 5e-3) break;
  }

  /* wrinkle directions (major principal strain) for display */
  const thetaW = new Float64Array(ne);
  for (let e = 0; e < ne; e++)
    thetaW[e] = 0.5 * Math.atan2(strains[e * 3 + 2], strains[e * 3] - strains[e * 3 + 1]);
  return { nhat: stresses, uhat: u, state, thetaW, tabNodes, cordT, iters, lastChanged, maxDD };
}

/* ============================================================================
 * OUT-OF-PLANE: div(nhat grad what) = -rhs, given a per-element tension field.
 * bc: 'corners' (pinned at the 4 boom attachments — the truth model) or
 *     'edges'   (classical drumhead — the uniform-tension reference).
 * opts.cordT  — edge-cord tension: adds tensioned-string stiffness T·w'' along
 *               the boundary (the catenary cord also carries transverse load).
 * opts.load   — 'uniform' (gravity/solar pressure, rhs = 1) or
 *               'slew' (rigid slew angular acceleration about the vertical
 *               centerline: rhs = xi - 1/2, antisymmetric billow).
 * Principal tensions are floored at nFloor so wrinkled/slack regions stay
 * SPD; physically the wrinkle direction still carries tension, and the
 * cross direction's true stiffness is small but nonzero.
 * ==========================================================================*/
function solveOutOfPlane(mesh, nfield, bc, nFloor, tabNodes, opts = {}) {
  const { nn, h, nodes, ne, elems } = mesh;
  const slew = opts.load === 'slew';
  const bw = nn + 1;
  const B = new Float64Array(nodes * (bw + 1));
  const f = new Float64Array(nodes);
  const detJw = (h / 2) * (h / 2);

  for (let e = 0; e < ne; e++) {
    /* floor principal tensions */
    let nxx = nfield[e * 3], nyy = nfield[e * 3 + 1], nxy = nfield[e * 3 + 2];
    const m = (nxx + nyy) / 2, r = Math.hypot((nxx - nyy) / 2, nxy);
    let p1 = m + r, p2 = m - r;
    const th = 0.5 * Math.atan2(2 * nxy, nxx - nyy);
    if (p1 < nFloor || p2 < nFloor) {
      p1 = Math.max(p1, nFloor); p2 = Math.max(p2, nFloor);
      const cs = Math.cos(th), sn = Math.sin(th);
      nxx = p1 * cs * cs + p2 * sn * sn;
      nyy = p1 * sn * sn + p2 * cs * cs;
      nxy = (p1 - p2) * cs * sn;
    }
    const Ke = new Float64Array(16), fe = new Float64Array(4);
    const ei = e % (nn - 1);                          // element column, for the slew arm
    for (const [xi, eta] of GP) {
      const g = shapeGrad(xi, eta, h);
      const rhs = slew ? ((ei + (xi + 1) / 2) * h - 0.5) : 1;
      for (let a = 0; a < 4; a++) {
        for (let b = 0; b < 4; b++)
          Ke[a * 4 + b] += (g.dx[a] * (nxx * g.dx[b] + nxy * g.dy[b]) +
                            g.dy[a] * (nxy * g.dx[b] + nyy * g.dy[b])) * detJw;
        fe[a] += g.N[a] * rhs * detJw;
      }
    }
    for (let a = 0; a < 4; a++) {
      const gi = elems[e * 4 + a];
      f[gi] += fe[a];
      for (let b = 0; b < 4; b++) {
        const gj = elems[e * 4 + b];
        if (gj <= gi) B[gi * (bw + 1) + (gi - gj)] += Ke[a * 4 + b];
      }
    }
  }

  /* edge cords: tensioned strings along the boundary, pinned at the tabs */
  if (opts.cordT > 0) {
    const k = opts.cordT / h;
    const seg = (n1, n2) => {
      B[n1 * (bw + 1)] += k; B[n2 * (bw + 1)] += k;
      const hi = Math.max(n1, n2), lo = Math.min(n1, n2);
      B[hi * (bw + 1) + (hi - lo)] -= k;
    };
    for (let i = 0; i < nn - 1; i++) {
      seg(i, i + 1);                                   // bottom
      seg((nn - 1) * nn + i, (nn - 1) * nn + i + 1);   // top
      seg(i * nn, (i + 1) * nn);                       // left
      seg(i * nn + nn - 1, (i + 1) * nn + nn - 1);     // right
    }
  }

  /* boundary conditions: the reinforced corner tabs are held flat by the
   * booms (whole patch pinned); the reference drumhead pins all edges */
  const pins = [];
  if (bc === 'corners')
    pins.push(...(tabNodes ?? [0, nn - 1, nn * (nn - 1), nn * nn - 1]));
  else for (let i = 0; i < nn; i++)
    pins.push(i, nn * (nn - 1) + i, i * nn, i * nn + nn - 1);
  for (const d of pins) { B[d * (bw + 1)] += 1e8; f[d] = 0; }

  bandFactor(B, nodes, bw);
  return bandSolve(B, nodes, bw, f);                  // what >= 0 = sag magnitude
}

/* ---------- field metrics: RMS slope coeff & P-V coeff (dimensionless) --- */
function fieldMetrics(mesh, what) {
  const { nn, h, ne, elems } = mesh;
  let s2 = 0, wmax = 0, wmin = Infinity;
  for (let i = 0; i < nn * nn; i++) {
    if (what[i] > wmax) wmax = what[i];
    if (what[i] < wmin) wmin = what[i];
  }
  const detJw = (h / 2) * (h / 2);
  let area = 0;
  for (let e = 0; e < ne; e++) {
    for (const [xi, eta] of GP) {
      const g = shapeGrad(xi, eta, h);
      let wx = 0, wy = 0;
      for (let a = 0; a < 4; a++) {
        wx += g.dx[a] * what[elems[e * 4 + a]];
        wy += g.dy[a] * what[elems[e * 4 + a]];
      }
      s2 += (wx * wx + wy * wy) * detJw;
      area += detJw;
    }
  }
  const mid = (nn - 1) / 2;
  return {
    rmsSlopeCoeff: Math.sqrt(s2 / area),              // S: slope_rms = S · sigma g L^2 / F
    pvCoeff: wmax - wmin,                             // PV: sag_pv  = PV · sigma g L^3 / F
    centerCoeff: what[mid * nn + mid],
  };
}

/* per-element |grad what| at element centers (for slope shading) */
function elementSlopes(mesh, what) {
  const { nn, h, ne, elems } = mesh;
  const out = new Float64Array(ne);
  for (let e = 0; e < ne; e++) {
    const g = shapeGrad(0, 0, h);
    let wx = 0, wy = 0;
    for (let a = 0; a < 4; a++) {
      wx += g.dx[a] * what[elems[e * 4 + a]];
      wy += g.dy[a] * what[elems[e * 4 + a]];
    }
    out[e] = Math.hypot(wx, wy);
  }
  return out;
}

/* ============================================================================
 * Top level: build the two dimensionless shapes once.
 * truth  — corner-loaded wrinkling tension field, corner-pinned sag
 * ref    — uniform tension nhat = I/sqrt(2), edge-fixed drumhead
 * ==========================================================================*/
function buildShapes(nn, nu, opts = {}) {
  const mesh = makeMesh(nn);
  const ip = solveInPlane(mesh, nu, opts);
  const nFloor = (opts.floorFrac ?? 0.02) * Math.SQRT1_2;
  const oop = { cordT: ip.cordT };
  const truth = solveOutOfPlane(mesh, ip.nhat, 'corners', nFloor, ip.tabNodes, oop);
  const slew = solveOutOfPlane(mesh, ip.nhat, 'corners', nFloor, ip.tabNodes,
                               { ...oop, load: 'slew' });
  const nuni = new Float64Array(mesh.ne * 3);
  for (let e = 0; e < mesh.ne; e++) {
    nuni[e * 3] = Math.SQRT1_2; nuni[e * 3 + 1] = Math.SQRT1_2;
  }
  const ref = solveOutOfPlane(mesh, nuni, 'edges', 0);
  return {
    mesh, inPlane: ip,
    truth: { what: truth, slopes: elementSlopes(mesh, truth), ...fieldMetrics(mesh, truth) },
    slew:  { what: slew,  ...fieldMetrics(mesh, slew) },
    ref:   { what: ref,   slopes: elementSlopes(mesh, ref),   ...fieldMetrics(mesh, ref) },
  };
}

/* ---------- analytic (uniform-tension, equivalent-circle) relations ------ */
const analytic = {
  sigma:   (mass, L) => mass / (L * L),                       // kg/m^2
  thick:   (mass, rho, L) => mass / (rho * L * L),            // m
  tension: (F, L) => F / (Math.SQRT2 * L),                    // N/m
  pressure:(sigma, g) => sigma * g,                           // Pa
  Rm:      (L) => L / Math.sqrt(Math.PI),                     // m
  rmsSlope:(p, L, N) => p * (L / Math.sqrt(Math.PI)) / (2 * Math.SQRT2 * N),
  sagC:    (p, L, N) => p * L * L / (4 * Math.PI * N),        // m
  Freq:    (sigma, g, L, target) => sigma * g * L * L / (2 * Math.sqrt(Math.PI) * target),
};

if (typeof module !== 'undefined') {
  module.exports = { makeMesh, solveInPlane, solveOutOfPlane, fieldMetrics,
                     elementSlopes, buildShapes, analytic };
}

