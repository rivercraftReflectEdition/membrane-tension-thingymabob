# Membrane Mirror — Gravity vs. Flatness

Interactive, single-file physics simulator comparing how flat a corner-tensioned
membrane mirror can be pulled **on the ground** (gravity bows it into a shallow
drumhead) versus **in space** (no transverse load — the same tension leaves it
flat). The point it makes: on the ground you spend hundreds of newtons per
corner just fighting gravity; in orbit that load is gone.

Everything is in **`index.html`** — no build, no backend, no external requests,
no browser storage. Open it locally or serve it from GitHub Pages
(Settings → Pages → deploy from branch → `main` / root).

State is shareable via query string, e.g. `?g=1.62&F=120` (keys: `F g L m tgt E rho`).

## Physics

Small-slope linear membrane theory (a film carries load only in tension):

```
N ∇²w = −p        p = σ·g,   σ = mass / L²
```

Two models run together and agree in the uniform-tension limit:

- **Corner-load FEM ("truth")** — plane-stress solve of the in-plane stress
  field from four outward diagonal corner forces, with a Miller–Hedgepeth
  wrinkling model (a membrane cannot carry compression; wrinkled elements
  collapse to uniaxial stiffness). The out-of-plane sag then solves
  `div(N·grad w) = −p` with that non-uniform tension field, pinned only at the
  corner tabs. Corners carry most of the load; free edges scallop and droop.
- **Uniform-tension reference ("ideal")** — the closed-form drumhead
  idealization: `N = F/(√2·L)`, `w_center = pL²/4πN`,
  `σ_slope = p·R_m/2√2·N` with `R_m = L/√π`, and
  `F_required = σ·g·L²/(2√π·σ_target)`.

Both are solved **once, dimensionless** at startup (in a worker); the stress
field shape depends only on Poisson's ratio, so every control change afterwards
is a pure rescale — `w ∝ σgL³/F`, `slope ∝ σgL²/F`.

### A note on the corner-force → tension mapping

`N = F/(√2·L)` is a convention, not an exact result: **no uniform isotropic
tension field can equilibrate four corner point loads** (a mid-line cut demands
`∫N dx = √2·F/L`; a diagonal cut demands `F/(√2·L)` — they differ by 2×). The
FEM resolves the real field: the interior carries ≈ 1.55× the nominal `N`
(the wrinkled edges shed load inward), so at the same pull the discrete film
sits somewhat flatter than the ideal-uniform reference, and pays for it with
scalloped, drooping free edges. Both numbers are shown; the reference keeps the
conventional formulas above.

In space the only transverse load is solar radiation pressure (9.1 µPa) — the
same solve scaled to it, so the space panel shows real (tiny) numbers and a
real required corner force instead of a bare zero.

The hero panels draw sag at a fixed ×50 vertical exaggeration (smoothly
compressed once the drawn amplitude would exceed ~24 % of the side), so the
membrane flattens continuously as tension rises — no rescale jumps.

### Out of scope in v1 (hooks noted in source)

Thermal wrinkling, boom buckling, catenary edge cords.

## Tests

```
node test/physics.test.js
```

The test extracts the physics core from between the `@physics-core` sentinels
in `index.html` — the tested code is byte-for-byte the shipped code — and
checks the spec sanity numbers (σ ≈ 1.56 g/m², t ≈ 1.1 µm, N ≈ 11.8 N/m,
RMS slope ≈ 3.9 mrad, sag ≈ 23 mm, F@1 mrad ≈ 970 N at the defaults), the
exact square-drumhead coefficient (0.10419), static-equilibrium flux across a
mid-line cut, and the mesh-converged FEM coefficients.
