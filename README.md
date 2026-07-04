# Membrane Mirror — Gravity vs. Flatness

Interactive, single-file physics simulator comparing how flat a corner-tensioned
membrane mirror can be pulled **on the ground** (gravity bows it into a shallow
drumhead) versus **in space** (no transverse load — the same tension leaves it
flat). The point it makes: on the ground you spend hundreds of newtons per
corner just fighting gravity; in orbit that load is gone.

The shipped page is **`index.html`** — a single self-contained file (no
backend, no external requests, no browser storage) served straight from
GitHub Pages (Settings → Pages → deploy from branch → `main` / root).

## Repo layout

```
src/template.html     markup + inlined brand mark
src/styles.css        theme variables (light + dark), all styling
src/physics-core.js   the FEM / analytic core (single source of truth)
src/app.js            controls, worker, rendering, charts, animation
build.js              node build.js → inlines src/ into index.html
test/physics.test.js  acceptance tests; also fails if index.html is stale
assets/               original logo SVG files
```

Edit under `src/`, run `node build.js`, commit both the source and the built
`index.html`. The test suite verifies the built page embeds exactly
`src/physics-core.js`.

Light theme is the default; the header button (or `?theme=dark`) flips it.
State is shareable via query string, e.g. `?g=1.62&F=120&phi=15&theme=dark`
(keys: `F g L m tgt E t phi grab slew theme`).

Mass and film thickness are independent inputs: **mass** is the total as-built
assembly (film + tapes + cords) and alone sets the gravity load σ = m/L²;
**thickness** enters only through E·t (in-plane stretch, the stiffening check,
film stress) and never changes the sag shape.

**Areal density** is a linked control, not extra state — σ = m/L² always.
Editing mass or density rewrites the other; dragging side length holds density
fixed and rescales mass, so ground slope follows the textbook σgL²/F ∝ L².
(At a fixed *lump* of mass the L cancels exactly — slope = S·m·g/F — which is
why side length used to look like it "did nothing" on the ground readouts.)

## Physics

Small-slope linear membrane theory (a film carries load only in tension):

```
N ∇²w = −p        p = σ·g,   σ = mass / L²
```

Two models run together and agree in the uniform-tension limit:

- **FEM ("truth")** — plane-stress solve of the in-plane stress field, with a
  Miller–Hedgepeth wrinkling model (a membrane cannot carry compression;
  wrinkled elements collapse to uniaxial stiffness). The out-of-plane sag then
  solves `div(N·grad w) = −p` with that non-uniform tension field, pinned only
  at the corner tabs.
- **Uniform-tension reference ("ideal")** — the closed-form drumhead
  idealization: `N = F/(√2·L)`, `w_center = pL²/4πN`,
  `σ_slope = p·R_m/2√2·N` with `R_m = L/√π`, and
  `F_required = σ·g·L²/(2√π·σ_target)`.

### Catenary edge cords (half-angle control)

At half-angle φ = 0 the boom force enters the film directly at the corner tabs
(corners carry most of the load; free edges scallop and droop). At φ > 0 each
edge is cut as a circular arc of half-angle φ with a cord along it:

- corner balance (two cord ends at 45° ∓ φ off the boom diagonal):
  `T_cord = F / (√2·(cos φ + sin φ))`
- Young–Laplace turns cord tension into outward normal edge traction:
  `q = T/R = 2·T·sin φ / L`, and the interior tension field becomes (nearly)
  uniform at `N ≈ q` — the FEM reproduces this closed form exactly.
- the cord also carries transverse load as a tensioned string along the edge.

Consequences the sim shows honestly: at **φ = 45°** the cords reproduce the
conventional `N = F/(√2·L)` exactly; at shallow angles most of the corner force
circulates in the cord and only the `sin φ` fraction tensions the film — a 5°
catenary is dramatically floppier than direct corner attach at the same pull
(shallow catenaries "too flat to pull out wrinkles"). The physics runs on the
square domain (the few-percent scallop sagitta `L·tan(φ/2)/2` is drawn, not
meshed), and the boom force is assumed to route entirely through the cords.

### Slew billow

Repointing in orbit: an overhead tracking pass at line-of-sight rate ω peaks at
angular acceleration α ≈ 0.65·ω². The tangential acceleration α·x loads the
film antisymmetrically about the slew axis (`rhs = ξ − ½`, one extra solve of
the same operator); its RMS slope adds to the orthogonal solar-pressure billow
in quadrature. The bottom chart sweeps ω 0–5 °/s; the animated preview rocks
the mirror through a time-compressed repoint cycle with the billow drawn
(hugely exaggerated) at an amplitude that honestly scales with rate².

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

### Corner attachment (grab width)

Load enters the film over a finite patch — physically the little bar/tab the
film bonds to, pulled through a self-centering pulley (whose equal split into
the two edges is exactly the model's symmetry assumption). The **grab width**
control sets that patch as a percentage of the side (default 3 %; the note
under it converts to mm at the current size — ~60 mm on a 2 m mini). A point
attachment would also be numerically ill-posed (corner slopes diverge with
mesh refinement), so the finite patch is physical and regularizing at once.
With cords (φ > 0) the corner detail largely washes out; at φ = 0 it matters.

### Defaults

L = 15 m, total mass = 1000 g (σ ≈ 4.44 g/m²), film t = 2.6 µm, F = 250 N/corner, φ = 10°,
grab 3 %, target 2.3 mrad, E = 4.4 GPa (effective coated-stack modulus; bare
CP1 is ~2.1). The classic sanity numbers (σ ≈ 1.56 g/m², N ≈ 11.8 N/m,
3.9 mrad, 23 mm, 970 N @ 1 mrad) belong to the ideal-uniform reference at
this spec point and are pinned by tests.

### Out of scope (hooks noted in source)

Thermal wrinkling, boom buckling, cord mass/elasticity, direct-attach + cord
hybrid load paths.

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
