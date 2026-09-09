# Is the CH4 model linear?

Short answer: **almost, but not exactly, and the distinction matters.**

* The model has no chemistry, no deposition and no concentration clipping, so
  CH4 is a passive conserved tracer and the model is *mathematically* an
  **affine** map of its inputs — `c = M·E + b(IC, BC)`, not `c = M·E`.
* Scaling emissions **and** initial/boundary conditions together by a constant
  scales the output by exactly that constant, bit for bit.
* But holding IC/BC fixed and varying emissions alone, the map is **not**
  additive. The PPM monotonicity limiter in the advection scheme is the sole
  nonlinearity, and it costs roughly **1% rms / a few % peak**, at every
  perturbation amplitude — the relative error does *not* shrink as the
  perturbation gets smaller.
* `ADJOINT_BWD` is the exact discrete adjoint at the forward state, so it *is*
  exactly linear in the adjoint forcing, and its gradient agrees with forward
  finite differences to ~0.1% for distributed perturbations.

So "the model is linear" is a good ~1% working approximation and a correct
statement about the *physics* we configured, but it is false as an exact
statement about the discretisation. Any argument that needs exact linearity
(a reusable Jacobian, superposition of source regions, scaling emissions and
expecting concentrations to scale) inherits a ~1% error, and scaling arguments
additionally need the IC/BC background subtracted first.

## Why the model is linear apart from one thing

`sciproc.F` in `driver/zenodo_cadj_fwd` calls, per sync step:

```
VDIFF → COUPLE → HADV → ZADV → HDIFF → DECOUPLE → CLDPROC → CHEM → AERO
```

`CLDPROC`, `CHEM`, `AERO`, `PHOT` and `AERO_DEPV` are all no-op stubs in this
build (see the module table in `CLAUDE.md`), and the CH4-only mechanism
declares `N_GC_DEPV = N_GC_DDEP = N_GC_SCAV = N_GC_WDEP = 0`, so the
deposition branches in `vdiffacm2.F` — including the `PLDV_HONO / CONC_NO2`
divisions, the only other concentration-dependent expressions in the file —
are never entered. There is also no positivity clipping anywhere on the CH4
path: the `CMIN` constants in `rdbcon.F` and `initscen.F` are fill values for
*missing* input species, not floors applied to computed concentrations.

That leaves four active operators, and three of them are exactly linear in the
CH4 field with coefficients that depend only on the meteorology:

* **`VDIFF`** (`vdiff/zenodo_fwd/vdiffacm2.F`) — ACM2 assembles a tridiagonal
  system whose coefficients come from `EDDYX` (met only; it is not passed
  `CGRID`), and emissions enter as an additive source term `EMIS = VDEMIS *
  DTS`. Linear in the state, affine in the emissions.
* **`HDIFF`** (`hdiff/zenodo_fwd/hdiff.F`) — explicit diffusion with `K` from
  wind deformation (`hcdiff3d.F`). Exactly linear, and its stencil is
  self-adjoint.
* **`COUPLE` / `DECOUPLE`** — multiply by air density × Jacobian. A diagonal
  scaling.

Sub-stepping is also concentration-independent: `advstep.F` and `syncstep.F`
never touch `CGRID`, and `HDIFF`'s `NSTEPS` comes from `hcdiff3d`'s stability
limit.

Yamartino's mass-conservation adjustment does not break this either.
`hadv.F` computes `ADJFAC` as the ratio of the air density before and after
the first advection direction, i.e. purely from the `ASPC` (density) slot; and
in `zadvyppm.F` the `FLX`-based velocity correction inside `VPPM` is applied
only when advecting the density `RJT`. Species are advected with
`CALL VPPM(NLAYS, DT(SS), DS, VEL(:,SS), CGRID(COL,ROW,:,A2C))` — no `FLX`
argument, so the velocities and sub-step count are fixed by the air density
and are identical across emission scenarios.

### The one nonlinearity: the PPM monotonicity limiter

`hadv/zenodo_dadj_fwd/hppm.F` and `vadv/zenodo_dadj_fwd/vppm.F` implement
Colella–Woodward PPM, whose reconstruction branches on the *tracer field
itself*:

```fortran
C0 = CON( I,S )   - CON( I-1,S )
C1 = CON( I+1,S ) - CON( I,S )
DC( I,S ) = 0.5 * ( C0 + C1 )
IF ( C0 * C1 .GT. 0.0 ) THEN                     ! van Leer slope limiter
   DC( I,S ) = SIGN( 1.0, DC( I,S ) )
             * MIN( ABS( DC( I,S ) ), 2.0*ABS( C0 ), 2.0*ABS( C1 ) )
ELSE
   DC( I,S ) = 0.0
END IF
...
IF ( ( CR( I ) - CON( I,S ) ) * ( CON( I,S ) - CL( I ) ) .GT. 0.0 ) THEN
   ... overshoot correction ...
ELSE                   ! local extremum: reconstruction flattened to a constant
   CL( I ) = CON( I,S )
   CR( I ) = CL( I )
END IF
```

These branches are *positively homogeneous* — every condition is a sign test
on a product of differences, and the `MIN`/`SIGN` construction scales linearly
— so `HPPM(λc) = λ·HPPM(c)` for `λ > 0`. They are **not additive**: which
branch a cell takes depends on the shape of the field, so
`HPPM(c₁+c₂) ≠ HPPM(c₁) + HPPM(c₂)`.

Two consequences worth internalising:

* The nonlinearity is a *kink*, not curvature. Difference quotients converge
  sub-linearly (`~d^0.6` measured), so neither a smaller step nor a central
  difference cleans it up the way it would for a smooth function.
* Because the branch conditions involve only *differences*, a spatially
  uniform background cancels out of them entirely — the non-additivity is
  proportional to the perturbation, not to the 1.8 ppm background, and the
  *relative* error is therefore amplitude-independent. Verified directly: with
  a standalone transcription of `HPPM`, absolute non-additivity was unchanged
  when the pedestal was varied over `1.77 → 0.885 → 0.177 → 0`, and scaled
  linearly with the perturbation amplitude.

### Why the backward model is exactly linear anyway

`ADJOINT_BWD` uses `hadv/zenodo_dadj_bwd` and `vadv/zenodo_dadj_bwd` — the
**discrete** adjoint (`dadj`), not the continuous one. `hadv_bwd.F` reads the
forward concentrations back from `ADJ_HADV_CHK`, replays the forward advection
sub-steps to reconstruct the intermediate states, and then calls
`HPPM_BWD(ni, n_spc_adv, con, conb, ...)` — Tapenade-generated code that takes
the limiter branches from `con`, the *forward* field, and propagates `conb`
through the transposed linearisation. `vppm_bwd.F`'s hand-written `ppm_bwd`
does the same, recomputing the branch tests from the forward field `cn`.

So the branch pattern is fixed by the forward trajectory and is independent of
the adjoint forcing. The backward model is therefore a genuine matrix-transpose
product: exactly linear in the forcing, and the exact gradient of the
linearised forward model.

Do not read the `hadv/yamo_cadj_bwd` sibling and conclude otherwise — it is a
*continuous* adjoint that applies the forward limited `HPPM` to the adjoint
field itself, and it is not part of this build.

## Measurements

All figures below are from the `tests/test-data` smoke-test bundle:
`au-test_v1`, 10 × 10 × 32, 2022-12-07, 24 h, hourly `CONC` output. `R(E)`
denotes the emission response `c(E) − c(0)`, obtained by differencing against
a zero-emission run. The prior in that bundle is large — peak response
0.399 ppmV on a 1.49–1.80 ppmV background, i.e. 23% — so amplitude-scaled
variants (×0.01, ×0.001) were run to check that the conclusions are not an
artefact of that.

The float32 quantisation of the `CONC` file (1.2 × 10⁻⁷ ppmV at 1.8 ppm)
contributes ≤ 0.1% to every figure below; it was budgeted explicitly rather
than assumed negligible.

### Homogeneity holds; additivity does not

| Test | max dev | rms dev |
| --- | --- | --- |
| `c(2E, 2·IC, 2·BC)` vs `2·c(E, IC, BC)` | **0** (bit-exact) | **0** |
| `R(0.5E)` vs `0.5·R(E)`, IC/BC fixed | 5.0% | 2.9% |
| `R(2E)` vs `2·R(E)` | 4.0% | 2.4% |
| `R(10E)` vs `10·R(E)` | 6.2% | 4.6% |

Deviations are relative to `max|reference|` and `rms|reference|` respectively.

### Non-additivity across source regions, by amplitude

Splitting the prior into two disjoint halves (`west`, `east`) and comparing
`R(west+east)` with `R(west) + R(east)`. The rms column is over cells carrying
> 5% of the peak response, so it is not diluted by empty cells:

| Amplitude | max dev | rms dev (signal cells) |
| --- | --- | --- |
| ×1 (peak 0.40 ppmV) | 2.60% | 1.61% |
| ×0.01 (peak 3.8 × 10⁻³ ppmV) | 1.62% | 0.54% |
| ×0.001 (peak 3.7 × 10⁻⁴ ppmV) | 1.28% | 0.63% |

**The relative error does not decay with amplitude.** The absolute residual
scales as `amplitude^1.12` — very close to linear — which is what rules out
the float32 floor as the explanation.

A smooth partition (a linear ramp in `w` and `1−w` rather than a step) is
*worse*, not better — 2.15% rms at full amplitude, 3.10% at ×0.01 — because
both halves then seed plumes across the whole domain.

Column means and domain totals are progressively more forgiving: for the
sharp split at full amplitude, 1.49% rms on the column mean and 0.09% on the
domain total per hour.

### Differentiability

Directional derivative about the prior, `g(d) = [c(E + d·E_west) − c(E)] / d`:

| `d` | `‖g(d) − g(d/3)‖_rms` |
| --- | --- |
| 0.3 → 0.1 | 1.6 × 10⁻⁴ |
| 0.1 → 0.03 | 7.0 × 10⁻⁵ |
| 0.03 → 0.01 | 4.9 × 10⁻⁵ |
| 0.01 → 0.003 | 9.8 × 10⁻⁵ ← float32 noise takes over |
| 0.003 → 0.001 | 2.0 × 10⁻⁴ |

The implied convergence order is ~0.6, i.e. **sub-linear**, and it stalls at
~0.4% rms around `d ≈ 0.01–0.03` before the `O(1/d)` quantisation noise
dominates. One-sided/central and forward/backward asymmetries decay at the
same sub-linear rate. This is the kink signature, and it is why the gradient
test in `openmethane`'s `tests/integration/fourdvar/test_grad_cmaq.py` has an
accuracy floor no choice of epsilon can beat.

### The adjoint is exactly linear in its forcing

Comparing `sens(F_A + F_B)` against `sens(F_A) + sens(F_B)`, and `sens(λF)`
against `λ·sens(F)`, on both `ADJ_LGRID` and `ADJ_LGRID_EM`:

| Forcing pair | max dev |
| --- | --- |
| `2F` and `10F` vs scaled `F` | ≤ 1 × 10⁻⁸ relative |
| `−F` vs `−sens(F)` | exactly 0 |
| two isolated cells, same sign | ≤ 4 × 10⁻⁸ relative |
| adjacent cells, opposite signs | ≤ 4 × 10⁻⁸ relative |
| uniform layer + embedded spike | ≤ 5 × 10⁻⁸ relative |
| smooth Gaussian bumps, all hours | ≤ 3 × 10⁻⁶ relative |
| **two independent random full-3D fields** | 8 × 10⁻⁶ relative |

The random-field case is the decisive one — it is maximally adversarial for a
limiter, and additivity still holds to the float32 rounding of the forcing
file itself.

### Adjoint gradient vs forward finite difference

Cost `J = Σ` layer-1 CH4 over the domain at the final hour; perturbation
distributed over all cells of layer 1 at every hour:

| `eps` | one-sided FD / adjoint | central FD / adjoint |
| --- | --- | --- |
| 0.3 | 0.999820 | 0.999761 |
| 0.03 | 1.003012 | 1.001205 |
| 0.003 | 1.014053 | 1.001266 |

0.02–0.13% agreement for a well-conditioned distributed perturbation, even at
`eps = 0.3`. The degradation at `eps = 0.003` is the quantisation floor, not
the model. This is consistent with the numbers recorded in `openmethane`'s
gradient test, which also measured 2–4% errors for *single-cell* perturbations
— the sharp-gradient case where the limiter is permanently engaged.

## What to say, and what not to say

Safe:

* CH4 here is a passive conserved tracer; there is no chemical or depositional
  nonlinearity, by construction.
* The model is affine in `(emissions, IC, BC)` and exactly homogeneous when all
  three are scaled together.
* The adjoint is the exact discrete adjoint at the forward state, so `∇J` is a
  true linear functional of the observation residuals.
* Superposition of source regions, and treating the emission–concentration
  relationship as a fixed Jacobian, are good to ~1% rms (a few % peak) — fine
  for most purposes, and worth stating as an approximation with that number
  attached.

Not safe:

* "The model is linear, therefore `c(E₁+E₂) = c(E₁) + c(E₂)`" — it is affine,
  so the background has to come out first, and even then additivity is only
  ~1%.
* "The response is exactly proportional to emissions" — only if IC and BC are
  scaled too.
* "Small perturbations are in the linear regime" — the *relative* error is
  amplitude-independent; shrinking the perturbation does not buy linearity, it
  only buys a worse signal-to-quantisation ratio.
* "The model is differentiable, so finite differences converge" — the limiter
  is a kink; convergence is sub-linear and floors around 0.4%.

## Reproducing this

Every measurement above came from `ADJOINT_FWD` / `ADJOINT_BWD` runs in the
`cmaq-adjoint` image against the checked-in test bundle, ~4 s per forward run
on the `au-test_v1` domain. The pattern:

```shell
# one scenario per directory, each holding its own emis.nc / FORCE.nc
docker run --rm \
  -v "$PWD/tests/test-data:/opt/data:ro" -v "$SCRATCH:/opt/work" \
  cmaq-adjoint /opt/work/runall.sh base zero west east
```

where `runall.sh` sources `/opt/data/fwd.env`, overrides `EMIS_1` (and
`INIT_*_1` / `BNDY_*_1` for the IC/BC scaling test) to point into the scenario
directory, and runs the binary in a fresh `RUN_DIR`. Backward runs need a
forward run in the same directory first, then `bwd.env` with `ADJ_FORCE`
overridden. `netCDF4` is not in the runtime image, so build the scenario
inputs and difference the outputs on the host.

Two traps that cost time here:

* When building a partitioned forcing, **add** to the field rather than
  assigning into it. Overwriting one cell of a broad field silently changes
  the total, and the resulting "non-additivity" is just the missing component.
  A deviation that is suspiciously identical across two very different
  amplitudes is the tell.
* Budget the float32 quantisation of `CONC` explicitly before believing any
  small-amplitude result. At ×0.001 amplitude it is 0.1% of the signal; at the
  amplitudes a naive finite-difference test would pick, it can be all of it.
