# Boundary condition sensitivity (`ADJ_LGRID_BC`)

`ADJOINT_BWD` can write the sensitivity of the cost function to the lateral
boundary condition, alongside the sensitivities to the initial concentration
(`ADJ_LGRID`) and to emissions (`ADJ_LGRID_EM`). A 4D-Var system can use it to
adjust the boundary condition directly, instead of standing in for the
boundary with emissions in the edge cells.

## Enabling it

Set `ADJ_LGRID_BC` to an output path in the backward run's environment:

```shell
ADJ_LGRID_BC="${RUN_DIR}/output/ADJ_LGRID_BC.20221207.nc"
```

If it is unset, nothing is computed or written, and the backward run's other
outputs are bit-for-bit what they were before this feature.

The backward run replays the forward run's horizontal advection with
`BNDY_GASC_1`, so it must point at the same BCON file the forward run used.

## The file

`ADJ_LGRID_BC` has the same structure as the `BNDY_GASC_1` input, so a caller
can line it up with that file record for record and cell for cell.

| property | value |
|---|---|
| file type | `BNDARY3`, `NTHIK = 1` |
| dimensions | `TSTEP × LAY × PERIM`, the same as `BNDY_GASC_1` |
| variables | one per `BNDY_GASC_1` variable, with the same names |
| units | `CF/ppmV`, cost function units per ppmV of boundary concentration |
| time axis | `BNDY_GASC_1`'s own `TSTEP`, from the record at or before the run start to the record at or after the run end |
| perimeter order | the I/O API standard S, E, N, W, each side with its corner |

Each record holds the total sensitivity of the cost to that BCON record. That
total is summed over every advection sub-step whose time interpolation uses
the record.

If `BNDY_GASC_1` is time-independent (`TSTEP = 0`), the output is too. It has
one record holding the sensitivity to that single boundary field over the
whole run.

Some entries are always zero:

- **The four corners.** Horizontal advection never reads them.
- **Outflow faces.** The forward model fills an outflow ghost cell by
  extrapolating from the interior (`ZFDBC`), not from BCON.
- **Variables no advected species reads.** These are kept so that the header
  matches `BNDY_GASC_1`. If several advected species share a BC surrogate,
  their sensitivities are summed into it.

The file is written once, at the end of the backward run. A BCON record is
read by the sub-steps on both sides of it, so its sensitivity is not complete
until the run reaches the record before it. That is unlike `ADJ_LGRID_EM`,
which is written one output step at a time.

## How it is computed

Horizontal advection is the only process in the active build that reads the
boundary condition.

- **Forward.** `XADV`/`YADV` fill the ghost cells on an inflow face with the
  coupled boundary value `BCFAC · BBUF · RHOJ_BND`. `BBUF` is `BNDY_GASC_1`
  interpolated in time to the mid-point of the advection sub-step. At a
  domain-edge PE, `HPPM` reads only the ghost cell next to the face, as the
  inflow flux.
- **Backward.** `HPPM_BWD` already accumulates the adjoint of that flux into
  the ghost cell's adjoint. Then:
  1. `XADV_BWD`/`YADV_BWD` return that ghost adjoint on forward-inflow faces.
     They return it as the change `HPPM_BWD` makes, because the ghost is
     seeded beforehand with a value that never reaches the interior.
  2. `RDBCON` (the backward model's copy) returns the coupling factor
     `BCFAC · RHOJ_BND`. `HADV_BWD` saves it, with each sub-step's start
     time, during its replay of the forward sub-steps. That replay calls
     `RDBCON` with exactly the times the forward run used.
  3. `BC_SENS_ACCUM` (module `ADJ_BC_SENS`) multiplies the ghost adjoint by
     the factor. It splits the result between the two BCON records that
     bracket the sub-step's mid-point, using `INTERPB`'s own linear weights.

Under MPI, each PE maps its local perimeter slots onto the global perimeter
and fills only the part of the boundary it owns. At the end of the run the
pieces are gathered to PE 0, which writes the file with plain I/O API
`WRITE3` (through `write3_bndy.f`). pario's `PWRITE3` only handles gridded
files. On a boundary file it logs a warning and returns failure, so the output
would silently go missing.

The code is in `cmaq/CCTM/hadv/zenodo_dadj_bwd/` (`ADJ_BC_SENS.F`,
`write3_bndy.f`, `hadv_bwd.F`, `xadv_bwd.F`, `yadv_bwd.F` and `rdbcon.F`).
`OP_ADJ_FILE` in `driver/zenodo_cadj_bwd/ADJOINT_FILES.F` creates the file
header. The forward build is untouched.

### Limitations

- Only gas-phase species are covered, because the output mirrors
  `BNDY_GASC_1`. Aerosol, non-reactive and tracer boundary files are not
  covered. The Open Methane CH4-only build advects a single gas species.
- `BNDY_GASC_1` must have `NTHIK = 1` and match the model grid. The run stops
  with an error otherwise.
- The model top has no boundary condition, so there is no top sensitivity to
  report.

## Validation

Measured on the test bundle's 10 × 10 × 32 domain for 2022-12-07. The forcing
is the bundle's `FORCE.20221207.nc`, a single unit cell in the surface layer
at 24:00. The hourly BCON is the bundle's
`cmaq/BCON.d01.au-test_v1.CH4only.hourly.nc`, described in `tests/README.md`.

**Closure.** The cost is linear in concentration. If transport also keeps a
uniform mixing ratio uniform, then
`Σ ADJ_LGRID(t0) + Σ ADJ_LGRID_BC = Σ forcing`:

| BCON | Σ `ADJ_LGRID` (00:00) | Σ `ADJ_LGRID_BC` | total | Σ forcing |
|---|---|---|---|---|
| bundle, single record | 0.0000 | 1.0001 | 1.0001 | 1 |
| hourly | 0.0000 | 1.0001 | 1.0001 | 1 |

After 24 hours on a 100 km domain, the initial condition has been flushed out
and nearly all of the sensitivity sits on the boundary. The sensitivity is on
the north and west faces, the inflow faces for this day's winds. The south and
east faces are within 1e-6 of zero, and the corners are exactly zero.

**Finite differences.** These perturb the hourly BCON on one edge, one band of
layers and one record, by ±δ ppmV. Each row compares the central difference of
the forward cost with `⟨ADJ_LGRID_BC, δ⟩ / δ`:

| case (edge, layers, record) | adjoint | FD, δ = 0.001 | FD, δ = 0.01 | FD, δ = 0.1 |
|---|---|---|---|---|
| N, 1–10, 22:00 | 0.5392 | 0.5419 (0.5%) | 0.5522 (2.3%) | 0.5487 (1.7%) |
| N, 1–10, 23:00 | 0.2028 | 0.2034 (0.3%) | 0.2179 (6.9%) | 0.2205 (8.0%) |
| N, 1–10, 21:00 | 0.1654 | 0.1638 (1.0%) | 0.1486 (11%) | 0.1179 (40%) |
| W, 1–10, 21:00 | 0.0364 | 0.0361 (0.8%) | 0.0356 (2.3%) | 0.0261 (39%) |
| all, all, 21:00 | 0.2199 | 0.2176 (1.1%) | 0.1956 (12%) | 0.1646 (34%) |
| N, 1–10, 24:00 | −0.00417 | −0.00399 (4.4%) | −0.00361 (15%) | −0.00273 (53%) |
| S, 1–10, 21:00 | 0 | 0 | 0 | 0 |

The error falls as δ shrinks, which is what forward nonlinearity looks like
(the PPM limiter), not an adjoint error. A record, edge or layer placed wrongly
would be off by a factor of several. Compare, for example, the 21:00 and 22:00
rows. The small negative sensitivity to the 24:00 record is confirmed by the
finite differences, so it is real model behaviour. One case is left out of
the table: layers 11–21 at 20:00. Its sensitivity (6.6e-5) is below what
float32 `CONC` output can resolve at δ = 0.001.

**Decomposition.** Single process and `mpirun -np 2` agree to a relative
2e-7. That is the same level at which `ADJ_LGRID` itself differs between the
two decompositions.

**Forward and existing outputs unchanged.** `ADJOINT_FWD`'s `CONC`, `ACONC`,
`CGRID`, checkpoints and `XFIRST` are identical to the v2.0.1 release image.
`ADJOINT_BWD`'s `ADJ_LGRID` and `ADJ_LGRID_EM` are also identical to v2.0.1,
with `ADJ_LGRID_BC` set and unset.

`tests/test-run.sh` covers what can be checked without netCDF tools in the
image:

- the backward run with `ADJ_LGRID_BC` set, against both the single-record
  and the hourly BCON;
- that the output is not empty;
- that single-process and two-process runs agree.

The closure and finite-difference checks above were run by hand.
