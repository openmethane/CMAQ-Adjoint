# The vertical advection flux diagnostic

`ADJOINT_FWD` reports what vertical advection moves across each layer
interface. This document describes what is measured, how to switch it on, how
to read it, and what it cannot tell you.

It exists for [openmethane/openmethane#256][256], which asks whether advective
flux through the model top is a material source of the drift measured in #248
and #242. Nothing here changes the solution — see
[Instrumentation, not physics](#instrumentation-not-physics).

[256]: https://github.com/openmethane/openmethane/issues/256

## What is being measured, and why in two directions

CMAQ gives methane no top boundary condition. `zadvyppm.F` carries the comment

```
C No boundary conditions are needed because VEL(1) = VEL(NLAYS+1) = 0
```

and then three lines later sets

```fortran
      VEL( NLAYS+1,SS ) = FLX( NLAYS+1 ) / RJT( NLAYS )
```

which is the column's air-mass budget residual, not zero. Air therefore crosses
the lid in both directions. What it carries is decided in `vppm.F`, and both
directions carry the same value:

- `PPM` assigns the top cell a zeroth-order intercept, `CM(NI+1) = CN(NI)`.
- The monotonicity test `(CR(I) - CN(I)) * (CN(I) - CL(I)) > 0` then fails with
  equality at `I = NI`, so the else branch sets `CL(NI) = CR(NI) = CN(NI)`.
- Cell `NLAYS` is consequently advected as a constant. Outflow through the top
  face and the explicit inflow term `FM(NI+1) = Y * CON(NI)` both carry
  `CON(NLAYS)`, the top layer's own concentration.

So within one call the interface is symmetric, and **the net flux is not the
quantity of interest**. The asymmetry is temporal: methane leaves at one time's
top-layer value and returns at a later one's, while the lateral boundary keeps
restoring the perimeter. A lid that exchanges hard with the two directions
nearly cancelling looks identical to a sealed lid if only the net is recorded.
Both directions are therefore kept apart, for methane and for air.

The air-mass pair gives the ventilation rate directly, which is the quantity
the order-of-magnitude argument in the #242 analysis turns on.

## Switching it on

| Output | Default | Controlled by |
| --- | --- | --- |
| Log lines, every sync step, domain integrated, every interface | always on | nothing |
| Gridded file, every output step, model top only | off | `CTM_VADV_TOPFLX` |

Set `CTM_VADV_TOPFLX` to the path of the file to write. Leave it unset and no
file is opened and no per-column accumulator is allocated.

In the Open Methane repository, `WRITE_VADV_TOPFLX` turns it on;
`cmaq_handle.setup_run` puts the path in the model's environment, and each
forward pass's file and log records are filed under
`output/vadv-topflx/<pass>/` so that the next pass does not overwrite them.

## The log lines

Both are written by rank 0 only, after a global sum across ranks, so one log
file carries the whole domain. Values are in **kg/s**, meaned over the sync
step. The timestamp is the **start** of the sync step.

`VADVTOP` is the headline: the model top interface alone.

```
 VADVTOP 2022341 000000 CH4      up  2.7783E+01 dn  6.4533E+01 net -3.6750E+01 kg/s
 VADVTOP 2022341 000000 AIR      up  3.3430E+07 dn  7.7666E+07 net -4.4236E+07 kg/s
```

`VADVFLX` carries the whole profile: `NLAYS+1` values, interface 1 at the
ground through interface `NLAYS+1` at the model top, one line per direction.

```
 VADVFLX 2022341 000000 CH4      UP kg/s  0.0000E+00  4.6357E+01  1.0347E+02 ...
 VADVFLX 2022341 000000 CH4      DN kg/s  0.0000E+00  4.2000E+01  9.3913E+01 ...
```

Interface 1 is always zero: the surface is impermeable, `VEL(1) = 0`.

Summing `UP - DN` at successive interfaces over a run closes the vertical mass
budget layer by layer, which is what #248 could not do because no fluxes were
archived.

Both records are fixed width, so a parser can take them by column:

| Field | Columns | Meaning |
| --- | --- | --- |
| tag | 2-8 | `VADVTOP` or `VADVFLX` |
| date | 10-16 | model date, `YYYYDDD` |
| time | 18-23 | model time, `HHMMSS`, start of the sync step |
| variable | 25-32 | species name, or `AIR`, blank padded |

On a `VADVTOP` line the three values sit in `1PE12.4` fields at columns 36-47,
51-62 and 67-78. On a `VADVFLX` line, columns 34-35 hold `UP` or `DN` and the
`NLAYS+1` values follow in consecutive `1PE12.4` fields from column 41.

Volume: six lines per sync step for the CH4-only mechanism, measured at 319
bytes a line, so about 17 MB of log for a month at a five-minute sync step and
8.5 MB at ten minutes.

## The gridded file

An I/O API gridded file on the model's own horizontal grid, one layer,
timestepped at the output step. Four variables for the CH4-only mechanism:

| Variable | Units | Meaning |
| --- | --- | --- |
| `CH4_UP` | kg/s | methane carried up through the model top |
| `CH4_DN` | kg/s | methane carried down through the model top |
| `AIR_UP` | kg/s | air carried up through the model top |
| `AIR_DN` | kg/s | air carried down through the model top |

Each is a **mean over the output interval**, per cell, and the timestamp is the
**start** of that interval — the same convention as `ACONC`. Net upward flux is
`_UP` minus `_DN`. Reading it is then:

```python
ds = xr.open_dataset(path, decode_times=False)
net = (ds["CH4_UP"] - ds["CH4_DN"]).sum(("COL", "ROW"))   # kg/s, per hour
gg_per_day = net * 86400.0 / 1.0e6
```

The spatial field is the point of writing a file rather than only a number: the
#255 boundary-condition staircase is injected on the perimeter, so a lid flux
driven by it should show up on the perimeter and not in the interior.

Size is about 2.3 GB a month per forward run on `aust10km`, which is 454 x 430
cells, at hourly output. Four variables of four bytes over 744 hours.

## How the mass conversion works

Advection runs on coupled concentrations. `couple/zenodo/couple.F` multiplies
gas, non-reactive and tracer species by `RHOJ`, and `rdbcon.F` records what
`RHOJ` is:

```
C   SqRDMT = Sq. Root [det ( metric tensor )]
C          = Vertical Jacobian / (map scale factor)**2
C   Air Density X SqRDMT = RHOJ
```

`VPPM` returns transports in units of the advected quantity times the sigma
thickness swept during the timestep. Because `RHOJ` already carries the map
scale factor, multiplying by the *nominal* cell area — and by nothing else —
gives an extensive quantity:

| Variable | Accumulator units | kg per unit |
| --- | --- | --- |
| gas, non-reactive, tracer | ppmV x kg air / m**2 | `1e-6 * XCELL * YCELL * MW / MWAIR` |
| air | kg / m**2 | `XCELL * YCELL` |

Two checks on that, both on the test bundle's 10 x 10 domain:

- `sum(DENSA_J * dsigma)` over the column is 9,936 kg/m**2 against a hydrostatic
  `(ps - ptop) / g` of 9,553 kg/m**2. The 4% gap is `1 / m**2` at 34S, which is
  the map scale factor the formula above relies on being already present.
- `CH4_DN / AIR_DN` is 8.32e-7, a mass mixing ratio, which is 1.51 ppmV. That is
  the top layer's own concentration, as the code says it must be.

Aerosol species are coupled through `JACOBM` in mass or number per volume
rather than through `RHOJ`, so this conversion does not apply to them. The
CH4-only mechanism this fork builds has none; if any were ever added they are
excluded from both outputs and a warning is logged, rather than being converted
wrongly.

## Instrumentation, not physics

`FLXUP` and `FLXDN` are optional `INTENT(OUT)` arguments of `VPPM`, assigned
immediately before it returns. Nothing reads them back, so requesting them
cannot perturb the solution.

Verified rather than asserted: a forward day run on the test bundle before and
after this change produces `CONC`, `ACONC`, `CGRID` and the vertical advection
checkpoint in which **every data value is bit-identical**. The files differ by
two bytes, which are the I/O API `CTIME` and `WTIME` wall-clock stamps.

## What it does not tell you

- **Only the advective term.** Vertical diffusion is sealed at the lid
  (`eddyx.F` sets `EDDYV(NLAYS) = 0`), so there is nothing to measure there, but
  horizontal advection and diffusion through the lateral boundary are not
  instrumented and the domain remains open. The vertical budget closes; the
  whole-domain budget still does not.
- **Nothing about whether the flux is wrong.** A lid flux is expected in any
  open domain. The question is whether it is large against the 0.22 Gg/day the
  #248 mass budget attributes to layers 27-31, and whether it falls when the
  #255 interpolation replaces the boundary-condition staircase.
- **The test-bundle numbers are not transferable.** That domain is 10 x 10 cells
  of 10 km. Column-integrated horizontal divergence across 100 km is a large
  fraction of the column, so the diagnosed lid velocity there is about 11 cm/s
  and the net air flux is 3.8% of the domain's air mass per day. Those figures
  say the instrument reads sensibly, and nothing about `aust10km`.

## Cost

Per sync step the reduction is one scalar `SUBST_GLOBAL_SUM` per interface per
direction per variable — 132 reductions for the CH4-only mechanism at 32
layers. The stenex no-op linked into the serial build has no array form, so
this cannot be a single call. The accumulation itself is one pass over
`NLAYS+1` per column per sub-step per advected variable.

Serial and two-rank runs of the test bundle agree on all 384 `VADVTOP` records
except two, which differ in the fifth significant figure — the summation order
of the reduction.
