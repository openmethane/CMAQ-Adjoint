# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project overview

CMAQ forward and backward adjoint model, forked from
[CMAQ 5.0 Adjoint](https://zenodo.org/records/3780216) (University of Colorado
Boulder development repo), modified to support CH4 modelling for the Open
Methane project. The bulk of the repository (`cmaq/`) is legacy CCTM-era
Fortran; there is no application-level Python/JS code, only build tooling.

## Commands

All building and testing is done through Docker — there is no supported
native build path in this repo (see "Building outside Docker" below).

```shell
make build          # build the docker image, tagged `cmaq-adjoint` (linux/amd64)
make build-aarch64   # same, for arm64
make run             # build then run the container, mounting the repo at /opt/project
make test             # build then run tests/test-run.sh inside the container
```

Run the compiled binaries directly:

```shell
docker run -it --rm cmaq-adjoint /opt/cmaq/bin/ADJOINT_FWD
docker run -it --rm cmaq-adjoint /opt/cmaq/bin/ADJOINT_BWD
```

Both binaries require input data and environment variables set by the caller
(not provided by this repo/image).

### Debugging the build

The Dockerfile is a multi-stage build; the `builder` stage (which has the
compiled source, libraries and build artifacts under `/opt/cmaq`) is discarded
in the final image, which makes compile failures hard to inspect. To debug:

```shell
docker build . --progress=plain --target builder -t cmaq-adjoint-builder
docker run -it --rm cmaq-adjoint-builder bash
```

### Testing

`make test` (or the `Run tests` CI step) mounts `tests/` at `/opt/tests` and
runs `/opt/tests/test-run.sh`, which:

1. Downloads and caches a test data bundle (MCIP sample output + matching
   initial/boundary conditions) from the Open Methane public data store into
   `tests/test-data/` (too large to commit).
2. Runs `ADJOINT_FWD` once single-process, then again under
   `mpirun -np 2` with `NPCOL_NPROW="2 1"`.

**This is a smoke test only** — it checks the binary compiles and runs to
completion on known inputs, not that the numerical output is correct. There is
no automated numerical regression test in this repo.

### Versioning and changelog

Version is managed by `uv version` (see `pyproject.toml`); changelog is built
by `towncrier` from fragment files in `changelog/`, named `<PR#>.<type>.md`
(`feature`/`improvement`/`fix`/`docs`/`deprecation`/`breaking`/`trivial` — see
`changelog/README.md`). Releases are cut via the manual `release.yaml`
GitHub Actions workflow (`workflow_dispatch`), which bumps the version,
runs `towncrier build`, tags, creates a GitHub release, then bumps to the next
`dev` pre-release version. Don't hand-edit `CHANGELOG.md` directly for a PR —
add a fragment file instead.

## Architecture

### The module-selection build scripts are the map of this codebase

This is a fork of the classic CMAQ CCTM source layout: `cmaq/CCTM/<component>/`
(e.g. `hadv`, `vadv`, `hdiff`, `vdiff`, `chem`, `aero`, `chkpnt`, `driver`,
`init`, `couple`, ...) each contain **many alternative implementations** as
sibling subdirectories (numerical schemes, adjoint/fwd/bwd/ddm3d variants,
no-op stubs, etc.) — most of these are unused legacy CMAQ science options
carried over from upstream and are dead weight for this project.

Which subdirectory is actually compiled for a given component is decided
entirely by two csh build scripts, not by any config in the source tree
itself:

- `scripts/bldit.adjoint.fwd.openmethane` → builds `ADJOINT_FWD`
- `scripts/bldit.adjoint.bwd.openmethane` → builds `ADJOINT_BWD`

Each has a block of `set Mod<Component> = <component>/<variant>` lines. As of
this writing the two builds select (the `zenodo_*` variants are the ones this
fork actually uses — the `yamo_*` / `acm2_inline_*` / `multiscale_*` siblings
are upstream legacy and are **not** built):

| Component | fwd (`ADJOINT_FWD`) | bwd (`ADJOINT_BWD`) |
| --- | --- | --- |
| `ModDriver` | `driver/zenodo_cadj_fwd` | `driver/zenodo_cadj_bwd` |
| `ModChkpnt` | `chkpnt/chkpnt_zenodo` | `chkpnt/chkpnt_zenodo` |
| `ModInit` | `init/zenodo_cadj_fwd` | `init/zenodo_cadj_bwd` |
| `ModCpl` | `couple/zenodo` | `couple/zenodo` |
| `ModHadv` | `hadv/zenodo_dadj_fwd` | `hadv/zenodo_dadj_bwd` |
| `ModVadv` | `vadv/zenodo_dadj_fwd` | `vadv/zenodo_dadj_bwd` |
| `ModHdiff` | `hdiff/zenodo_fwd` | `hdiff/zenodo_bwd` |
| `ModVdiff` | `vdiff/zenodo_fwd` | `vdiff/zenodo_bwd` |
| `ModPhot` | `phot/phot_noop` | `phot/phot_noop` |
| `ModChem` | `chem/chem_noop` | `chem/chem_bwd_noop` |
| `ModAero` | `aero/aero_noop` | `aero/aero_bwd_noop` |
| `ModAdepv` | `aero_depv/aero_depv_noop` | `aero_depv/aero_depv_noop` |
| `ModCloud` | `cloud/cloud_noop` | `cloud/cloud_bwd_noop` |
| `ModPa` | `procan/pa` | `procan/pa` |
| `ModUtil` | `util/util_adj` | `util/util_adj` |

Mechanism: `mech/cb05cl_ae5_aq_CH4only` — one gas species (`CH4`), no
aerosol or non-reactive species, no dry or wet deposition.

Re-derive the table from the scripts rather than trusting it if the build
behaviour surprises you:

```shell
grep -E "^set Mod|^ *set Mod" scripts/bldit.adjoint.{fwd,bwd}.openmethane
```

**Before reading or editing any file under `cmaq/CCTM/`, check these two
scripts first** to confirm the module is actually part of the active build
for the binary you care about (fwd vs bwd) — otherwise you may be looking at
an unused legacy variant that has no effect on program behavior. The fwd and
bwd scripts can (and do) select different variants of the same component.

### Linearity of the CH4 model

`docs/ch4-linearity.md` settles this with code analysis plus direct
experiment: CH4 is a passive conserved tracer (no chemistry, no deposition, no
clipping), so the model is **affine** in `(emissions, IC, BC)` and exactly
homogeneous when all three are scaled together — but the PPM monotonicity
limiter in `hppm.F`/`vppm.F` makes the emission-to-concentration map
**non-additive at ~1% rms / a few % peak, independent of perturbation
amplitude**. `ADJOINT_BWD` is the exact discrete adjoint at the forward state
and *is* exactly linear in its forcing. Read that document before making any
argument that relies on linearity.

`cmaq/ICL` holds the global includes (species tables, mechanism definitions,
etc.) referenced by `INCLUDE SUBST_*` across the Fortran source, independent
of the per-component module selection above.

### Docker build pipeline

`Dockerfile` uses `ghcr.io/openmethane/cmaq:stable` (the
[openmethane/CMAQ](https://github.com/openmethane/CMAQ) repo) as both the
`builder` base and the final runtime base — that upstream repo is the
reference for the underlying CMAQ build system (`config.cmaq`), and for the
I/O API / netCDF / MPICH toolchain this project assumes is already present.
The `builder` stage:

1. Installs build tooling (gfortran, mpich dev headers, netcdf dev headers, csh, m4).
2. Copies `cmaq/` in as the CCTM model source tree.
3. Runs `bldit.adjoint.fwd.openmethane` then `bldit.adjoint.bwd.openmethane`
   to produce `/opt/cmaq/bin/ADJOINT_FWD` and `/opt/cmaq/bin/ADJOINT_BWD`.

The final image copies only `/opt/cmaq/bin` out of the builder and installs
just the runtime shared libs (`libnetcdff7`, `mpich`) — no compiler toolchain.

CI (`.github/workflows/build_docker.yaml`) builds and pushes this image to
`ghcr.io/openmethane/cmaq-adjoint` on pushes to `main`, on `v*` tags, and on
PRs, and runs the same `tests/test-run.sh` smoke test against the built image.
