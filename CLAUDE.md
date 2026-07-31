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

Each has a block of `set Mod<Component> = <component>/<variant>` lines, e.g.
(bwd build, as of this writing):

```
set ModDriver = driver/yamo_adj_bwd
set ModChkpnt = chkpnt/chkpnt_ioapi
set ModInit   = init/adj_bwd
set ModCpl    = couple/gencoor
set ModHadv   = hadv/yamo_cadj_bwd
set ModVadv   = vadv/vyamo_cadj_bwd
set ModHdiff  = hdiff/multiscale_adj_bwd
set ModVdiff  = vdiff/acm2_inline_adj_bwd
set ModPhot   = phot/phot_noop
set ModChem   = chem/chem_bwd_noop
set ModAero   = aero/aero_bwd_noop
set ModAdepv  = aero_depv/aero_depv_noop
set ModCloud  = cloud/cloud_bwd_noop
set ModPa     = procan/pa
set ModUtil   = util/util_adj
```

**Before reading or editing any file under `cmaq/CCTM/`, check these two
scripts first** to confirm the module is actually part of the active build
for the binary you care about (fwd vs bwd) — otherwise you may be looking at
an unused legacy variant that has no effect on program behavior. The fwd and
bwd scripts can (and do) select different variants of the same component.

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
