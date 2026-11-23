# CMAQ Adjoint

CMAQ forward and backward adjoint, based on
[CMAQ 5.0 Adjoint](https://zenodo.org/records/3780216),
forked from the
[development repo](https://adjoint.colorado.edu:8080) at University of
Colorado Boulder.

Additional changes have been made to support CH4 modelling for the Open Methane
project.

## Package

The adjoint is published as one Debian package, `openmethane-cmaq-adjoint`,
attached to each [GitHub release](https://github.com/openmethane/CMAQ-Adjoint/releases)
together with a `SHA256SUMS` file. It is built for Debian 13 (trixie), and its
dependencies name trixie packages.

| Package | Contents |
| --- | --- |
| `openmethane-cmaq-adjoint` | `/opt/cmaq/bin/ADJOINT_FWD` and `/opt/cmaq/bin/ADJOINT_BWD`. Depends on the libraries they link, which `dpkg-shlibdeps` generates, and `mpich` for `mpirun`. |

**Images are no longer published.** The last image is
`ghcr.io/openmethane/cmaq-adjoint:2.1.0`, and it and earlier ones stay
available. The `CMAQ_ADJOINT_VERSION` environment variable went with them:
`dpkg-query -W openmethane-cmaq-adjoint` gives the installed version.

### Installing

Pin the package by URL and checksum, for example in a Dockerfile. In the file
name of a pre-release, GitHub shows the `~` of the version as `.`, so
`2.2.0~rc1` is `openmethane-cmaq-adjoint_2.2.0.rc1_amd64.deb`.

```dockerfile
FROM scratch AS cmaq-adjoint-deb
ADD --checksum=sha256:<sum> https://github.com/openmethane/CMAQ-Adjoint/releases/download/v<ver>/openmethane-cmaq-adjoint_<ver>_amd64.deb /

FROM debian:trixie-slim
RUN --mount=type=bind,from=cmaq-adjoint-deb,target=/debs \
    apt-get update && apt-get install -y /debs/*.deb
```

Then run the adjoint with:

```shell
# run the fwd adjoint
/opt/cmaq/bin/ADJOINT_FWD

# run the bwd adjoint
/opt/cmaq/bin/ADJOINT_BWD
```

These binaries will need input data and a number of environment variables,
which should be familiar to users of the adjoint.

### Testing an unmerged change downstream

Define the package as a named stage, as above, and install it from a bind
mount on that stage. BuildKit lets a build replace a named stage with
`--build-context <stage>=<source>`, without editing the Dockerfile.

Locally, from the pull request's CI artifact `debs`:

```shell
# <branch> is this repository's pull request branch
RUN_ID=$(gh run list -R openmethane/CMAQ-Adjoint --branch <branch> --workflow build.yaml \
  --status success --limit 1 --json databaseId -q '.[0].databaseId')
gh run download "$RUN_ID" -R openmethane/CMAQ-Adjoint -n debs -D debs/cmaq-adjoint
docker build --build-context cmaq-adjoint-deb=./debs/cmaq-adjoint -t openmethane:dev .
```

In downstream CI, or in AWS, the package needs a public URL. Run the
`Publish dev build` workflow on the pull request's branch:

```shell
gh workflow run dev-build.yaml -R openmethane/CMAQ-Adjoint --ref <branch>
```

It publishes a pre-release tagged `dev-pr<PR>-<sha7>` with the package, and
writes the `ADD --checksum=…` line to the job summary. Pin that in the
downstream pull request temporarily, and move the pin to a real release before
it merges. The pre-release is deleted when this repository's pull request
closes.

This repository is also downstream of
[openmethane/CMAQ](https://github.com/openmethane/CMAQ), whose
`openmethane-cmaq-dev` package it builds against. To build against an unmerged
one, download that pull request's `debs` artifact and build with
`--build-context cmaq-dev=./debs/cmaq`. To run CI against it, pin a dev
pre-release of CMAQ in the `Dockerfile`'s `CMAQ_DEV_URL` and `CMAQ_DEV_SHA256`.

### Building

`make build` builds the package in Docker and exports it to `dist/`, and
`make test` installs it in a clean `debian:trixie-slim` and runs the tests. See
[docs/building.md](docs/building.md) for the steps, and for building without
Docker.

## Development

See [docs/development.md](docs/development.md) for debugging the build,
running the tests and how pull requests are reviewed.

# Citations

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3780216.svg)](https://doi.org/10.5281/zenodo.3780216)

This repo contains the work of many contributors. Attribution for specific
features and changes can be found in the change log present in the headers of
some individual files, and in the authors in git history.
