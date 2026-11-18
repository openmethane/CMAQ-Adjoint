# CMAQ Adjoint

CMAQ forward and backward adjoint, based on
[CMAQ 5.0 Adjoint](https://zenodo.org/records/3780216),
forked from the
[development repo](https://adjoint.colorado.edu:8080) at University of
Colorado Boulder.

Additional changes have been made to support CH4 modelling for the Open Methane
project.

## Build and run

The easiest way to run the CMAQ Adjoint is using
[Docker](https://www.docker.com/get-started/).

To build the docker container, tagged `cmaq-adjoint`:

```shell
make build
```

Then run the adjoint in the container with:

```shell
# run the fwd adjoint
docker run -it --rm cmaq-adjoint /opt/cmaq/bin/ADJOINT_FWD

# run the bwd adjoint
docker run -it --rm cmaq-adjoint /opt/cmaq/bin/ADJOINT_BWD
```

These binaries will need input data and a number of environment variables,
which should be familiar to users of the adjoint.

### Running without docker

This project uses [openmethane/CMAQ](https://github.com/openmethane/CMAQ) as
the base docker image, which provides a build environment pre-configured for
compiling and running CMAQ-based tools.

If you wish to build and run the forward and backward adjoint in a different
environment, that repo provides a good reference for the packages, dependencies
and configuration which will need to be provided.

## Development

See [docs/development.md](docs/development.md) for debugging the build,
running the tests and how pull requests are reviewed.

# Citations

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3780216.svg)](https://doi.org/10.5281/zenodo.3780216)

This repo contains the work of many contributors. Attribution for specific
features and changes can be found in the change log present in the headers of
some individual files, and in the authors in git history.
