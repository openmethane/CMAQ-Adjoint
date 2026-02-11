# CMAQ Adjoint

CMAQ forward and backward adjoint, based on
[CMAQ 5.0 Adjoint](https://zenodo.org/records/3780216),
forked from the
[development repo](https://adjoint.colorado.edu:8080) at University of
Colorado Boulder.

Additional changes have been made to support CH4 modeling for the Open Methane
project.

## Build and run

The easiest way to run the CMAQ Adjoint is using
[Docker](https://www.docker.com/get-started/).

To build the docker container, tagged `cmaq-adj`:

```shell
make build
```

Then run the adjoint in the container with:

```shell
# run the fwd adjoint
docker run -it --rm cmaq-adj /opt/cmaq/bin/ADJOINT_FWD

# run the bwd adjoint
docker run -it --rm cmaq-adj /opt/cmaq/bin/ADJOINT_BWD
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

The final docker container only includes the forward and backward adjoint
binaries in `/opt/cmaq/bin`. To debug the `builder` image:

```shell
# build the "builder" image from the Dockerfile
docker build . --progress=plain --target builder -t cmaq-adj-builder

# run an interactive bash terminal in the builder
docker run -it --rm cmaq-adj-builder bash
```

The `builder` will include the source code, libraries and build artifacts in
`/opt/cmaq`.

## Testing

A simple test script can be invoked with:

```shell
make test
```

This mounts the `tests` folder under `/opt/tests` and runs the script under
`/opt/tests/test-run.sh`. The test script downloads a bundle of test data and
a file describing necessary environment variables, and then runs the forward
adjoint on this input.

**Note:** this does **not** perform a numerical test, it simply ensures that
the compiled binary will successfully run to completion on known inputs.

# Citations

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3780216.svg)](https://doi.org/10.5281/zenodo.3780216)
