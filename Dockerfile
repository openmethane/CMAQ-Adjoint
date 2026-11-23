# syntax=docker/dockerfile:1
#
# Builds ADJOINT_FWD and ADJOINT_BWD and packages them as the Debian package
# openmethane-cmaq-adjoint (see README.md). No image from this Dockerfile is
# published. Export the package with
#   docker build --target deb --output type=local,dest=dist .
#
# Each binary is one RUN, and runs the same scripts/build.sh command as a
# native build (docs/building.md).

# The openmethane-cmaq-dev package the build runs against. To update, change
# the version and checksum together: the checksum is in the SHA256SUMS of the
# release at https://github.com/openmethane/CMAQ/releases
ARG CMAQ_VERSION=1.1.0
ARG CMAQ_DEV_SHA256=23cf52d92e5f2ac96ac3f48b7fad7e00887307261360447298e2aa78daf6fdd8
ARG CMAQ_DEV_URL=https://github.com/openmethane/CMAQ/releases/download/v${CMAQ_VERSION}/openmethane-cmaq-dev_${CMAQ_VERSION}_amd64.deb

# PEP 440 version converted to a Debian version, see scripts/debian-version.sh
ARG DEB_VERSION=0.0.0~development

# ---------------------------------------------------------------------------
# openmethane-cmaq-dev, as a file to install from. Replace this stage with
#   --build-context cmaq-dev=<directory holding the .deb>
# to build against an unmerged version.
FROM scratch AS cmaq-dev
ARG CMAQ_DEV_URL
ARG CMAQ_DEV_SHA256
ADD --checksum=sha256:${CMAQ_DEV_SHA256} ${CMAQ_DEV_URL} /

# ---------------------------------------------------------------------------
# The build environment is only openmethane-cmaq-dev, and what it depends on
FROM debian:trixie-slim AS builder

ARG DEBIAN_FRONTEND=noninteractive

RUN --mount=type=bind,from=cmaq-dev,target=/debs <<EOT
apt-get update -qy
apt-get install -qyy -o APT::Install-Recommends=false /debs/openmethane-cmaq-dev_*.deb
apt-get clean
rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
EOT

# Where config.cmaq (from openmethane-cmaq-dev) puts the binaries
ENV BIN_DIR=/opt/cmaq/bin

WORKDIR /opt/cmaq-adjoint

# The sources are kept apart from anything that CMAQ installs
COPY cmaq /opt/cmaq-adjoint/cmaq
ENV ADJOINT_SRC=/opt/cmaq-adjoint/cmaq

COPY scripts/build.sh /opt/cmaq-adjoint/scripts/build.sh
COPY scripts/bldit.adjoint.fwd.openmethane /opt/cmaq-adjoint/scripts/bldit.adjoint.fwd.openmethane
RUN /opt/cmaq-adjoint/scripts/build.sh fwd

COPY scripts/bldit.adjoint.bwd.openmethane /opt/cmaq-adjoint/scripts/bldit.adjoint.bwd.openmethane
RUN /opt/cmaq-adjoint/scripts/build.sh bwd

# ---------------------------------------------------------------------------
# openmethane-cmaq-adjoint: the two binaries, where the build put them
FROM builder AS package

RUN <<EOT
apt-get update -qy
apt-get install -qyy -o APT::Install-Recommends=false dpkg-dev
EOT

ARG DEB_VERSION
COPY scripts/build-deb.sh /opt/cmaq-adjoint/scripts/build-deb.sh
RUN <<EOT
set -e
/opt/cmaq-adjoint/scripts/build-deb.sh \
    openmethane-cmaq-adjoint "$DEB_VERSION" \
    "Open Methane CMAQ adjoint: the forward and backward models for CH4" \
    /out \
    $BIN_DIR/ADJOINT_FWD \
    $BIN_DIR/ADJOINT_BWD \
    -- mpich
EOT

# ---------------------------------------------------------------------------
# Only the package, for export with --target deb
FROM scratch AS deb
COPY --from=package /out/*.deb /

# ---------------------------------------------------------------------------
# openmethane-cmaq-adjoint installed into a clean image, which must resolve
# every library. Run the tests in it with
#   docker run --rm -v ./tests:/opt/tests <image> /opt/tests/test-run.sh
FROM debian:trixie-slim AS test

ARG DEBIAN_FRONTEND=noninteractive

RUN --mount=type=bind,from=package,source=/out,target=/debs <<EOT
apt-get update -qy
apt-get install -qyy -o APT::Install-Recommends=false /debs/openmethane-cmaq-adjoint_*.deb
apt-get clean
rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
EOT

# Before wget is installed, so that nothing it brings in hides a missing dependency
COPY tests/check-libs.sh /opt/tests/check-libs.sh
RUN /opt/tests/check-libs.sh openmethane-cmaq-adjoint

# Only the tests need these: test-run.sh downloads the test data
RUN <<EOT
apt-get update -qy
apt-get install -qyy -o APT::Install-Recommends=false ca-certificates wget
apt-get clean
rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
EOT

CMD ["/opt/tests/test-run.sh"]
