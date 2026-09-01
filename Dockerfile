FROM ghcr.io/openmethane/cmaq:1.0.2 AS builder

ARG DEBIAN_FRONTEND=noninteractive

RUN <<EOT
apt-get update -qy
apt-get install -qyy \
    -o APT::Install-Recommends=false \
    -o APT::Install-Suggests=false \
    build-essential \
    m4 \
    csh \
    gfortran \
    libnetcdff-dev \
    libmpich-dev

apt-get clean
rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
EOT

# Define environment variables needed for CMAQ and lib Makefiles
ENV CMAQ_DIRNAME=CCTM
ENV BIN_DIR=/opt/cmaq/bin

WORKDIR /opt/cmaq

# Build a modified version of CMAQ-Adjoint in ch4 only mode
COPY cmaq /opt/cmaq/models
COPY scripts/bldit.adjoint.fwd.openmethane /opt/scripts/bldit.adjoint.fwd.openmethane
ARG FORCE_BUILD=unknown
RUN cd /opt/scripts && ./bldit.adjoint.fwd.openmethane
COPY scripts/bldit.adjoint.bwd.openmethane /opt/scripts/bldit.adjoint.bwd.openmethane
RUN cd /opt/scripts && ./bldit.adjoint.bwd.openmethane

# Then, use a final image without extra packages for our runtime environment
FROM ghcr.io/openmethane/cmaq:1.0.2

# These will be overwritten in GHA due to https://github.com/docker/metadata-action/issues/295
# These must be duplicated in .github/workflows/build_docker.yaml
LABEL org.opencontainers.image.title="CMAQ Adjoint"
LABEL org.opencontainers.image.description="CMAQ forward and backward adjoint"
LABEL org.opencontainers.image.authors="Peter Rayner <peter.rayner@superpowerinstitute.com.au>, Jared Lewis <jared.lewis@climate-resource.com>"
LABEL org.opencontainers.image.vendor="The Superpower Institute"

RUN <<EOT
apt-get update -qy
apt-get install -qyy \
    -o APT::Install-Recommends=false \
    -o APT::Install-Suggests=false \
    ca-certificates \
    libnetcdff7 \
    mpich \
    wget

apt-get clean
rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
EOT

# CMAQ_ADJOINT_VERSION will be overridden in release builds with semver vX.Y.Z
ARG CMAQ_ADJOINT_VERSION=development
# Make the $CMAQ_ADJOINT_VERSION available as an env var inside the container
ENV CMAQ_ADJOINT_VERSION=$CMAQ_ADJOINT_VERSION

COPY --from=builder /opt/cmaq/bin /opt/cmaq/bin

WORKDIR /opt/cmaq

ENTRYPOINT ["/bin/bash"]
