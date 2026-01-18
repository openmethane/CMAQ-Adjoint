# Install required dependencies from conda
FROM continuumio/miniconda3 as conda

# Install and package up the conda environment
# Creates a standalone environment in /opt/venv
COPY environment.yml /opt/environment.yml
RUN conda env create -f /opt/environment.yml
RUN conda install -c conda-forge conda-pack
RUN conda-pack -n cmaq_adj -o /tmp/env.tar && \
  mkdir /opt/venv && cd /opt/venv && \
  tar xf /tmp/env.tar && \
  rm /tmp/env.tar

# We've put venv in same path it'll be in final image,
# so now fix up paths:
RUN /opt/venv/bin/conda-unpack

# Build dependencies and adjoint from source
FROM debian:bookworm-slim as builder

ARG DEBIAN_FRONTEND=noninteractive
ENV PATH=/opt/venv/bin:$PATH
ENV LD_LIBRARY_PATH=/opt/venv/bin:$LD_LIBRARY_PATH

RUN <<EOT
apt-get update -qy
apt-get install -qyy \
    -o APT::Install-Recommends=false \
    -o APT::Install-Suggests=false \
    ca-certificates \
    build-essential \
    m4 \
    csh \
    wget

apt-get clean
rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
EOT

COPY --from=conda /opt/venv /opt/venv

WORKDIR /opt/cmaq

# Build ioapi
COPY templates/ioapi /opt/cmaq/templates/ioapi
COPY scripts/common.sh /opt/cmaq/scripts/common.sh
COPY scripts/build_ioapi.sh /opt/cmaq/scripts/build_ioapi.sh
RUN bash /opt/cmaq/scripts/build_ioapi.sh

# Build a modified version of CMAQ in ch4 only mode
COPY BLDMAKE_git /opt/cmaq/BLDMAKE_git
COPY CCTM /opt/cmaq/CCTM
COPY ICL /opt/cmaq/ICL
COPY pario /opt/cmaq/pario
COPY stenex /opt/cmaq/stenex
COPY scripts/build_all.sh /opt/cmaq/scripts/build_all.sh
COPY scripts/bldit.adjoint.fwd.openmethane /opt/cmaq/scripts/bldit.adjoint.fwd.openmethane
COPY scripts/bldit.adjoint.bwd.openmethane /opt/cmaq/scripts/bldit.adjoint.bwd.openmethane
RUN bash /opt/cmaq/scripts/build_all.sh

# Then, use a final image without extra packages for our runtime environment
FROM debian:bookworm-slim

# These will be overwritten in GHA due to https://github.com/docker/metadata-action/issues/295
# These must be duplicated in .github/workflows/build_docker.yaml
LABEL org.opencontainers.image.title="CMAQ Adjoint"
LABEL org.opencontainers.image.description="CMAQ forward and backward adjoint"
LABEL org.opencontainers.image.authors="Peter Rayner <peter.rayner@superpowerinstitute.com.au>, Jared Lewis <jared.lewis@climate-resource.com>"
LABEL org.opencontainers.image.vendor="The Superpower Institute"

# CMAQ_ADJ_VERSION will be overridden in release builds with semver vX.Y.Z
ARG CMAQ_ADJ_VERSION=development
# Make the $CMAQ_ADJ_VERSION available as an env var inside the container
ENV CMAQ_ADJ_VERSION=$CMAQ_ADJ_VERSION

ENV TZ=Etc/UTC
ENV CMAQ_VERSION="5.0.2"
ENV PATH=/opt/venv/bin:$PATH
ENV LD_LIBRARY_PATH=/opt/venv/bin:$LD_LIBRARY_PATH

COPY --from=conda /opt/venv /opt/venv

COPY --from=builder /opt/cmaq /opt/cmaq

RUN <<EOT
apt-get update -qy
apt-get install -qyy \
    -o APT::Install-Recommends=false \
    -o APT::Install-Suggests=false \
    ca-certificates \
    wget

apt-get clean
rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*
EOT

WORKDIR /opt/cmaq

ENTRYPOINT ["/bin/bash"]
