# Install required dependencies from conda
FROM continuumio/miniconda3:25.3.1-1 as conda

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

# Define environment variables needed for CMAQ and lib Makefiles
ENV ROOT=/opt/cmaq
env CMAQ_DIRNAME=CCTM
ENV IOAPI_DIR=$ROOT/ioapi

WORKDIR /opt/cmaq

# Build ioapi
COPY lib/ioapi /opt/ioapi
RUN /opt/ioapi/build_ioapi.sh

# Build pario and stenex libraries
COPY lib/pario /opt/pario
COPY lib/stenex /opt/stenex
COPY scripts/build_libs.sh /opt/scripts/build_libs.sh
RUN /opt/scripts/build_libs.sh

# Build the bldmake tool, needed by CMAQ build scripts
COPY tools/BLDMAKE_git /opt/BLDMAKE_git
RUN <<EOT
cd /opt/BLDMAKE_git
make
mv bldmake /usr/local/bin # ensure bldmake is in the PATH
EOT

# Build a modified version of CMAQ-Adjoint in ch4 only mode
COPY cmaq /opt/cmaq
COPY scripts/build_adj.sh /opt/cmaq/scripts/build_adj.sh
COPY scripts/bldit.adjoint.fwd.openmethane /opt/cmaq/scripts/bldit.adjoint.fwd.openmethane
COPY scripts/bldit.adjoint.bwd.openmethane /opt/cmaq/scripts/bldit.adjoint.bwd.openmethane
RUN /opt/cmaq/scripts/build_adj.sh

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

COPY --from=builder /opt/cmaq/BLD_fwd_CH4only/ADJOINT_FWD /opt/cmaq/bin/ADJOINT_FWD
COPY --from=builder /opt/cmaq/BLD_bwd_CH4only/ADJOINT_BWD /opt/cmaq/bin/ADJOINT_BWD

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
