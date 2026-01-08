FROM continuumio/miniconda3 as conda

COPY environment.yml /opt/environment.yml
RUN conda env create -f /opt/environment.yml

# Install conda-pack:
RUN conda install -c conda-forge conda-pack

# Use conda-pack to create a standalone enviornment
# in /venv:
RUN conda-pack -n cmaq -o /tmp/env.tar && \
  mkdir /opt/venv && cd /opt/venv && tar xf /tmp/env.tar && \
  rm /tmp/env.tar

# We've put venv in same path it'll be in final image,
# so now fix up paths:
RUN /opt/venv/bin/conda-unpack

FROM debian:bookworm as build

ARG DEBIAN_FRONTEND=noninteractive
ENV CMAQ_VERSION="5.0.2"
ENV PATH=/opt/venv/bin:$PATH
ENV LD_LIBRARY_PATH=/opt/venv/bin:$LD_LIBRARY_PATH

RUN apt-get update && \
    apt-get install -y build-essential m4 csh wget && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /opt/cmaq

COPY --from=conda /opt/venv /opt/venv

# Build ioapi
COPY templates/ioapi /opt/cmaq/templates/ioapi
COPY scripts/common.sh /opt/cmaq/scripts/common.sh
COPY scripts/build_00_ioapi.sh /opt/cmaq/scripts/build_00_ioapi.sh
RUN bash /opt/cmaq/scripts/build_00_ioapi.sh

# Build a modified version of CMAQ in ch4 only mode
COPY templates/cmaq /opt/cmaq/templates/cmaq
COPY src/CMAQv5.0.2_notpollen /opt/cmaq/CMAQv5.0.2_notpollen
COPY scripts/build_10_cmaq.sh /opt/cmaq/scripts/build_10_cmaq.sh
RUN bash /opt/cmaq/scripts/build_10_cmaq.sh

# Build the CMAQ adjoint
COPY scripts/build_11_cmaq_adj.sh /opt/cmaq/scripts/build_11_cmaq_adj.sh
COPY src/cmaq_adj /opt/cmaq/cmaq_adj
RUN bash /opt/cmaq/scripts/build_11_cmaq_adj.sh


FROM debian:bookworm AS runtime

MAINTAINER Jared Lewis <jared.lewis@climate-resource.com>

ENV TZ=Etc/UTC
ENV CMAQ_VERSION="5.0.2"
ENV PATH=/opt/venv/bin:$PATH
ENV LD_LIBRARY_PATH=/opt/venv/bin:$LD_LIBRARY_PATH

WORKDIR /opt/cmaq
COPY --from=conda /opt/venv /opt/venv
COPY --from=build /opt/cmaq /opt/cmaq

RUN apt-get update && \
    apt-get install -y make csh wget && \
    rm -rf /var/lib/apt/lists/*

ENTRYPOINT ["/bin/bash"]
