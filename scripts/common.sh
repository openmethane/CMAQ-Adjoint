#!/usr/bin/env bash

# Setup the build environment for iotools
ARCH=$(uname -m)
BLD_OS=$(/bin/uname -s)$(/bin/uname -r | cut -d. -f1)
export ARCH
export BLD_OS
export BIN=${BLD_OS}_${ARCH}gfort
export BASEDIR=$PWD
export CPLMODE=nocpl
export LD_LIBRARY_PATH="/opt/venv/lib:${LD_LIBRARY_PATH}"
export PATH="/opt/venv/bin:${PATH}"
export ROOT=$PWD
export CMAQ_DIRNAME=CCTM
export CONDA_INC="/opt/venv/include"
export CONDA_LIB="/opt/venv/lib"
export MPICH_DIR=${CONDA_LIB}
export MPICH_INC=${CONDA_INC}
export IOAPI_DIR=${ROOT}/ioapi
