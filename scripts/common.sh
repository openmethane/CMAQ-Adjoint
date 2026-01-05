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
export CMAQ_REPO=/home/peter/work/openmethane-beta/cmaq-zenodo/
export CMAQ_DIRNAME=CCTM
export CONDA_INC=${HOME}/anaconda3/envs/cmaq/include
export CONDA_LIB=${HOME}/anaconda3/envs/cmaq/lib
export MPICH_DIR=/home/peter/anaconda3/lib/
export MPICH_INC=/home/peter/anaconda3/include
export IOAPI_DIR=${ROOT}/ioapi
export MCIP_DIR=/home/peter/work/openmethane-beta/CMAQ-EPA/PREP/mcip
