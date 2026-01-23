#!/usr/bin/env bash
#
# Copyright 2026 The Superpower Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

set -e

ROOT="${PWD:-/opt/cmaq}"
CMAQ_DIRNAME="${CMAQ_DIRNAME:-CCTM}"
IOAPI_DIR="${IOAPI_DIR:-$ROOT/ioapi}"

export CONDA_INC="/opt/venv/include"
export CONDA_LIB="/opt/venv/lib"
export MPICH_DIR=${CONDA_LIB}
export MPICH_INC=${CONDA_INC}

# build the bldmake tool
pushd BLDMAKE_git
make
popd

# use bldmake to build forward and backward adjoint
pushd ${ROOT}
csh scripts/bldit.adjoint.fwd.openmethane
csh scripts/bldit.adjoint.bwd.openmethane
popd
