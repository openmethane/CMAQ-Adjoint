#!/bin/bash
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

# Build ioapi
IOAPI_DIR="${IOAPI_DIR:-/opt/cmaq/ioapi}"
IOAPI_VERSION="3.1"
IOAPI_TAR="ioapi-${IOAPI_VERSION}.tar.gz"

# path to this script and its parent folder
LIB_IOAPI_DIR="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

# fetch ioapi source and extract it
wget -nv https://www.cmascenter.org/ioapi/download/ioapi-${IOAPI_VERSION}.tar.gz  -O ${IOAPI_TAR}
mkdir -p "${IOAPI_DIR}"
tar -xzvf ${IOAPI_TAR} -C "${IOAPI_DIR}"
rm ${IOAPI_TAR}

cd "${IOAPI_DIR}"

# Copy in Makefile template
# This uses the correct configuration for the built NetCDF library
cp "${LIB_IOAPI_DIR}"/Makefile Makefile
cp "${LIB_IOAPI_DIR}"/Makeinclude.* ioapi/

# Set the binary architecture in $BIN, this will be used by the ioapi Makefile
ARCH=$(uname -m)
BLD_OS=$(/bin/uname -s)$(/bin/uname -r | cut -d. -f1)
export BIN=${BLD_OS}_${ARCH}gfort
mkdir -p $BIN

# Build ioapi
make configure
make
