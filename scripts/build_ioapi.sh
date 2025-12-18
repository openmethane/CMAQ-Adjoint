#!/bin/bash
#
# Copyright 2024 The Superpower Institute
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
. scripts/common.sh
# Build ioapi
IOAPI_VERSION="3.1"

cd $ROOT
wget -nv https://www.cmascenter.org/ioapi/download/ioapi-${IOAPI_VERSION}.tar.gz  -O ioapi-${IOAPI_VERSION}.tar.gz
mkdir -p ioapi
tar -xzvf  ioapi-${IOAPI_VERSION}.tar.gz -C ioapi
pushd ioapi

# Copy in Makefile template
# This uses the correct configuration for the built NetCDF library
cp $ROOT/templates/ioapi/Makefile Makefile
cp $ROOT/templates/ioapi/Makeinclude.* ioapi/
# The following makefiles will be overridden by the .nocpl.sed templates
# but need to exist for the configure step to work
#touch ioapi/Makefile m3tools/Makefile

mkdir -p $BIN

make configure
make

popd

# Clean up
rm ioapi-${IOAPI_VERSION}.tar.gz
