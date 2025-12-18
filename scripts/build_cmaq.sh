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

set -e
pushd ${ROOT}
# make some libraries
for item in pario stenex/se; do
  pushd $item
  echo "Building $item"
  cp ${ROOT}/templates/${item}/makefile.gcc .
  make -f makefile.gcc
  popd
done

pushd BLDMAKE_git
make
popd


