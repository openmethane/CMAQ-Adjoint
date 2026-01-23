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

LIB_ROOT="${LIB_ROOT:-/opt}"

cd "${LIB_ROOT}"

# build required libraries
for LIB in pario stenex/se; do
   pushd $LIB
   echo "Building $LIB"
   make -f makefile.gcc
   popd
done
