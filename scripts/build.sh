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
# Build one of the adjoint binaries: build.sh <fwd|bwd|all>
#
# This only dispatches to the bldit script for each binary, so that the
# Dockerfile, CI and a native build all run the same commands. The build needs
# the libraries and config.cmaq from openmethane-cmaq-dev, see docs/building.md.

set -euo pipefail

SCRIPTS_DIR="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"

run_step() {
  case "$1" in
    fwd) cd "$SCRIPTS_DIR" && ./bldit.adjoint.fwd.openmethane ;;
    bwd) cd "$SCRIPTS_DIR" && ./bldit.adjoint.bwd.openmethane ;;
    *) echo "Unknown binary: $1 (expected fwd, bwd or all)" >&2; exit 2 ;;
  esac
}

if [ "$#" -ne 1 ]; then
  echo "Usage: $0 <fwd|bwd|all>" >&2
  exit 2
fi

if [ "$1" = "all" ]; then
  for step in fwd bwd; do
    echo "=== build.sh $step ==="
    ( run_step "$step" )
  done
else
  run_step "$1"
fi
