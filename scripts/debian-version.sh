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
# Convert a PEP 440 version to a Debian version, so that pre-releases sort
# before the release:
#   1.1.0 -> 1.1.0   1.1.1.dev1 -> 1.1.1~dev1   1.2.0rc1 -> 1.2.0~rc1
# With no argument, reads the version from pyproject.toml.

set -euo pipefail

if [ "$#" -gt 1 ]; then
  echo "Usage: $0 [pep440-version]" >&2
  exit 2
fi

if [ "$#" -eq 1 ]; then
  version="$1"
else
  root="$( cd -- "$(dirname "$0")/.." >/dev/null 2>&1 ; pwd -P )"
  version="$(sed -n 's/^version = "\(.*\)"$/\1/p' "$root/pyproject.toml")"
fi

sed -E 's/\.?(dev|a|b|rc)([0-9]+)/~\1\2/g' <<<"${version#v}"
