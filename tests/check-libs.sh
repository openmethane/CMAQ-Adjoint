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
# Fail if any ELF file in the installed openmethane-cmaq package has a
# shared library that cannot be found.

set -euo pipefail

PACKAGE="${1:-openmethane-cmaq}"

echo "Depends of ${PACKAGE}: $(dpkg-query -W -f='${Depends}' "${PACKAGE}")"

status=0
count=0
while read -r file; do
  [ -f "$file" ] || continue
  head -c4 "$file" | grep -q '^.ELF' || continue
  count=$((count + 1))
  if ldd "$file" | grep 'not found'; then
    echo "FAILED: $file has missing shared libraries"
    status=1
  fi
done < <(dpkg -L "${PACKAGE}")

if [ "$count" -eq 0 ]; then
  echo "FAILED: no ELF files found in ${PACKAGE}"
  exit 1
fi
[ "$status" -eq 0 ] && echo "PASSED: all ${count} ELF files in ${PACKAGE} resolve their libraries"
exit "$status"
