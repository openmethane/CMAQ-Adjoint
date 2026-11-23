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
# Package files that a build produced, where it produced them, as a .deb.
# Nothing is moved, renamed or rewritten.
#
# Usage: build-deb.sh [--no-shlibdeps] <name> <version> <description> \
#            <outdir> <file-or-dir>... [-- <extra Depends>...]
#
# Runtime dependencies are generated from the ELF files with dpkg-shlibdeps,
# unless --no-shlibdeps is given (for packages whose Depends: is a toolchain).

set -euo pipefail

shlibdeps=1
if [ "${1:-}" = "--no-shlibdeps" ]; then shlibdeps=0; shift; fi
if [ "$#" -lt 5 ]; then
  sed -n '/^# Usage/,/^$/p' "$0" >&2
  exit 2
fi
name="$1" version="$2" description="$3" outdir="$4"
shift 4

files=()
while [ "$#" -gt 0 ] && [ "$1" != "--" ]; do files+=("$1"); shift; done
[ "${1:-}" = "--" ] && shift
extra_deps=("$@")

arch="$(dpkg --print-architecture)"
work="$(mktemp -d)"
trap 'rm -rf "$work"' EXIT
pkg="$work/pkg"
mkdir -p "$pkg/DEBIAN" "$outdir"

# Copy each path to the same absolute path under the package root
for f in "${files[@]}"; do
  [ -e "$f" ] || { echo "build-deb.sh: $f does not exist" >&2; exit 1; }
  mkdir -p "$pkg$(dirname "$f")"
  cp -a "$f" "$pkg$f"
done

deps=()
if [ "$shlibdeps" = 1 ]; then
  mapfile -t elfs < <(find "$pkg" -type f -exec sh -c \
    'head -c4 "$1" | grep -q "^.ELF" && echo "$1"' _ {} \;)
  if [ "${#elfs[@]}" -gt 0 ]; then
    # dpkg-shlibdeps must run in a directory with a debian/control
    mkdir -p "$work/src/debian"
    printf 'Source: %s\n\nPackage: %s\nArchitecture: any\n' "$name" "$name" > "$work/src/debian/control"
    shlibs="$(cd "$work/src" && dpkg-shlibdeps -O "${elfs[@]}" | sed -n 's/^shlibs:Depends=//p')"
    [ -z "$shlibs" ] || deps+=("$shlibs")
  fi
fi
deps+=("${extra_deps[@]}")
depends="$(IFS=,; echo "${deps[*]}" | sed -E 's/ *, */, /g')"

cat > "$pkg/DEBIAN/control" <<CONTROL
Package: $name
Version: $version
Architecture: $arch
Maintainer: Open Methane <noreply@openmethane.invalid>
Section: science
Priority: optional
Depends: $depends
Description: $description
CONTROL

deb="$outdir/${name}_${version}_${arch}.deb"
dpkg-deb --root-owner-group --build "$pkg" "$deb" >/dev/null
echo "Built $deb"
echo "Depends of $name: $depends"
