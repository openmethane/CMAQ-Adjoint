#!/usr/bin/env python3
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
"""Build the ADJ_FORCE input used by the backward smoke test.

This is a development-time tool, not part of the test run. The file it
produces is shipped in the test-data bundle, and this script exists so that
file is reproducible rather than an opaque binary of unknown provenance.

It needs netCDF4, which the runtime image deliberately does not carry. Run it
on a host with py4dvar's virtualenv available, e.g.

    /path/to/openmethane/.venv/bin/python tests/make-force-file.py \\
        --template /path/to/openmethane/tests/test-data/templates/force_template.nc \\
        --output tests/test-data/FORCE.20221207.nc

The forcing is deliberately a single non-zero cell in the *final* record
rather than a broad field. The backward driver walks output steps from the run
end towards the run start, so the final record is the first one it reads. If
the forcing were spread over every record, the test could not tell a correct
run from one that reads the record one output step out of step -- it would
pick up a plausible non-zero field either way. With one populated record at
the run end, an off-by-one read misses it entirely and LGRID stays identically
zero, which the test asserts against.
"""

import argparse
import pathlib
import shutil
import sys

import netCDF4


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--template",
        required=True,
        type=pathlib.Path,
        help="py4dvar force_template.nc to base the file on",
    )
    parser.add_argument(
        "--output", required=True, type=pathlib.Path, help="file to write"
    )
    parser.add_argument(
        "--species", default="CH4", help="forcing variable name (default: CH4)"
    )
    parser.add_argument(
        "--record",
        type=int,
        default=-1,
        help="time record index to populate (default: -1, the last record)",
    )
    parser.add_argument("--layer", type=int, default=0, help="layer index (default: 0)")
    parser.add_argument("--row", type=int, default=5, help="row index (default: 5)")
    parser.add_argument("--col", type=int, default=5, help="column index (default: 5)")
    parser.add_argument(
        "--value", type=float, default=1.0, help="forcing value (default: 1.0)"
    )
    args = parser.parse_args(argv)

    if not args.template.is_file():
        parser.error(f"template not found: {args.template}")

    args.output.parent.mkdir(parents=True, exist_ok=True)
    shutil.copyfile(args.template, args.output)

    with netCDF4.Dataset(args.output, "a") as ds:
        if args.species not in ds.variables:
            parser.error(
                f"{args.species} not in {args.template} "
                f"(has: {', '.join(ds.variables)})"
            )
        var = ds.variables[args.species]
        var[:] = 0.0
        var[args.record, args.layer, args.row, args.col] = args.value

        tflag = ds.variables["TFLAG"][args.record, 0, :]
        date, time = int(tflag[0]), int(tflag[1])

    print(f"wrote {args.output}")
    print(
        f"  {args.species}[record {args.record}, layer {args.layer}, "
        f"row {args.row}, col {args.col}] = {args.value}"
    )
    print(f"  that record is stamped {date}:{time:06d}")
    print("  all other values are zero")
    return 0


if __name__ == "__main__":
    sys.exit(main())
