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
"""Build the hourly BCON input used by the boundary sensitivity test.

This is a development-time tool, not part of the test run. The file it
produces ships in the test data bundle as
cmaq/BCON.d01.au-test_v1.CH4only.hourly.nc, and this script exists so that
file is reproducible rather than an opaque binary of unknown provenance.

It needs netCDF4, which the runtime image deliberately does not carry, e.g.

    uv run --no-project --with netCDF4 tests/make-hourly-bcon.py \\
        --template tests/test-data/cmaq/BCON.d01.au-test_v1.CH4only.nc \\
        --output tests/test-data/cmaq/BCON.d01.au-test_v1.CH4only.hourly.nc

The test bundle's BCON has a single time-independent record (TSTEP = 0), so
it never exercises the backward model's interpolation between BCON records.
Production BCON is hourly: one record per hour of the day plus one, the same
time axis as the meteorology. This writes that shape for the test day,
2022-12-07, 25 records from 00:00 to 24:00.

The values are the template's, scaled by a factor that differs by hour, by
perimeter edge and by layer. That keeps them physically sensible while
making a record, edge or layer that ADJ_LGRID_BC puts in the wrong place show
up as a different value, rather than as the same number in a different slot.
"""

import argparse
import math
import pathlib
import sys

import netCDF4
import numpy as np

SDATE = 2022341  # 2022-12-07, the test day
STIME = 0
TSTEP = 10000  # one hour, HHMMSS
NREC = 25


def edge_index(ncols, nrows):
    """Edge (0=S, 1=E, 2=N, 3=W) of each NTHIK=1 perimeter slot.

    I/O API orders the perimeter S, E, N, W, each side carrying the corner at
    its counter-clockwise end: S is NCOLS+1 long, E NROWS+1, N NCOLS+1 and
    W NROWS+1.
    """
    return np.concatenate(
        [
            np.full(ncols + 1, 0),
            np.full(nrows + 1, 1),
            np.full(ncols + 1, 2),
            np.full(nrows + 1, 3),
        ]
    )


def scale(nrec, nlays, edges):
    """Multiplicative factor, shape (nrec, nlays, nperim)."""
    hour = np.arange(nrec)[:, None, None]
    layer = np.arange(nlays)[None, :, None]
    edge = edges[None, None, :]
    return (
        1.0
        + 0.01 * np.sin(2.0 * math.pi * hour / 24.0)
        + 0.005 * (edge + 1)
        + 0.0002 * layer
    )


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--template",
        required=True,
        type=pathlib.Path,
        help="the test bundle's single-record BCON file",
    )
    parser.add_argument(
        "--output", required=True, type=pathlib.Path, help="file to write"
    )
    args = parser.parse_args(argv)

    if not args.template.is_file():
        parser.error(f"template not found: {args.template}")

    args.output.parent.mkdir(parents=True, exist_ok=True)

    with netCDF4.Dataset(args.template) as src, netCDF4.Dataset(
        args.output, "w", format="NETCDF3_CLASSIC"
    ) as dst:
        if src.FTYPE != 2 or src.NTHIK != 1:
            parser.error(f"{args.template} is not an NTHIK=1 boundary file")

        names = [v for v in src.variables if v != "TFLAG"]
        nlays = len(src.dimensions["LAY"])
        nperim = len(src.dimensions["PERIM"])
        edges = edge_index(src.NCOLS, src.NROWS)
        assert edges.size == nperim

        dst.createDimension("TSTEP", None)
        dst.createDimension("DATE-TIME", 2)
        dst.createDimension("LAY", nlays)
        dst.createDimension("VAR", len(names))
        dst.createDimension("PERIM", nperim)

        attrs = {k: src.getncattr(k) for k in src.ncattrs()}
        for k in ("NCO", "_NCProperties"):
            attrs.pop(k, None)
        # I/O API keeps the file's update description (UPDSC) in HISTORY and
        # DESC3 fails without it
        attrs["HISTORY"] = (
            f"tests/make-hourly-bcon.py from {args.template.name}\n"
        )
        attrs.update(SDATE=np.int32(SDATE), STIME=np.int32(STIME),
                     TSTEP=np.int32(TSTEP))
        attrs["FILEDESC"] = (
            "Hourly BCON for the CMAQ-Adjoint boundary sensitivity test, "
            "derived from the test bundle's single-record BCON by "
            "tests/make-hourly-bcon.py"
        ).ljust(80 * 60)
        dst.setncatts(attrs)

        tflag = dst.createVariable("TFLAG", "i4", ("TSTEP", "VAR", "DATE-TIME"))
        tflag.setncatts({k: src["TFLAG"].getncattr(k)
                         for k in src["TFLAG"].ncattrs()})
        for t in range(NREC):
            # every record of the test day is on 2022341 apart from 24:00,
            # which I/O API stamps as 00:00 the next day
            date, hour = (SDATE, t) if t < 24 else (SDATE + 1, 0)
            tflag[t, :, 0] = date
            tflag[t, :, 1] = hour * 10000

        factor = scale(NREC, nlays, edges)
        for name in names:
            var = dst.createVariable(name, "f4", ("TSTEP", "LAY", "PERIM"))
            var.setncatts({k: src[name].getncattr(k)
                           for k in src[name].ncattrs()})
            var[:] = (src[name][0] * factor).astype("f4")

    print(f"wrote {args.output}")
    print(f"  {NREC} records from {SDATE}:{STIME:06d} every {TSTEP:06d}")
    print(f"  variables: {', '.join(names)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
