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

# fail on any command exiting with non-zero status
set -e

# paths to the adjoint binaries being tested
ADJ_FWD_BIN="${ADJ_FWD_BIN:-/opt/cmaq/bin/ADJOINT_FWD}"
ADJ_BWD_BIN="${ADJ_BWD_BIN:-/opt/cmaq/bin/ADJOINT_BWD}"

# url to the test data bundle
TEST_DATA_URL="${TEST_DATA_URL:-https://openmethane.s3.amazonaws.com/tests/cmaq-adjoint/cmaq-adjoint-test-data.20261028.tar.gz}"

# path to this script and its parent folder
TESTS_ROOT="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
TEST_RUNS_ROOT="/tmp/cmaq-adjoint-tests"
mkdir -p "${TEST_RUNS_ROOT}"

# fetch test data too large to fit in the repo
TEST_DATA_DIR="${TESTS_ROOT}/test-data"
if [ ! -d "${TEST_DATA_DIR}" ]; then
  echo "Fetching test data"
  wget -c -nv "${TEST_DATA_URL}" \
    -O "${TEST_RUNS_ROOT}/test-data.tar.gz"
  mkdir -p "${TEST_DATA_DIR}"
  tar -xf "${TEST_RUNS_ROOT}/test-data.tar.gz" -C "${TEST_DATA_DIR}"
  rm "${TEST_RUNS_ROOT}/test-data.tar.gz"
fi
DATA_DIR="${TEST_DATA_DIR}"

# hourly BCON for the test day (see tests/make-hourly-bcon.py). The bundle's
# own BCON has a single time-independent record.
HOURLY_BCON="${TEST_DATA_DIR}/cmaq/BCON.d01.au-test_v1.CH4only.hourly.nc"

# The bundle gained fwd.env (renamed from test.env), bwd.env and the adjoint
# forcing file when the backward test was added, and the hourly BCON when the
# boundary sensitivity test was. Check for them up front so a stale cached
# bundle gives a clear message rather than a confusing model failure part way
# through.
for REQUIRED in fwd.env bwd.env FORCE.20221207.nc \
    cmaq/BCON.d01.au-test_v1.CH4only.hourly.nc; do
  if [ ! -f "${DATA_DIR}/${REQUIRED}" ]; then
    echo "ERROR: ${DATA_DIR}/${REQUIRED} is missing."
    echo "The cached test data bundle is out of date. Delete"
    echo "${TEST_DATA_DIR} to re-fetch it, or see tests/README.md to rebuild."
    exit 1
  fi
done

FAILURES=0

fail () {
  echo "FAIL: $1"
  FAILURES=$((FAILURES + 1))
}

prepare_test_run () {
  TEST_RUN_DIR=$1
  if [ -d "${TEST_RUN_DIR}" ]; then
    echo "Clearing previous run"
    rm -rf "${TEST_RUN_DIR}"
  fi
  mkdir -p "${TEST_RUN_DIR}/output"
  mkdir -p "${TEST_RUN_DIR}/chkpnt"

  RUN_DIR="${TEST_RUN_DIR}"
  DATA_DIR="${TEST_DATA_DIR}"
}

read_environment_vars () {
  ENV_FILE=$1
  # set environment variables which will control the adjoint run
  set -o allexport && source "${ENV_FILE}" && set +o allexport
}

# assert the run reached the end of its time loop rather than stopping early
assert_normal_completion () {
  LOG=$1
  PROGRAM=$2
  if ! grep -aq "Normal Completion of program ${PROGRAM}" "${LOG}"; then
    fail "${PROGRAM} did not report normal completion in ${LOG}"
  fi
}

# a run can complete "normally" while producing garbage, so check the log for
# non-finite values -- these show up in the model's own diagnostic output
assert_no_nan () {
  LOG=$1
  if grep -aqiE "\bnan\b|\-?infinity|\*{6,}" "${LOG}"; then
    fail "non-finite values in ${LOG}"
    grep -aniE "\bnan\b|\-?infinity|\*{6,}" "${LOG}" | head -5
  fi
}

# IO/API classic netCDF stores its data section as raw big-endian floats, so a
# file whose variables are entirely zero is almost all zero bytes -- only the
# header carries content. Counting non-zero bytes therefore distinguishes "the
# model wrote real numbers" from "the model wrote a correctly shaped field of
# zeros" without needing a netCDF-aware tool in the runtime image.
#
# The threshold is deliberately loose. Measured on this domain: a run with
# all-zero forcing produces ~6 kB of non-zero bytes (header only), a run with
# the real forcing ~245 kB. 20 kB sits well clear of both. This is a smoke
# test, not a numerical check -- it distinguishes "sensitivities were
# computed" from "the field is empty", nothing finer.
#
# Worth noting that the zero-forcing run still reports normal completion, so
# this check, not assert_normal_completion, is what catches a backward run
# that silently fails to pick the forcing up.
assert_has_nonzero_data () {
  FILE=$1
  MIN_BYTES=${2:-20000}
  if [ ! -f "${FILE}" ]; then
    fail "${FILE} was not created"
    return
  fi
  NONZERO=$(tr -d '\0' < "${FILE}" | wc -c)
  echo "  $(basename "${FILE}"): ${NONZERO} non-zero bytes"
  if [ "${NONZERO}" -lt "${MIN_BYTES}" ]; then
    fail "$(basename "${FILE}") looks empty (${NONZERO} non-zero bytes," \
         "expected at least ${MIN_BYTES})"
  fi
}

# IO/API boundary files are small, so the 20 kB default above would reject
# a correct one. Measured on this domain: the ADJ_LGRID_BC header alone is
# ~5.8 kB of non-zero bytes, and the bundle's single-record BCON adds ~3.3 kB
# of data on top (only the slots on inflow faces are non-zero), ~9.2 kB in
# all. 7.5 kB sits between the two. The hourly BCON gives ~72 kB.
ADJ_LGRID_BC_MIN_BYTES=7500

# Print the byte offset at which a classic (CDF1) or 64-bit-offset (CDF2)
# netCDF file's data section starts, i.e. the smallest "begin" of any
# variable, by walking the header. The runtime image has no netCDF tools,
# but perl-base is always there. An IO/API file holds only record variables,
# so everything from this offset to the end of the file is data.
nc_data_offset () {
  perl -e '
    open(my $f, "<:raw", $ARGV[0]) or die "$ARGV[0]: $!\n";
    local $/; my $h = <$f>; my $p = 4;
    my $offsz = (substr($h, 3, 1) eq "\x02") ? 8 : 4;
    sub i32 { my $v = unpack("N", substr($h, $p, 4)); $p += 4; $v }
    sub off { my $v = unpack($offsz == 8 ? "Q>" : "N", substr($h, $p, $offsz));
              $p += $offsz; $v }
    sub name { my $n = i32(); my $s = substr($h, $p, $n);
               $p += ($n + 3) & ~3; $s }
    # attribute list: tag, count, then name, type, nelems, padded values
    sub atts { i32(); for (1 .. i32()) { name(); my $t = i32(); my $n = i32();
                 $p += ($n * (0, 1, 1, 2, 4, 4, 8)[$t] + 3) & ~3 } }
    i32();                                           # numrecs
    i32(); for (1 .. i32()) { name(); i32() }        # dimensions
    atts();                                          # global attributes
    my $begin;
    i32(); for (1 .. i32()) {                        # variables
      name(); $p += 4 * i32(); atts(); i32(); i32();
      my $b = off(); $begin = $b if !defined $begin || $b < $begin }
    print $begin;
  ' "$1"
}

# Compare the data sections of two ADJ_LGRID_BC files. The backward model is
# not bit-reproducible across decompositions -- ADJ_LGRID itself differs
# between one and two processes at the ~1e-7 level -- so allow a difference
# of up to 1e-5 of the largest value. A perimeter slot gathered to the wrong
# place would be wrong by O(1).
#
# The data section is read from each file's own header, so the check follows
# the test bundle's grid, layers and variables. It is NREC records of TFLAG
# plus the variables, as big-endian floats. TFLAG is read as floats too,
# which is harmless since it is identical in both. od prints a non-finite
# float as nan or inf, which awk would read as 0, so those fail explicitly.
assert_bc_close () {
  FILE_A=$1
  FILE_B=$2
  if [ ! -f "${FILE_A}" ] || [ ! -f "${FILE_B}" ]; then
    fail "cannot compare ${FILE_A} and ${FILE_B}: missing file"
    return
  fi
  BYTES_A=$(( $(stat -c %s "${FILE_A}") - $(nc_data_offset "${FILE_A}") ))
  BYTES_B=$(( $(stat -c %s "${FILE_B}") - $(nc_data_offset "${FILE_B}") ))
  if [ "${BYTES_A}" -ne "${BYTES_B}" ] || [ "${BYTES_A}" -le 0 ]; then
    fail "$(basename "${FILE_B}") data sections differ in size" \
         "(${BYTES_A} and ${BYTES_B} bytes)"
    return
  fi
  echo "  comparing ${BYTES_A} bytes of data"
  if ! paste \
      <(tail -c "${BYTES_A}" "${FILE_A}" | od -An -v -t f4 --endian=big -w4) \
      <(tail -c "${BYTES_B}" "${FILE_B}" | od -An -v -t f4 --endian=big -w4) |
    awk 'tolower($0) ~ /nan|inf/ { bad++; next }
         { d = $1 - $2; if (d < 0) d = -d; a = ($1 < 0) ? -$1 : $1
           if (d > maxd) maxd = d; if (a > maxa) maxa = a }
         END { printf "  max |value| %g, max |difference| %g, non-finite %d\n",
                      maxa, maxd, bad
               exit !(bad == 0 && maxa > 0 && maxd <= 1e-5 * maxa) }'; then
    fail "$(basename "${FILE_B}") differs from the single-process run" \
         "or holds non-finite values"
  fi
}

use_hourly_bcon () {
  export BNDY_GASC_1="${HOURLY_BCON}"
  export BNDY_AERO_1="${HOURLY_BCON}"
  export BNDY_NONR_1="${HOURLY_BCON}"
  export BNDY_TRAC_1="${HOURLY_BCON}"
}

run_fwd () {
  pushd "${RUN_DIR}" > /dev/null
  ${MPI_PREFIX} ${ADJ_FWD_BIN}
  popd > /dev/null
}

run_bwd () {
  pushd "${RUN_DIR}" > /dev/null
  ${MPI_PREFIX} ${ADJ_BWD_BIN}
  popd > /dev/null
}

# --------------------------------------------------------------------------
echo "=== 001: forward, single process ==="
prepare_test_run "${TEST_RUNS_ROOT}/001-fwd-sp"
read_environment_vars "${TEST_DATA_DIR}/fwd.env"
MPI_PREFIX=""
run_fwd
cat "${LOGFILE}"
assert_normal_completion "${LOGFILE}" "DRIVER_FWD"

# --------------------------------------------------------------------------
echo "=== 002: forward, two processes ==="
prepare_test_run "${TEST_RUNS_ROOT}/002-fwd-mp"
read_environment_vars "${TEST_DATA_DIR}/fwd.env"
export NPCOL_NPROW="2 1"
MPI_PREFIX="mpirun -np 2"
run_fwd
cat "${LOGFILE}"
assert_normal_completion "${LOGFILE}" "DRIVER_FWD"

# --------------------------------------------------------------------------
# The backward run consumes the forward run's checkpoints, XFIRST state and
# concentration output, so it has to follow a forward run in the same
# directory rather than starting from a clean one.
echo "=== 003: forward then backward, single process ==="
prepare_test_run "${TEST_RUNS_ROOT}/003-bwd-sp"
read_environment_vars "${TEST_DATA_DIR}/fwd.env"
MPI_PREFIX=""
run_fwd
FWD_LOGFILE="${LOGFILE}"
assert_normal_completion "${FWD_LOGFILE}" "DRIVER_FWD"

read_environment_vars "${TEST_DATA_DIR}/bwd.env"
# not in the bundle's bwd.env, which predates it
export ADJ_LGRID_BC="${RUN_DIR}/output/ADJ_LGRID_BC.20221207.nc -v"
run_bwd
cat "${LOGFILE}"
assert_normal_completion "${LOGFILE}" "DRIVER_BWD"
assert_no_nan "${LOGFILE}"

# ADJ_FORCE holds a single non-zero cell in its final record, which is the
# first record the backward run reads. If the sensitivities come back all
# zero, the forcing never reached LGRID -- most likely because the record was
# read at the wrong time.
assert_has_nonzero_data "${RUN_DIR}/output/ADJ_LGRID.20221207.nc"
assert_has_nonzero_data "${RUN_DIR}/output/ADJ_LGRID_EM.20221207.nc"
assert_has_nonzero_data "${RUN_DIR}/output/ADJ_LGRID_BC.20221207.nc" \
  "${ADJ_LGRID_BC_MIN_BYTES}"

# --------------------------------------------------------------------------
# The same as 003 with an hourly BCON, which takes the backward model's
# boundary sensitivity through the interpolation between BCON records that
# the bundle's single-record BCON never exercises. The forward run must use
# the same BCON, since the backward run replays it.
echo "=== 004: forward then backward, hourly BCON, single process ==="
prepare_test_run "${TEST_RUNS_ROOT}/004-bwd-bc-sp"
read_environment_vars "${TEST_DATA_DIR}/fwd.env"
use_hourly_bcon
MPI_PREFIX=""
run_fwd
assert_normal_completion "${LOGFILE}" "DRIVER_FWD"

read_environment_vars "${TEST_DATA_DIR}/bwd.env"
export ADJ_LGRID_BC="${RUN_DIR}/output/ADJ_LGRID_BC.20221207.nc -v"
run_bwd
cat "${LOGFILE}"
assert_normal_completion "${LOGFILE}" "DRIVER_BWD"
assert_no_nan "${LOGFILE}"
assert_has_nonzero_data "${RUN_DIR}/output/ADJ_LGRID_BC.20221207.nc" \
  "${ADJ_LGRID_BC_MIN_BYTES}"
BC_SP="${RUN_DIR}/output/ADJ_LGRID_BC.20221207.nc"

# --------------------------------------------------------------------------
# 004 again on two processes. Each PE accumulates the part of the perimeter
# on its own edge of the domain and PE 0 gathers and writes them, so the
# output should match 004 to round-off. A column split puts the west and
# east faces on different PEs and divides the south and north faces between
# them; a row split does the opposite, so between them every face's local
# to global slot mapping is checked.
run_bc_two_processes () {
  NAME=$1
  prepare_test_run "${TEST_RUNS_ROOT}/${NAME}"
  read_environment_vars "${TEST_DATA_DIR}/fwd.env"
  use_hourly_bcon
  # after fwd.env, which sets NPCOL_NPROW="1 1"
  export NPCOL_NPROW=$2
  MPI_PREFIX="mpirun -np 2"
  run_fwd
  assert_normal_completion "${LOGFILE}" "DRIVER_FWD"

  read_environment_vars "${TEST_DATA_DIR}/bwd.env"
  export ADJ_LGRID_BC="${RUN_DIR}/output/ADJ_LGRID_BC.20221207.nc -v"
  run_bwd
  cat "${LOGFILE}"
  assert_normal_completion "${LOGFILE}" "DRIVER_BWD"
  assert_no_nan "${LOGFILE}"
  assert_has_nonzero_data "${RUN_DIR}/output/ADJ_LGRID_BC.20221207.nc" \
    "${ADJ_LGRID_BC_MIN_BYTES}"
  assert_bc_close "${BC_SP}" "${RUN_DIR}/output/ADJ_LGRID_BC.20221207.nc"
}

echo "=== 005: forward then backward, hourly BCON, two processes by column ==="
run_bc_two_processes 005-bwd-bc-mp-col "2 1"

echo "=== 006: forward then backward, hourly BCON, two processes by row ==="
run_bc_two_processes 006-bwd-bc-mp-row "1 2"

# --------------------------------------------------------------------------
if [ "${FAILURES}" -ne 0 ]; then
  echo
  echo "${FAILURES} check(s) failed"
  exit 1
fi

echo
echo "All checks passed"

# clean up
rm -rf "${TEST_RUNS_ROOT}"
