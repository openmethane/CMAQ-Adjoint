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
TEST_DATA_URL="${TEST_DATA_URL:-https://openmethane.s3.amazonaws.com/tests/cmaq-adjoint/cmaq-adjoint-test-data.20260804.tar.gz}"

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

# The bundle gained fwd.env (renamed from test.env), bwd.env and the adjoint
# forcing file when the backward test was added. Check for them up front so a
# stale cached bundle gives a clear message rather than a confusing model
# failure part way through.
for REQUIRED in fwd.env bwd.env FORCE.20221207.nc; do
  if [ ! -f "${DATA_DIR}/${REQUIRED}" ]; then
    echo "ERROR: ${DATA_DIR}/${REQUIRED} is missing."
    echo "The test data bundle predates the backward test. Delete"
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
