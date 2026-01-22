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

# path to the adjoint fwd binary being tested
ADJ_FWD_BIN="${ADJ_FWD_BIN:-/opt/cmaq/bin/ADJOINT_FWD}"

# url to the test data bundle
TEST_DATA_URL="${TEST_DATA_URL:-https://openmethane.s3.amazonaws.com/tests/cmaq-adj/cmaq-adj-test-data.tar.gz}"

# path to this script and its parent folder
TESTS_ROOT="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
TEST_RUNS_ROOT="/tmp/docker-cmaq-adj-tests"
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

# run a single process test
prepare_test_run "${TEST_RUNS_ROOT}/001-sp"
read_environment_vars "${TEST_DATA_DIR}/test.env"
pushd "${RUN_DIR}"
${ADJ_FWD_BIN}
cat "${LOGFILE}"
popd

# run a multi process test
prepare_test_run "${TEST_RUNS_ROOT}/002-mp"
read_environment_vars "${TEST_DATA_DIR}/test.env"
pushd "${RUN_DIR}"
NPCOL_NPROW="2 1" mpirun -np 2 ${ADJ_FWD_BIN}
cat "${LOGFILE}"
popd

# clean up
rm -rf "${TEST_RUNS_ROOT}"
