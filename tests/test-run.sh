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

# path to this script and its parent folder
TESTS_ROOT="$( cd -- "$(dirname "$0")" >/dev/null 2>&1 ; pwd -P )"
TEST_RUNS_ROOT="/tmp/docker-cmaq-adj-tests"
mkdir -p "${TEST_RUNS_ROOT}"

# fetch test data too large to fit in the repo
TEST_DATA_DIR="${TESTS_ROOT}/test-data"
if [ ! -d "${TEST_DATA_DIR}" ]; then
  echo "Fetching test data"
  wget -c -nv https://openmethane.s3.amazonaws.com/tests/cmaq-adj/cmaq-adj-test-data.tar.gz \
    -O "${TEST_RUNS_ROOT}/test-data.tar.gz"
  mkdir -p "${TEST_DATA_DIR}"
  tar -xf "${TEST_RUNS_ROOT}/test-data.tar.gz" -C "${TEST_DATA_DIR}"
  rm "${TEST_RUNS_ROOT}/test-data.tar.gz"
fi

# run a single threaded test
DATA_DIR="${TEST_DATA_DIR}"
RUN_DIR="${TEST_RUNS_ROOT}/001-sp"
if [ -d "${RUN_DIR}" ]; then
  echo "Clearing previous run"
  rm -rf "${RUN_DIR}"
fi
mkdir -p "${RUN_DIR}/output"
mkdir -p "${RUN_DIR}/chkpnt"
set -o allexport && source "${TEST_DATA_DIR}/test.env" && set +o allexport

/opt/cmaq/BLD_fwd_CH4only/ADJOINT_FWD

cat "${LOGFILE}"

# clean up
rm -rf "${TEST_RUNS_ROOT}"
