#!/bin/bash
set -eo pipefail

JOB_NUM=$1
INPUT_FILE=$2
OUTPUT_FILE=$3
COMPACT_FILE=$4

NEVENTS=1000
SKIP=$(( JOB_NUM * NEVENTS ))

echo "=== job ${JOB_NUM} on $(hostname) at $(date) ==="
echo "compact: ${COMPACT_FILE}   input: ${INPUT_FILE}   skip: ${SKIP}"

set +e
source /opt/setup_mucoll.sh
set -e

which ddsim

ddsim \
    --inputFile      "${INPUT_FILE}" \
    --steeringFile   sim_steer_GEN_CONDOR.py \
    --compactFile    "${COMPACT_FILE}" \
    --skipNEvents    "${SKIP}" \
    --numberOfEvents "${NEVENTS}" \
    --random.seed    $(( JOB_NUM + 1 )) \
    --outputFile     "${OUTPUT_FILE}"

echo "=== done; wrote ${OUTPUT_FILE} ==="
ls -lh "${OUTPUT_FILE}"
