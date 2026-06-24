#!/bin/bash
set -eo pipefail


# args
JOB_NUM=$1
INPUT_FILE=$2
OUTPUT_FILE=$3
COMPACT_FILE=$4


source ~/mucoll/.setup_shower_profile.sh

echo "Using compact file: ${COMPACT_FILE}"
echo "Job ${Job_NUM}"
echo "Input: ${INPUT_FILE}"
echo "Outpit: ${OUTPUT_FILE}"
which ddsim

NEVENTS = 1000
SKIP=$((JOB_NUM * NEVENTS))

ddsim \
	--inputFile ${INPUT_FILE} \
	--steeringFile sim_steer_condor.py \
	--outputFile ${OUTPUT_FILE} \
	--numberOfEvents ${NEVENTS} \
	--compactFile ${COMPACT_FILE} \
	--skipNEvents ${SKIP} \
	--random.seed $((JOB_NUM + 1)) \\
