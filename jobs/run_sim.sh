#!/bin/bash

# args
JOB_NUM=$1
INPUT_FILE=$2
OUTPUT_FILE=$3

source ~/mucoll/.enterApptainer.sh

source ~/mucoll/.setup_shower_profile.sh

COMPACT_FILE=$4
echo "Using compact file: ${COMPACT_FILE}"

sed 's|SIM.compactFile = *|SIM.compactFile = \"${COMPACT_FILE}\"|' ~/mucoll/SteeringMacros/k4Reco/sim_steer_condor.py
echo "Job ${Job_NUM}"
echo "Input: ${INPUT_FILE}"
echo "Outpit: ${OUTPUT_FILE}"


ddsim \
	--inputFile ${INPUT_FILE} \
	--steeringFile sim_steer_condor.py \
	--outputFile ${OUTPUT_FILE} \
	--numberOfEvents 1000 \
	--compactFile ${COMPACT_FILE}
