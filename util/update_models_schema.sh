#!/usr/bin/env bash
# Generate nextflow_schema with updated basecaller enumerations
#
# This script uses `nextflow config` to obtain the basecaller container,
# creates JSON arrays of the models using the container's list-models script
# and injects them with jq to create nextflow_schema.json.new.
set -euo pipefail

TARGET=$1
ENGINE=$2

if ! command -v nextflow &> /dev/null
then
    # we should be in CI, nextflow is installed right here
    NEXTFLOW="./nextflow"
else
    NEXTFLOW=`which nextflow`
fi

# work out how to inspect the container contents
DORADO_CONTAINER=$(${NEXTFLOW} config -flat | grep "process.'withLabel:wf_basecalling'.container" | awk -F'= ' '{print $2}' | sed "s,',,g")
echo "# DORADO_CONTAINER=${DORADO_CONTAINER}"
if [ "$ENGINE" = "simg" ]; then
    CMD_PREFIX="singularity exec docker://${DORADO_CONTAINER}"
else
    CMD_PREFIX="docker run ${DORADO_CONTAINER}"
fi

# Convert model lists to JSON arrays
SIMPLEX_MODELS=$(${CMD_PREFIX} list-models --simplex --only-names | jq -Rn '[inputs]')
MODBASE_MODELS=$(${CMD_PREFIX} list-models --modbase --only-names | jq -Rn '[inputs]')

# Inject JSON arrays to relevant schema enum
jq \
    -j \
    --indent 4 \
    --argjson simplex_models "${SIMPLEX_MODELS}" \
    --argjson modbase_models "${MODBASE_MODELS}" \
    '(.definitions.basecalling_options.properties.basecaller_cfg.enum) = $simplex_models |
    (.definitions.basecalling_options.properties.remora_cfg.enum) = $modbase_models' \
    ${TARGET}/nextflow_schema.json > ${TARGET}/nextflow_schema.json.new

echo "# Updated schema generated, you should inspect it before adopting it!"
echo "diff ${TARGET}/nextflow_schema.json ${TARGET}/nextflow_schema.json.new"
echo "mv ${TARGET}/nextflow_schema.json.new ${TARGET}/nextflow_schema.json"
