#!/usr/bin/env bash
# Generate nextflow_schema with updated basecaller enumerations
#
# This script uses `nextflow config` to obtain the basecaller container,
# creates JSON arrays of the models using the container's list-models script
# and injects them with jq to create nextflow_schema.json.new.
set -euo pipefail

TARGET=$1
DORADO_CONTAINER=$(nextflow config -flat | grep "process.'withLabel:wf_basecalling'.container" | awk -F'= ' '{print $2}' | sed "s,',,g")
echo "# DORADO_CONTAINER=${DORADO_CONTAINER}"

# Convert model lists to JSON arrays
SIMPLEX_MODELS=$(singularity exec "docker://${DORADO_CONTAINER}" list-models --simplex --only-names | jq -Rn '[inputs]')
MODBASE_MODELS=$(singularity exec "docker://${DORADO_CONTAINER}" list-models --modbase --only-names | jq -Rn '[inputs]')

# Inject JSON arrays to relevant schema enum
jq \
    --indent 4 \
    --argjson simplex_models "${SIMPLEX_MODELS}" \
    --argjson modbase_models "${MODBASE_MODELS}" \
    '(.definitions.basecalling_options.properties.basecaller_cfg.enum) = $simplex_models |
    (.definitions.basecalling_options.properties.remora_cfg.enum) = $modbase_models' \
    ${TARGET}/nextflow_schema.json > ${TARGET}/nextflow_schema.json.new

# Remove newline at end of file
truncate -s -1 ${TARGET}/nextflow_schema.json.new

echo "# Updated schema generated, you should inspect it before adopting it!"
echo "diff ${TARGET}/nextflow_schema.json.new ${TARGET}/nextflow_schema.json"
echo "mv ${TARGET}/nextflow_schema.json.new ${TARGET}/nextflow_schema.json"
