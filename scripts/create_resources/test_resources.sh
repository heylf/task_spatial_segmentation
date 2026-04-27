#!/bin/bash

# get the root of the directory
REPO_ROOT=$(git rev-parse --show-toplevel)

# ensure that the command below is run from the root of the repository
cd "$REPO_ROOT"

# # remove this when you have implemented the script
# echo "TODO: replace the commands in this script with the sequence of components that you need to run to generate test_resources."
# echo "  Inside this script, you will need to place commands to generate example files for each of the 'src/api/file_*.yaml' files."
# exit 1

set -e

DATASET_ID=mouse_brain_combined

RAW_DATA=resources_test/common
DATASET_DIR=resources_test/task_spatial_segmentation/$DATASET_ID

mkdir -p $DATASET_DIR

# process dataset
viash run src/data_processors/process_dataset/config.vsh.yaml -- \
  --input_sp $RAW_DATA/2023_10x_mouse_brain_xenium_rep1/dataset.zarr \
  --input_sc $RAW_DATA/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad \
  --output_spatial_dataset $DATASET_DIR/spatial_dataset.zarr \
  --output_scrnaseq_reference $DATASET_DIR/scrnaseq_reference.h5ad

# run one method
viash run src/methods/cellpose/config.vsh.yaml -- \
    --input $DATASET_DIR/spatial_dataset.zarr \
    --output $DATASET_DIR/prediction.h5ad

# run one metric
# TODO: implement this!
# viash run src/metrics/ari/config.vsh.yaml -- \
#     --input_prediction $DATASET_DIR/prediction.h5ad \
#     --input_scrnaseq_reference $DATASET_DIR/scrnaseq_reference.h5ad \
#     --output $DATASET_DIR/score.h5ad

# write manual state.yaml. this is not actually necessary but you never know it might be useful
cat > $DATASET_DIR/state.yaml << HERE
id: $DATASET_ID
spatial_dataset: spatial_dataset.zarr
scrnaseq_reference: scrnaseq_reference.h5ad
prediction: prediction.h5ad
score: score.h5ad
HERE

# only run this if you have access to the openproblems-data bucket
aws s3 sync --profile op \
  "$DATASET_DIR" s3://openproblems-data/resources_test/task_spatial_segmentation/mouse_brain_combined/ \
  --delete --dryrun
