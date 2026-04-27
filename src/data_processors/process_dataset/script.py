import random
import anndata as ad
import spatialdata as sd
import os
import shutil

## VIASH START
par = {
    'input_sp': 'resources_test/common/2023_10x_mouse_brain_xenium_rep1/dataset.zarr',
    'input_sc': 'resources_test/common/2023_yao_mouse_brain_scrnaseq_10xv2/dataset.h5ad',
    'output_spatial_dataset': 'resources_test/task_spatial_segmentation/mouse_brain_combined/output_spatial_dataset.zarr',
    'output_scrnaseq_reference': 'resources_test/task_spatial_segmentation/mouse_brain_combined/output_scrnaseq_reference.h5ad',
    'span': 0.3,
    'seed': 123,
    'n_top_genes': 3000
}
## VIASH END

# set seed if need be
if par["seed"]:
    print(f">> Setting seed to {par['seed']}")
    random.seed(par["seed"])

print(">> Load data", flush=True)
sc_data = ad.read_h5ad(par["input_sc"])

print(">> Processing sc_data", flush=True)

# TODO: process the single-cell dataset

print(f"single cell data: {sc_data}")

print(">> Writing data", flush=True)
sc_data.write_h5ad(par["output_scrnaseq_reference"], compression="gzip")

# read input_sp
print(">> Read spatial data", flush=True)
sp_data = sd.read_zarr(par["input_sp"])

print(">> Processing spatial data", flush=True)
# TODO: process the spatial dataset

print(f"spatial data: {sp_data}")

print(">> Writing spatial data", flush=True)
# remove directory if it exists
if os.path.exists(par["output_spatial_dataset"]):
    shutil.rmtree(par["output_spatial_dataset"])
sp_data.write(par["output_spatial_dataset"], overwrite=True)
