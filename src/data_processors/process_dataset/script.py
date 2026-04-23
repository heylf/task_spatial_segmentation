import sys
import random
import numpy as np
import anndata as ad
import scanpy as sc
import openproblems as op
import json

## VIASH START
par = {
    'input_sp': 'resources_test/task_spatial_segmentation/mouse_brain_combined/common_ist.zarr',
    'input_sc': 'resources_test/task_spatial_segmentation/mouse_brain_combined/common_scrnaseq.h5ad',
    'output_spatial_dataset': 'resources_test/task_spatial_segmentation/mouse_brain_combined/output_spatial_dataset.zarr',
    'output_scrnaseq': 'resources_test/task_spatial_segmentation/mouse_brain_combined/output_scrnaseq.h5ad',
    'method': 'xenium',
    'seed': 123,
    'config': 'task_spatial_segmentation/src/data_processors/process_dataset/config/config_default.json'
}

meta = {
    'resources_dir': 'target/executable/data_processors/process_dataset',
    'config': 'target/executable/data_processors/process_dataset/.config.vsh.yaml'
}
## VIASH END

# import helper functions
sys.path.append(meta['resources_dir'])

config = op.project.read_viash_config(meta["config"])

# set seed if need be
if par["seed"]:
    print(f">> Setting seed to {par['seed']}")
    random.seed(par["seed"])

print(">> Load data", flush=True)
adata = ad.read_h5ad(par["input_sc"])
print("input_sc:", adata)

print(f">> Process {par['method']} data")

if par['config']:
    print(f">> Perform standard data preprocessing")
    with open(par['config'], "r") as f:
        config = json.load(f)

    # Add config to params
    for key, value in config.items():
        setattr(par, key, value)

    adata.layers["counts"] = adata.X.copy()
    
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)
    adata.layers['normlog'] = adata.X
    
    sc.pp.highly_variable_genes(
        adata,
        flavor="seurat_v3",
        layer="counts",
        span=par['span'],
        n_top_genes=par['n_top_genes']
    )

    adata.var.sort_values("means")
    sc.pp.scale(adata, zero_center=False)
    adata.layers['normlogscale'] = adata.X
    
    adata.X = adata.layers['counts']

    # cell area normalization
    sc.pp.calculate_qc_metrics(adata, inplace=True)
    for x in ['transcript_counts', 'n_genes_by_counts']:
        adata.obs[f'canorm_{x}'] = adata.obs[f'{x}'] / adata.obs['cell_area']

print(">> Writing data", flush=True)
adata.write_h5ad(par["output_scrnaseq"])
