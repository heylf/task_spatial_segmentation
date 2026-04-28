
import random
import anndata as ad
import scanpy as sc
import pandas as pd

## VIASH START
par = {
    'input': 'resources_test/task_spatial_segmentation/mouse_brain_combined/output_scrnaseq_reference.h5ad',
    'output': 'resources_test/task_spatial_segmentation/mouse_brain_combined/method_prediction.h5ad',
    'label': 'cell_type',
    'n_neighbors': 20,
    'min_dist': 0.1,
    'spread': 1.2,
    'resolution': 1.0,
    'seed': 123
}
## VIASH END

# set seed if need be
if par["seed"]:
    print(f">> Setting seed to {par['seed']}")
    random.seed(par["seed"])

print('>> Reading input files', flush=True)
input = ad.read_h5ad(par['input'])

print('>> Perform Leiden clustering', flush=True)
sc.pp.neighbors(input, n_neighbors=par['n_neighbors'], random_state=par['seed'])
sc.tl.umap(input, min_dist=par['min_dist'], spread=par['spread'], random_state=par['seed'])
sc.tl.leiden(input, resolution=par['resolution'], key_added=par["label"], random_state=par['seed'])

print(">> Write output AnnData to file", flush=True)
output = ad.AnnData(
    obs=pd.DataFrame(input.obs[par["label"]]),
    uns={
        "dataset_id": input.uns["dataset_id"],
        "normalization_id": input.uns["normalization_id"],
        #"method_id": input.uns["method_id"], #TODO
    },
)

output.write_h5ad(par['output'], compression='gzip')
