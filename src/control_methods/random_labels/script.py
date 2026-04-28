
import anndata as ad
import random
import pandas as pd

## VIASH START
par = {
    "input": "resources_test/task_spatial_segmentation/mouse_brain_combined/output_scrnaseq_reference.h5ad",
    "output": "resources_test/task_spatial_segmentation/mouse_brain_combined/reference_prediction.h5ad",
    "seed": 123,
    "label": "cell_type"
}
meta = {
    "name": "random_labels",
}
## VIASH END

if par["seed"]:
    print(f">> Setting seed to {par['seed']}")
    random.seed(par["seed"])

print("Load input data", flush=True)
input = ad.read_h5ad(par["input"])

print("Create random labels", flush=True)
input.obs[par["label"]] = [random.randint(1, 10) for _ in range(input.n_obs)]

print("Create output AnnData", flush=True)
output = ad.AnnData(
    obs=pd.DataFrame(input.obs[par["label"]]),
    uns={
        "dataset_id": input.uns["dataset_id"],
        "normalization_id": input.uns["normalization_id"],
        "method_id": meta["name"],
    },
)

print("Write output to file", flush=True)
output.write_h5ad(par["output"], compression="gzip")
