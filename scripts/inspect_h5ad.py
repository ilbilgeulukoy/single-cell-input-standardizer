import sys
import anndata as ad
import numpy as np
from scipy.sparse import issparse


def main(filepath: str) -> None:
    adata = ad.read_h5ad(filepath)

    print("=== AnnData object ===")
    print(adata)

    print("\n=== Shape ===")
    print("Cells / n_obs:", adata.n_obs)
    print("Genes / n_vars:", adata.n_vars)

    print("\n=== obs columns ===")
    print(list(adata.obs.columns))

    print("\n=== var columns ===")
    print(list(adata.var.columns))

    print("\n=== First obs rows ===")
    print(adata.obs.head())

    print("\n=== First var rows ===")
    print(adata.var.head())

    print("\n=== Sample counts ===")
    if "sample_id" in adata.obs.columns:
        print(adata.obs["sample_id"].value_counts())

    print("\n=== Patient counts ===")
    if "patient_id" in adata.obs.columns:
        print(adata.obs["patient_id"].value_counts())

    print("\n=== NaN check ===")
    if issparse(adata.X):
        print("NaN in X:", np.isnan(adata.X.data).sum())
    else:
        print("NaN in X:", np.isnan(adata.X).sum())

    print("NaN in obs:", adata.obs.isna().values.any())
    print("NaN in var:", adata.var.isna().values.any())


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python scripts/inspect_h5ad.py path/to/file.h5ad")
        sys.exit(1)

    main(sys.argv[1])
