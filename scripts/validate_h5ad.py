from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import anndata as ad
from scipy import sparse

from src.standard_schema import STANDARD_OBS_COLUMNS


def validate_h5ad(path: Path) -> dict[str, Any]:
    adata = ad.read_h5ad(path)

    obs_columns = list(adata.obs.columns)
    missing_standard_obs_columns = [
        column for column in STANDARD_OBS_COLUMNS if column not in obs_columns
    ]

    placeholder_counts = {}
    for column in STANDARD_OBS_COLUMNS:
        if column in adata.obs.columns:
            values = adata.obs[column].astype(str)
            placeholder_counts[column] = int(values.str.contains("REVIEW_REQUIRED", na=False).sum())

    result = {
        "h5ad_path": str(path),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "shape": [int(adata.n_obs), int(adata.n_vars)],
        "is_sparse": bool(sparse.issparse(adata.X)),
        "nnz": int(adata.X.nnz) if sparse.issparse(adata.X) else None,
        "density": (
            float(adata.X.nnz / (adata.n_obs * adata.n_vars))
            if sparse.issparse(adata.X) and adata.n_obs and adata.n_vars
            else None
        ),
        "obs_columns": obs_columns,
        "var_columns": list(adata.var.columns),
        "standard_obs_columns": STANDARD_OBS_COLUMNS,
        "missing_standard_obs_columns": missing_standard_obs_columns,
        "duplicate_obs_names": int(adata.obs_names.duplicated().sum()),
        "duplicate_var_names": int(adata.var_names.duplicated().sum()),
        "review_required_counts": placeholder_counts,
        "technical_validation_passed": (
            len(missing_standard_obs_columns) == 0
            and int(adata.obs_names.duplicated().sum()) == 0
        ),
        "curation_ready": sum(placeholder_counts.values()) == 0,
        "passed": (
            len(missing_standard_obs_columns) == 0
            and int(adata.obs_names.duplicated().sum()) == 0
        ),
    }

    return result


def main() -> None:
    parser = argparse.ArgumentParser(description="Validate a standardized AnnData h5ad output.")
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--output-json", required=True)

    args = parser.parse_args()

    result = validate_h5ad(Path(args.h5ad))

    output_json = Path(args.output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    output_json.write_text(json.dumps(result, indent=2), encoding="utf-8")

    print("Validated:", args.h5ad)
    print("Shape:", tuple(result["shape"]))
    print("Sparse:", result["is_sparse"])
    print("Missing standard obs columns:", result["missing_standard_obs_columns"])
    print("Duplicate obs names:", result["duplicate_obs_names"])
    print("Duplicate var names:", result["duplicate_var_names"])
    print("Review-required counts:", result["review_required_counts"])
    print("Technical validation passed:", result["technical_validation_passed"])
    print("Curation ready:", result["curation_ready"])
    print("Passed:", result["passed"])


if __name__ == "__main__":
    main()
