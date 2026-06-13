from __future__ import annotations

import json
import re
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp
from anndata import AnnData
from scipy.sparse import issparse


DATASET_ID = 'GSE_TEST_10X_H5'
PAPER_ID = 'tenx_h5_example'

DOWNLOADED_DIR = Path('data/test_tenx_h5_template')
ANNOTATED_DIR = Path('generated_cases/tenx_h5_example/output')

OUTPUT_H5AD = Path('generated_cases/tenx_h5_example/output/tenx_h5_example_standardized.h5ad')
SUMMARY_JSON = Path('generated_cases/tenx_h5_example/output/tenx_h5_example_standardized.summary.json')

H5_FILE_PATTERN = '*.h5'
SAMPLE_ID_REGEX = '(GSM[0-9]+)'

EXPECTED = {'n_cells': 'REVIEW_REQUIRED',
 'n_genes': 'REVIEW_REQUIRED',
 'obs_columns': ['sample_id',
                 'patient_id',
                 'dataset_id',
                 'cancer_type',
                 'tumor_site',
                 'metastasis_site',
                 'tumor_treatment',
                 'cancer_site_origin',
                 'tumour_grade',
                 'tumour_stage',
                 'histological_subtype',
                 'patient_ethnicity']}
USE_EXTERNAL_GENE_MAPPING = False

STANDARD_OBS_COLUMNS = [
    "sample_id",
    "patient_id",
    "dataset_id",
    "cancer_type",
    "tumor_site",
    "metastasis_site",
    "tumor_treatment",
    "cancer_site_origin",
    "tumour_grade",
    "tumour_stage",
    "histological_subtype",
    "patient_ethnicity",
]

CURATED_SAMPLE_METADATA = {'GSM1': {'patient_id': 'REVIEW_REQUIRED',
          'dataset_id': 'GSE_TEST_10X_H5',
          'cancer_type': 'REVIEW_REQUIRED',
          'tumor_site': 'REVIEW_REQUIRED',
          'metastasis_site': 'Unknown',
          'tumor_treatment': 'REVIEW_REQUIRED',
          'cancer_site_origin': 'REVIEW_REQUIRED',
          'tumour_grade': 'Unknown',
          'tumour_stage': 'Unknown',
          'histological_subtype': 'REVIEW_REQUIRED',
          'patient_ethnicity': 'Unknown'},
 'GSM2': {'patient_id': 'REVIEW_REQUIRED',
          'dataset_id': 'GSE_TEST_10X_H5',
          'cancer_type': 'REVIEW_REQUIRED',
          'tumor_site': 'REVIEW_REQUIRED',
          'metastasis_site': 'Unknown',
          'tumor_treatment': 'REVIEW_REQUIRED',
          'cancer_site_origin': 'REVIEW_REQUIRED',
          'tumour_grade': 'Unknown',
          'tumour_stage': 'Unknown',
          'histological_subtype': 'REVIEW_REQUIRED',
          'patient_ethnicity': 'Unknown'}}


def write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")


def ensure_folders() -> None:
    DOWNLOADED_DIR.mkdir(parents=True, exist_ok=True)
    ANNOTATED_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY_JSON.parent.mkdir(parents=True, exist_ok=True)


def decode_h5_array(values) -> list[str]:
    result = []
    for value in values:
        if isinstance(value, bytes):
            result.append(value.decode("utf-8"))
        else:
            result.append(str(value))
    return result


def parse_sample_id_from_path(path: Path) -> str:
    match = re.search(SAMPLE_ID_REGEX, path.name)
    if match is None:
        raise ValueError(f"Could not parse sample_id from h5 filename: {path.name}")
    return match.group(1)


def add_standard_obs_metadata(adata: AnnData, sample_id: str) -> AnnData:
    if sample_id not in CURATED_SAMPLE_METADATA:
        raise KeyError(
            f"Sample {sample_id} is missing from CURATED_SAMPLE_METADATA. "
            "Add it to the recipe before running this script."
        )

    metadata = CURATED_SAMPLE_METADATA[sample_id]

    adata.obs["sample_id"] = sample_id
    for column in STANDARD_OBS_COLUMNS:
        if column == "sample_id":
            continue
        adata.obs[column] = metadata.get(column, "Unknown")

    adata.obs = adata.obs[STANDARD_OBS_COLUMNS].copy()
    return adata


def discover_h5_files() -> list[Path]:
    files = sorted([path for path in DOWNLOADED_DIR.rglob(H5_FILE_PATTERN) if path.is_file()])
    if not files:
        raise FileNotFoundError(
            f"No h5 files found under {DOWNLOADED_DIR} with pattern {H5_FILE_PATTERN}"
        )
    return files


def detect_matrix_group(handle: h5py.File):
    if "matrix" in handle:
        return handle["matrix"], "modern_10x_h5_v3_matrix_group"

    candidates = []
    required = {"barcodes", "data", "indices", "indptr", "shape"}

    for key in handle.keys():
        item = handle[key]
        if isinstance(item, h5py.Group) and required.issubset(set(item.keys())):
            candidates.append(key)

    if len(candidates) == 1:
        group_name = candidates[0]
        return handle[group_name], f"old_10x_h5_v2_genome_group:{group_name}"

    if len(candidates) > 1:
        raise ValueError(f"Multiple possible old-style genome groups found: {candidates}")

    raise ValueError("Could not detect a supported 10x h5 layout.")


def read_feature_metadata(matrix_group) -> pd.DataFrame:
    if "features" in matrix_group:
        features = matrix_group["features"]
        gene_ids = decode_h5_array(features["id"][:])
        gene_symbols = decode_h5_array(features["name"][:])
    else:
        gene_ids = decode_h5_array(matrix_group["genes"][:])
        gene_symbols = decode_h5_array(matrix_group["gene_names"][:])

    return pd.DataFrame(
        {
            "dataset_gene_id": gene_ids,
            "dataset_gene_symbol": gene_symbols,
        }
    )


def aggregate_duplicate_genes(adata: AnnData) -> AnnData:
    gene_symbols = pd.Index(adata.var["dataset_gene_symbol"].astype(str))

    if not gene_symbols.duplicated().any():
        adata.var_names = gene_symbols
        return adata

    unique_symbols, inverse = np.unique(gene_symbols, return_inverse=True)

    mapping = sp.csc_matrix(
        (np.ones(len(inverse)), (inverse, np.arange(len(inverse)))),
        shape=(len(unique_symbols), len(inverse)),
    )

    if not issparse(adata.X):
        x = sp.csr_matrix(adata.X)
    else:
        x = adata.X.tocsr()

    new_x = x @ mapping.T

    new_var = (
        adata.var.assign(dataset_gene_symbol=gene_symbols)
        .groupby("dataset_gene_symbol", sort=False)
        .first()
    )
    new_var = new_var.loc[unique_symbols].copy()
    new_var.index = unique_symbols

    aggregated = AnnData(X=new_x, obs=adata.obs.copy(), var=new_var)
    aggregated.obs_names = adata.obs_names.copy()
    aggregated.var_names = pd.Index(unique_symbols).astype(str)
    return aggregated


def read_10x_h5_as_adata(path: Path, sample_id: str) -> tuple[AnnData, dict]:
    with h5py.File(path, "r") as handle:
        matrix_group, layout_label = detect_matrix_group(handle)

        data = matrix_group["data"][:]
        indices = matrix_group["indices"][:]
        indptr = matrix_group["indptr"][:]
        shape = tuple(matrix_group["shape"][:])

        barcodes = decode_h5_array(matrix_group["barcodes"][:])
        features = read_feature_metadata(matrix_group)

    gene_by_cell = sp.csc_matrix((data, indices, indptr), shape=shape)
    x = gene_by_cell.T.tocsr()

    if x.shape[0] != len(barcodes):
        raise ValueError(
            f"Barcode count does not match matrix cells for {path}: "
            f"matrix cells={x.shape[0]}, barcodes={len(barcodes)}"
        )

    if x.shape[1] != features.shape[0]:
        raise ValueError(
            f"Feature count does not match matrix genes for {path}: "
            f"matrix genes={x.shape[1]}, features={features.shape[0]}"
        )

    adata = AnnData(X=x)
    adata.obs_names = pd.Index([f"{sample_id}_{barcode}" for barcode in barcodes]).astype(str)
    adata.var_names = features["dataset_gene_symbol"].astype(str).values

    adata.var["dataset_gene_id"] = features["dataset_gene_id"].astype(str).values
    adata.var["dataset_gene_symbol"] = features["dataset_gene_symbol"].astype(str).values

    adata.obs_names_make_unique()
    adata = aggregate_duplicate_genes(adata)

    inspection = {
        "path": str(path),
        "layout_label": layout_label,
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "duplicated_gene_symbols": int(features["dataset_gene_symbol"].duplicated().sum()),
    }

    return adata, inspection


def check_missing_samples(sample_ids: list[str]) -> None:
    missing = sorted([sample_id for sample_id in sample_ids if sample_id not in CURATED_SAMPLE_METADATA])
    if missing:
        print("Missing samples in CURATED_SAMPLE_METADATA:", missing)
        raise KeyError("Some observed sample IDs are missing from the curated metadata dictionary.")


def construct_adata_from_10x_h5_files() -> tuple[AnnData, list[dict]]:
    h5_files = discover_h5_files()
    sample_ids = [parse_sample_id_from_path(path) for path in h5_files]

    check_missing_samples(sample_ids)

    list_adata = []
    inspections = []

    for path, sample_id in zip(h5_files, sample_ids):
        adata, inspection = read_10x_h5_as_adata(path, sample_id)
        adata = add_standard_obs_metadata(adata, sample_id)

        list_adata.append(adata)
        inspections.append(inspection)

        print(sample_id, adata.shape, inspection["layout_label"])

    if not list_adata:
        raise ValueError("No AnnData objects were created.")

    print("Loaded h5 samples:", len(list_adata))

    adata_combined = ad.concat(list_adata, axis=0, join="outer", index_unique=None)
    adata_combined.obs_names_make_unique()

    return adata_combined, inspections


def summarize_adata(adata: AnnData, inspections: list[dict]) -> dict:
    if issparse(adata.X):
        nnz = int(adata.X.nnz)
    else:
        nnz = int(np.count_nonzero(adata.X))

    review_required_counts = {}
    for column in STANDARD_OBS_COLUMNS:
        values = adata.obs[column].astype(str)
        review_required_counts[column] = int(values.str.contains("REVIEW_REQUIRED", na=False).sum())

    return {
        "dataset_id": DATASET_ID,
        "paper_id": PAPER_ID,
        "output_h5ad": str(OUTPUT_H5AD),
        "shape": list(adata.shape),
        "n_cells": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "nnz": nnz,
        "obs_columns": list(adata.obs.columns),
        "var_columns": list(adata.var.columns),
        "h5_inspections": inspections,
        "review_required_counts": review_required_counts,
        "technical_validation_passed": True,
        "curation_ready": sum(review_required_counts.values()) == 0,
        "use_external_gene_mapping": USE_EXTERNAL_GENE_MAPPING,
    }


def main() -> None:
    ensure_folders()

    adata, inspections = construct_adata_from_10x_h5_files()

    if USE_EXTERNAL_GENE_MAPPING:
        raise NotImplementedError(
            "External gene mapping is not enabled in this minimal 10x h5 generated template yet."
        )

    print("Skipping external gene mapping backend.")

    adata.write_h5ad(OUTPUT_H5AD)

    summary = summarize_adata(adata, inspections)
    write_json(SUMMARY_JSON, summary)

    print("Wrote h5ad:", OUTPUT_H5AD)
    print("Wrote summary:", SUMMARY_JSON)
    print("AnnData shape:", adata.shape)
    print("Obs columns:", list(adata.obs.columns))
    print("Curation ready:", summary["curation_ready"])


if __name__ == "__main__":
    main()
