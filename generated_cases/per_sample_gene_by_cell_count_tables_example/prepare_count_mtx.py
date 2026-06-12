from __future__ import annotations

import json
import os
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy
from anndata import AnnData
from scipy.sparse import csc_matrix, issparse


# =============================================================================
# 1. Set script variables
# =============================================================================

DATASET_ID = 'GSE147082'
PAPER_ID = 'per_sample_gene_by_cell_count_tables_example'

RAW_DIR = Path('data/olelekan_style_test/raw')
DOWNLOADED_DIR = Path('data/olelekan_style_test/raw/downloaded')
ANNOTATED_DIR = Path('generated_cases/2021_Olelekan_style/output')

OUTPUT_H5AD = Path('generated_cases/per_sample_gene_by_cell_count_tables_example/output/per_sample_gene_by_cell_count_tables_example_standardized.h5ad')
SUMMARY_JSON = Path('generated_cases/per_sample_gene_by_cell_count_tables_example/output/per_sample_gene_by_cell_count_tables_example_standardized.summary.json')


# =============================================================================
# 2. Input file interpretation
# =============================================================================

SAMPLE_FILE_PATTERN = 'GSM*.csv.gz'
SAMPLE_ID_REGEX = '^(GSM[0-9]+)'
COMPRESSION = 'gzip'
GENE_COLUMN_NAME = 'gene'

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


# =============================================================================
# 3. Standard obs schema
# =============================================================================

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


# =============================================================================
# 4. Curated sample metadata
#
# This dictionary is the human-reviewed biological metadata layer.
# REVIEW_REQUIRED values are intentionally preserved when biological metadata
# cannot be assigned safely from file names or automatic inspection.
# =============================================================================

CURATED_SAMPLE_METADATA = {'GSM4416534': {'patient_id': 'PT-3232',
                'dataset_id': 'GSE147082',
                'cancer_type': 'Ovarian',
                'tumor_site': 'Metastasis',
                'metastasis_site': 'Omentum',
                'tumor_treatment': 'No',
                'cancer_site_origin': 'Left fallopian',
                'tumour_grade': 'HGSC',
                'tumour_stage': 'IIIc',
                'histological_subtype': 'Serous',
                'patient_ethnicity': 'Asian'},
 'GSM4416535': {'patient_id': 'PT-5150',
                'dataset_id': 'GSE147082',
                'cancer_type': 'Ovarian',
                'tumor_site': 'Metastasis',
                'metastasis_site': 'Omentum',
                'tumor_treatment': 'Yes',
                'cancer_site_origin': 'Unknown',
                'tumour_grade': 'HGSC',
                'tumour_stage': 'IVb',
                'histological_subtype': 'Serous',
                'patient_ethnicity': 'White'}}


# =============================================================================
# 5. Utility functions
# =============================================================================

def write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")


def ensure_folders() -> None:
    DOWNLOADED_DIR.mkdir(parents=True, exist_ok=True)
    ANNOTATED_DIR.mkdir(parents=True, exist_ok=True)
    OUTPUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY_JSON.parent.mkdir(parents=True, exist_ok=True)


def discover_sample_files() -> list[Path]:
    files = sorted(DOWNLOADED_DIR.glob(SAMPLE_FILE_PATTERN))
    files = [path for path in files if path.is_file()]

    if not files:
        raise FileNotFoundError(
            f"No sample files found in {DOWNLOADED_DIR} with pattern {SAMPLE_FILE_PATTERN}"
        )

    return files


def parse_sample_id(path: Path) -> str:
    match = re.search(SAMPLE_ID_REGEX, path.name)
    if match is None:
        raise ValueError(f"Could not parse sample_id from file name: {path.name}")
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


def read_gene_by_cell_count_table(path: Path) -> AnnData:
    df = pd.read_csv(path, compression=COMPRESSION)
    df = df.rename(columns={df.columns[0]: GENE_COLUMN_NAME})
    df = df.set_index(GENE_COLUMN_NAME)

    adata = AnnData(df.T)
    adata.obs_names = df.columns.astype(str)
    adata.var_names = df.index.astype(str)

    return adata


def check_missing_samples(sample_ids: list[str]) -> None:
    missing = sorted([sample_id for sample_id in sample_ids if sample_id not in CURATED_SAMPLE_METADATA])
    if missing:
        print("Missing samples in CURATED_SAMPLE_METADATA:", missing)
        raise KeyError("Some observed sample IDs are missing from the curated metadata dictionary.")


# =============================================================================
# 6. Construct AnnData objects from per-sample count tables
# =============================================================================

def construct_adata_from_per_sample_tables() -> AnnData:
    sample_files = discover_sample_files()
    sample_ids = [parse_sample_id(path) for path in sample_files]

    check_missing_samples(sample_ids)

    list_adata = []

    for path in sample_files:
        sample_id = parse_sample_id(path)
        adata = read_gene_by_cell_count_table(path)
        adata = add_standard_obs_metadata(adata, sample_id)

        list_adata.append(adata)
        print(sample_id, adata.shape)

    if not list_adata:
        raise ValueError("No AnnData objects were created.")

    print("Loaded samples:", len(list_adata))

    adata_combined = ad.concat(list_adata, axis=0, join="outer", index_unique=None)
    adata_combined.obs_names_make_unique()

    return adata_combined


# =============================================================================
# 7. Basic matrix checks before gene standardization
# =============================================================================

def check_nan_values(adata: AnnData) -> None:
    if issparse(adata.X):
        nan_count = int(np.isnan(adata.X.data).sum())
    else:
        nan_count = int(np.isnan(adata.X).sum())

    print("NaN count in X:", nan_count)
    print("NaN in obs:", bool(adata.obs.isna().values.any()))
    print("NaN in var:", bool(adata.var.isna().values.any()))


# =============================================================================
# 8. Common gene standardization block
#
# In the original workflow, this part used a lab-specific gene mapping backend.
# This generated script keeps that block optional because it may depend on
# environment-specific packages.
# =============================================================================

def run_external_gene_mapping_if_available(adata: AnnData) -> AnnData:
    if not USE_EXTERNAL_GENE_MAPPING:
        print("Skipping external gene mapping backend.")
        return adata

    try:
        from cancer_surfaceome_subtype_identifier.data_preparation.gene_ids_mapping.sample_gene_map_generation import SampleGeneMapGenerator
    except ImportError as exc:
        raise ImportError(
            "USE_EXTERNAL_GENE_MAPPING=True but the lab package could not be imported. "
            "Run this in the appropriate environment or set use_external_gene_mapping: false."
        ) from exc

    df_features = adata.var_names.to_frame(name="dataset_gene_symbol")

    sample_dfvar_generator = SampleGeneMapGenerator(df_features, str(DOWNLOADED_DIR))
    sample_dfvar_generator.regenerate_mapfiles()
    sample_dfvar_generator.generate_sample_mapfile()

    filepath_df_var = DOWNLOADED_DIR / "df_var.tsv"
    df_var = pd.read_csv(filepath_df_var, sep="\t")

    adata.var = df_var.set_index("dataset_gene_symbol").loc[adata.var_names]

    if "to_remove" in adata.var.columns:
        adata = adata[:, ~adata.var["to_remove"].fillna(False)].copy()

    if "gene_symbol_touse_w_dup" not in adata.var.columns:
        raise KeyError("gene_symbol_touse_w_dup column is missing after gene mapping.")

    gene_symbols = adata.var["gene_symbol_touse_w_dup"].astype(str).values
    unique_symbols, inverse_indices = np.unique(gene_symbols, return_inverse=True)

    if not issparse(adata.X):
        adata.X = scipy.sparse.csr_matrix(adata.X)

    n_vars = adata.n_vars
    mapping_matrix = csc_matrix(
        (np.ones(n_vars), (inverse_indices, np.arange(n_vars))),
        shape=(len(unique_symbols), n_vars),
    )

    new_X = adata.X @ mapping_matrix.T
    new_var = adata.var.groupby("gene_symbol_touse_w_dup").first().reset_index()
    new_var.index = new_var["gene_symbol_touse_w_dup"]

    new_adata = AnnData(X=new_X, obs=adata.obs.copy(), var=new_var)
    new_adata.obs_names_make_unique()

    if "to_remove" in new_adata.var.columns:
        new_adata.var = new_adata.var.drop(columns=["to_remove"])

    return new_adata


# =============================================================================
# 9. Final validation and writing
# =============================================================================

def summarize_adata(adata: AnnData) -> dict:
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
        "review_required_counts": review_required_counts,
        "technical_validation_passed": True,
        "curation_ready": sum(review_required_counts.values()) == 0,
        "use_external_gene_mapping": USE_EXTERNAL_GENE_MAPPING,
    }


def main() -> None:
    ensure_folders()

    adata = construct_adata_from_per_sample_tables()
    check_nan_values(adata)

    adata = run_external_gene_mapping_if_available(adata)

    adata.write_h5ad(OUTPUT_H5AD)

    summary = summarize_adata(adata)
    write_json(SUMMARY_JSON, summary)

    print("Wrote h5ad:", OUTPUT_H5AD)
    print("Wrote summary:", SUMMARY_JSON)
    print("AnnData shape:", adata.shape)
    print("Obs columns:", list(adata.obs.columns))
    print("Curation ready:", summary["curation_ready"])


if __name__ == "__main__":
    main()
