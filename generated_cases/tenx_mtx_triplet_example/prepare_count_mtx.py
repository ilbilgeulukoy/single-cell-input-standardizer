from __future__ import annotations

import gzip
import json
import re
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy
from anndata import AnnData
from scipy.io import mmread
from scipy.sparse import csc_matrix, issparse


# =============================================================================
# 1. Set script variables
# =============================================================================

DATASET_ID = 'GSE_TEST_10X_MTX'
PAPER_ID = 'tenx_mtx_triplet_example'

DOWNLOADED_DIR = Path('data/test_tenx_mtx_standardizer')
ANNOTATED_DIR = Path('generated_cases/tenx_mtx_triplet_example/output')

OUTPUT_H5AD = Path('generated_cases/tenx_mtx_triplet_example/output/tenx_mtx_triplet_example_standardized.h5ad')
SUMMARY_JSON = Path('generated_cases/tenx_mtx_triplet_example/output/tenx_mtx_triplet_example_standardized.summary.json')


# =============================================================================
# 2. Input file interpretation
# =============================================================================

SAMPLE_ID_REGEX = '(GSM[0-9]+)'
MATRIX_PATTERN = '*matrix.mtx*'
FEATURES_PATTERN = '*features.tsv*'
GENES_PATTERN = '*genes.tsv*'
BARCODES_PATTERN = '*barcodes.tsv*'

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
# =============================================================================

CURATED_SAMPLE_METADATA = {'GSM1': {'patient_id': 'REVIEW_REQUIRED',
          'dataset_id': 'GSE_TEST_10X_MTX',
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
          'dataset_id': 'GSE_TEST_10X_MTX',
          'cancer_type': 'REVIEW_REQUIRED',
          'tumor_site': 'REVIEW_REQUIRED',
          'metastasis_site': 'Unknown',
          'tumor_treatment': 'REVIEW_REQUIRED',
          'cancer_site_origin': 'REVIEW_REQUIRED',
          'tumour_grade': 'Unknown',
          'tumour_stage': 'Unknown',
          'histological_subtype': 'REVIEW_REQUIRED',
          'patient_ethnicity': 'Unknown'}}


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


def read_lines(path: Path) -> list[str]:
    if path.name.endswith(".gz"):
        with gzip.open(path, "rt") as handle:
            return [line.rstrip("\n") for line in handle]
    return path.read_text().splitlines()


def parse_sample_id_from_path(path: Path) -> str:
    match = re.search(SAMPLE_ID_REGEX, str(path))
    if match is None:
        raise ValueError(f"Could not parse sample_id from path: {path}")
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


# =============================================================================
# 6. Discover 10x Matrix Market triplets
# =============================================================================

def find_first_matching_file(folder: Path, patterns: list[str]) -> Path | None:
    for pattern in patterns:
        matches = sorted([path for path in folder.glob(pattern) if path.is_file()])
        if matches:
            return matches[0]
    return None


def discover_10x_triplets() -> list[dict]:
    candidate_folders = [DOWNLOADED_DIR]
    candidate_folders.extend([path for path in DOWNLOADED_DIR.rglob("*") if path.is_dir()])

    triplets = []

    for folder in candidate_folders:
        matrix_path = find_first_matching_file(folder, [MATRIX_PATTERN])
        barcodes_path = find_first_matching_file(folder, [BARCODES_PATTERN])
        features_path = find_first_matching_file(folder, [FEATURES_PATTERN, GENES_PATTERN])

        if matrix_path and barcodes_path and features_path:
            sample_id = parse_sample_id_from_path(matrix_path)
            triplets.append(
                {
                    "sample_id": sample_id,
                    "folder": str(folder),
                    "matrix_path": matrix_path,
                    "features_path": features_path,
                    "barcodes_path": barcodes_path,
                }
            )

    unique = {}
    for triplet in triplets:
        key = (
            str(triplet["matrix_path"]),
            str(triplet["features_path"]),
            str(triplet["barcodes_path"]),
        )
        unique[key] = triplet

    triplets = list(unique.values())

    if not triplets:
        raise FileNotFoundError(f"No complete 10x triplets found under {DOWNLOADED_DIR}.")

    return sorted(triplets, key=lambda item: item["sample_id"])


# =============================================================================
# 7. Read 10x triplets into AnnData
# =============================================================================

def read_features_table(path: Path) -> pd.DataFrame:
    rows = []

    for line in read_lines(path):
        if not line:
            continue
        rows.append(line.split("\t"))

    if not rows:
        raise ValueError(f"Empty features/genes file: {path}")

    max_len = max(len(row) for row in rows)
    normalized = [row + [""] * (max_len - len(row)) for row in rows]
    df = pd.DataFrame(normalized)

    if df.shape[1] == 1:
        df.columns = ["gene_symbol"]
        df["gene_id"] = df["gene_symbol"]
    else:
        df = df.rename(columns={0: "gene_id", 1: "gene_symbol"})

    return df


def read_10x_mtx_triplet(matrix_path: Path, features_path: Path, barcodes_path: Path) -> AnnData:
    matrix = mmread(str(matrix_path)).tocsr()

    barcodes = read_lines(barcodes_path)
    features = read_features_table(features_path)

    if matrix.shape[0] == len(features) and matrix.shape[1] == len(barcodes):
        x = matrix.T.tocsr()
    elif matrix.shape[0] == len(barcodes) and matrix.shape[1] == len(features):
        x = matrix.tocsr()
    else:
        raise ValueError(
            "Matrix dimensions do not match features/barcodes. "
            f"matrix={matrix.shape}, features={len(features)}, barcodes={len(barcodes)}"
        )

    adata = AnnData(X=x)
    adata.obs_names = pd.Index(barcodes).astype(str)
    adata.var_names = features["gene_symbol"].astype(str).values

    adata.var["dataset_gene_id"] = features["gene_id"].astype(str).values
    adata.var["dataset_gene_symbol"] = features["gene_symbol"].astype(str).values

    adata.obs_names_make_unique()
    adata.var_names_make_unique()

    return adata


def check_missing_samples(sample_ids: list[str]) -> None:
    missing = sorted([sample_id for sample_id in sample_ids if sample_id not in CURATED_SAMPLE_METADATA])
    if missing:
        print("Missing samples in CURATED_SAMPLE_METADATA:", missing)
        raise KeyError("Some observed sample IDs are missing from the curated metadata dictionary.")


def construct_adata_from_10x_triplets() -> AnnData:
    triplets = discover_10x_triplets()
    sample_ids = [triplet["sample_id"] for triplet in triplets]

    check_missing_samples(sample_ids)

    list_adata = []

    for triplet in triplets:
        sample_id = triplet["sample_id"]

        adata = read_10x_mtx_triplet(
            matrix_path=triplet["matrix_path"],
            features_path=triplet["features_path"],
            barcodes_path=triplet["barcodes_path"],
        )

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
# 8. Optional external gene mapping backend
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

    adata = construct_adata_from_10x_triplets()
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
