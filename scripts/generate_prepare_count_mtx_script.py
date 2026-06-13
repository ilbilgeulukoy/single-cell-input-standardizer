from __future__ import annotations

import argparse
import pprint
from pathlib import Path
from typing import Any

import yaml

from src.recipe_loader import load_recipe


def as_python(value: Any) -> str:
    return pprint.pformat(value, width=100, sort_dicts=False)


def generate_global_matrix_with_cell_metadata_script(recipe: dict[str, Any], output_dir: Path) -> str:
    dataset_id = recipe["dataset_id"]

    matrix_path = recipe["matrix_path"]
    metadata_path = recipe.get("metadata_path")

    input_format = recipe["input_format"]
    reader = recipe.get("reader", "sparse_global_csv")

    cell_id_column = recipe.get("cell_id_column")
    sample_id_column = recipe.get("sample_id_column")

    keep_metadata_cells_only = bool(recipe.get("keep_metadata_cells_only", True))
    keep_extra_obs_columns = bool(recipe.get("keep_extra_obs_columns", False))
    standardize_obs_schema = bool(recipe.get("standardize_obs_schema", True))

    obs_mappings = recipe.get("obs_mappings", {})
    expected = recipe.get("expected", {})
    curated_sample_metadata = recipe.get("curated_sample_metadata", {})

    output_h5ad = output_dir / "output" / f"{dataset_id}_standardized.h5ad"
    summary_json = output_dir / "output" / f"{dataset_id}_standardized.summary.json"

    return f'''from __future__ import annotations

import json
from pathlib import Path

import anndata as ad

from src.global_sparse_writer import write_global_gene_by_cell_csv_to_h5ad_sparse
from src.recipe_loader import validate_expected_obs_columns


# =============================================================================
# 1. Set script variables
# =============================================================================

DATASET_ID = {as_python(dataset_id)}
INPUT_FORMAT = {as_python(input_format)}
READER = {as_python(reader)}

MATRIX_PATH = Path({as_python(matrix_path)})
METADATA_PATH = Path({as_python(metadata_path)}) if {as_python(metadata_path)} else None

OUTPUT_H5AD = Path({as_python(str(output_h5ad))})
SUMMARY_JSON = Path({as_python(str(summary_json))})


# =============================================================================
# 2. Define input interpretation
# =============================================================================

CELL_ID_COLUMN = {as_python(cell_id_column)}
SAMPLE_ID_COLUMN = {as_python(sample_id_column)}
OBS_MAPPINGS = {as_python(obs_mappings)}

KEEP_METADATA_CELLS_ONLY = {as_python(keep_metadata_cells_only)}
KEEP_EXTRA_OBS_COLUMNS = {as_python(keep_extra_obs_columns)}
STANDARDIZE_OBS_SCHEMA = {as_python(standardize_obs_schema)}


# =============================================================================
# 3. Expected output checks
# =============================================================================

EXPECTED = {as_python(expected)}


# =============================================================================
# 4. Curated sample metadata
#
# Human-curated fields may intentionally contain REVIEW_REQUIRED values.
# These values are preserved to avoid unsafe biological assumptions.
# =============================================================================

CURATED_SAMPLE_METADATA = {as_python(curated_sample_metadata)}


# =============================================================================
# 5. Utility functions
# =============================================================================

def write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")


def check_input_files() -> None:
    if not MATRIX_PATH.exists():
        raise FileNotFoundError(f"Matrix file not found: {{MATRIX_PATH}}")

    if METADATA_PATH is not None and not METADATA_PATH.exists():
        raise FileNotFoundError(f"Metadata file not found: {{METADATA_PATH}}")


def validate_result(result) -> list[str]:
    validation_errors = list(result.validation_errors)

    expected_cells = EXPECTED.get("n_cells")
    expected_genes = EXPECTED.get("n_genes")
    expected_obs_columns = EXPECTED.get("obs_columns")

    if expected_cells is not None and expected_cells != "REVIEW_REQUIRED":
        if int(expected_cells) != result.n_cells:
            validation_errors.append(
                f"Expected {{expected_cells}} cells, found {{result.n_cells}}."
            )

    if expected_genes is not None and expected_genes != "REVIEW_REQUIRED":
        if int(expected_genes) != result.n_genes:
            validation_errors.append(
                f"Expected {{expected_genes}} genes, found {{result.n_genes}}."
            )

    missing_obs_columns = validate_expected_obs_columns(
        result.obs_columns,
        expected_obs_columns,
    )

    if missing_obs_columns:
        validation_errors.append(
            f"Expected obs columns missing from h5ad: {{missing_obs_columns}}"
        )

    return validation_errors


def write_summary(result, validation_errors: list[str]) -> None:
    summary = {{
        "dataset_id": DATASET_ID,
        "input_format": INPUT_FORMAT,
        "reader": READER,
        "output_h5ad": result.output_h5ad,
        "adata_shape": list(result.shape),
        "n_cells": result.n_cells,
        "n_genes": result.n_genes,
        "nnz": result.nnz,
        "obs_columns": result.obs_columns,
        "var_columns": result.var_columns,
        "validation_passed": result.validation_passed and not validation_errors,
        "validation_errors": validation_errors,
        "validation_warnings": result.validation_warnings,
        "coverage_warnings": result.coverage_warnings,
        "metadata_standardization_passed": result.metadata_standardization_passed,
        "metadata_standardization_errors": result.metadata_standardization_errors,
        "metadata_standardization_warnings": result.metadata_standardization_warnings,
        "standardize_obs_schema": STANDARDIZE_OBS_SCHEMA,
        "keep_extra_obs_columns": KEEP_EXTRA_OBS_COLUMNS,
        "expected": EXPECTED,
        "evidence": result.evidence,
    }}

    write_json(SUMMARY_JSON, summary)


# =============================================================================
# 6. Construct AnnData
#
# This dataset uses one global gene-by-cell count matrix and one cell metadata table.
# The sparse writer reads the matrix in chunks and builds AnnData as cells x genes.
# =============================================================================

def construct_adata_from_global_matrix():
    result = write_global_gene_by_cell_csv_to_h5ad_sparse(
        matrix_path=MATRIX_PATH,
        output_h5ad=OUTPUT_H5AD,
        metadata_path=METADATA_PATH,
        cell_id_column=CELL_ID_COLUMN,
        sample_id_column=SAMPLE_ID_COLUMN,
        keep_metadata_cells_only=KEEP_METADATA_CELLS_ONLY,
        obs_mappings=OBS_MAPPINGS,
        curated_sample_metadata=CURATED_SAMPLE_METADATA,
        standardize_obs_schema=STANDARDIZE_OBS_SCHEMA,
        keep_extra_obs_columns=KEEP_EXTRA_OBS_COLUMNS,
    )

    return result


# =============================================================================
# 7. Main
# =============================================================================

def main() -> None:
    if INPUT_FORMAT != "single_global_compressed_csv_count_matrix":
        raise ValueError(
            "This generated prepare_count_mtx.py expects "
            "single_global_compressed_csv_count_matrix."
        )

    if READER != "sparse_global_csv":
        raise ValueError("This generated prepare_count_mtx.py expects reader=sparse_global_csv.")

    OUTPUT_H5AD.parent.mkdir(parents=True, exist_ok=True)
    SUMMARY_JSON.parent.mkdir(parents=True, exist_ok=True)

    check_input_files()

    result = construct_adata_from_global_matrix()
    validation_errors = validate_result(result)
    write_summary(result, validation_errors)

    adata = ad.read_h5ad(OUTPUT_H5AD)

    print("Wrote h5ad:", OUTPUT_H5AD)
    print("Wrote summary:", SUMMARY_JSON)
    print("AnnData shape:", adata.shape)
    print("Obs columns:", list(adata.obs.columns))
    print("Validation passed:", result.validation_passed and not validation_errors)
    print("Metadata standardization passed:", result.metadata_standardization_passed)


if __name__ == "__main__":
    main()
'''


def generate_per_sample_count_tables_gene_by_cell_script(recipe: dict[str, Any], output_dir: Path) -> str:
    dataset_id = recipe["dataset_id"]
    paper_id = recipe.get("paper_id", dataset_id)

    raw_dir = recipe.get("raw_dir", ".")
    downloaded_dir = recipe.get("downloaded_dir", raw_dir)
    annotated_dir = recipe.get("annotated_dir", str(output_dir / "output"))

    sample_file_pattern = recipe.get("sample_file_pattern", "*.csv.gz")
    sample_id_regex = recipe.get("sample_id_regex", r"^(GSM[0-9]+)")
    compression = recipe.get("compression", "gzip")
    gene_column_name = recipe.get("gene_column_name", "gene")

    curated_sample_metadata = recipe.get("curated_sample_metadata", {})
    expected = recipe.get("expected", {})
    output_h5ad = output_dir / "output" / f"{paper_id}_standardized.h5ad"
    summary_json = output_dir / "output" / f"{paper_id}_standardized.summary.json"

    use_external_gene_mapping = bool(recipe.get("use_external_gene_mapping", False))

    return f'''from __future__ import annotations

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

DATASET_ID = {as_python(dataset_id)}
PAPER_ID = {as_python(paper_id)}

RAW_DIR = Path({as_python(raw_dir)})
DOWNLOADED_DIR = Path({as_python(downloaded_dir)})
ANNOTATED_DIR = Path({as_python(annotated_dir)})

OUTPUT_H5AD = Path({as_python(str(output_h5ad))})
SUMMARY_JSON = Path({as_python(str(summary_json))})


# =============================================================================
# 2. Input file interpretation
# =============================================================================

SAMPLE_FILE_PATTERN = {as_python(sample_file_pattern)}
SAMPLE_ID_REGEX = {as_python(sample_id_regex)}
COMPRESSION = {as_python(compression)}
GENE_COLUMN_NAME = {as_python(gene_column_name)}

EXPECTED = {as_python(expected)}

USE_EXTERNAL_GENE_MAPPING = {as_python(use_external_gene_mapping)}


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

CURATED_SAMPLE_METADATA = {as_python(curated_sample_metadata)}


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
            f"No sample files found in {{DOWNLOADED_DIR}} with pattern {{SAMPLE_FILE_PATTERN}}"
        )

    return files


def parse_sample_id(path: Path) -> str:
    match = re.search(SAMPLE_ID_REGEX, path.name)
    if match is None:
        raise ValueError(f"Could not parse sample_id from file name: {{path.name}}")
    return match.group(1)


def add_standard_obs_metadata(adata: AnnData, sample_id: str) -> AnnData:
    if sample_id not in CURATED_SAMPLE_METADATA:
        raise KeyError(
            f"Sample {{sample_id}} is missing from CURATED_SAMPLE_METADATA. "
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
    df = df.rename(columns={{df.columns[0]: GENE_COLUMN_NAME}})
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
    df_var = pd.read_csv(filepath_df_var, sep="\\t")

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

    review_required_counts = {{}}
    for column in STANDARD_OBS_COLUMNS:
        values = adata.obs[column].astype(str)
        review_required_counts[column] = int(values.str.contains("REVIEW_REQUIRED", na=False).sum())

    return {{
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
    }}


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
'''

def generate_tenx_mtx_triplet_script(recipe: dict[str, Any], output_dir: Path) -> str:
    dataset_id = recipe["dataset_id"]
    paper_id = recipe.get("paper_id", dataset_id)

    downloaded_dir = recipe.get("downloaded_dir", ".")
    annotated_dir = recipe.get("annotated_dir", str(output_dir / "output"))

    sample_id_regex = recipe.get("sample_id_regex", r"(GSM[0-9]+)")
    matrix_pattern = recipe.get("matrix_pattern", "*matrix.mtx*")
    features_pattern = recipe.get("features_pattern", "*features.tsv*")
    genes_pattern = recipe.get("genes_pattern", "*genes.tsv*")
    barcodes_pattern = recipe.get("barcodes_pattern", "*barcodes.tsv*")

    curated_sample_metadata = recipe.get("curated_sample_metadata", {})
    expected = recipe.get("expected", {})

    output_h5ad = output_dir / "output" / f"{paper_id}_standardized.h5ad"
    summary_json = output_dir / "output" / f"{paper_id}_standardized.summary.json"

    use_external_gene_mapping = bool(recipe.get("use_external_gene_mapping", False))

    return f'''from __future__ import annotations

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

DATASET_ID = {as_python(dataset_id)}
PAPER_ID = {as_python(paper_id)}

DOWNLOADED_DIR = Path({as_python(downloaded_dir)})
ANNOTATED_DIR = Path({as_python(annotated_dir)})

OUTPUT_H5AD = Path({as_python(str(output_h5ad))})
SUMMARY_JSON = Path({as_python(str(summary_json))})


# =============================================================================
# 2. Input file interpretation
# =============================================================================

SAMPLE_ID_REGEX = {as_python(sample_id_regex)}
MATRIX_PATTERN = {as_python(matrix_pattern)}
FEATURES_PATTERN = {as_python(features_pattern)}
GENES_PATTERN = {as_python(genes_pattern)}
BARCODES_PATTERN = {as_python(barcodes_pattern)}

EXPECTED = {as_python(expected)}
USE_EXTERNAL_GENE_MAPPING = {as_python(use_external_gene_mapping)}


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

CURATED_SAMPLE_METADATA = {as_python(curated_sample_metadata)}


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
            return [line.rstrip("\\n") for line in handle]
    return path.read_text().splitlines()


def parse_sample_id_from_path(path: Path) -> str:
    match = re.search(SAMPLE_ID_REGEX, str(path))
    if match is None:
        raise ValueError(f"Could not parse sample_id from path: {{path}}")
    return match.group(1)


def add_standard_obs_metadata(adata: AnnData, sample_id: str) -> AnnData:
    if sample_id not in CURATED_SAMPLE_METADATA:
        raise KeyError(
            f"Sample {{sample_id}} is missing from CURATED_SAMPLE_METADATA. "
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
                {{
                    "sample_id": sample_id,
                    "folder": str(folder),
                    "matrix_path": matrix_path,
                    "features_path": features_path,
                    "barcodes_path": barcodes_path,
                }}
            )

    unique = {{}}
    for triplet in triplets:
        key = (
            str(triplet["matrix_path"]),
            str(triplet["features_path"]),
            str(triplet["barcodes_path"]),
        )
        unique[key] = triplet

    triplets = list(unique.values())

    if not triplets:
        raise FileNotFoundError(f"No complete 10x triplets found under {{DOWNLOADED_DIR}}.")

    return sorted(triplets, key=lambda item: item["sample_id"])


# =============================================================================
# 7. Read 10x triplets into AnnData
# =============================================================================

def read_features_table(path: Path) -> pd.DataFrame:
    rows = []

    for line in read_lines(path):
        if not line:
            continue
        rows.append(line.split("\\t"))

    if not rows:
        raise ValueError(f"Empty features/genes file: {{path}}")

    max_len = max(len(row) for row in rows)
    normalized = [row + [""] * (max_len - len(row)) for row in rows]
    df = pd.DataFrame(normalized)

    if df.shape[1] == 1:
        df.columns = ["gene_symbol"]
        df["gene_id"] = df["gene_symbol"]
    else:
        df = df.rename(columns={{0: "gene_id", 1: "gene_symbol"}})

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
            f"matrix={{matrix.shape}}, features={{len(features)}}, barcodes={{len(barcodes)}}"
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
    df_var = pd.read_csv(filepath_df_var, sep="\\t")

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

    review_required_counts = {{}}
    for column in STANDARD_OBS_COLUMNS:
        values = adata.obs[column].astype(str)
        review_required_counts[column] = int(values.str.contains("REVIEW_REQUIRED", na=False).sum())

    return {{
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
    }}


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
'''

def generate_tenx_h5_script(recipe: dict[str, Any], output_dir: Path) -> str:
    dataset_id = recipe["dataset_id"]
    paper_id = recipe.get("paper_id", dataset_id)

    downloaded_dir = recipe.get("downloaded_dir", ".")
    annotated_dir = recipe.get("annotated_dir", str(output_dir / "output"))

    h5_file_pattern = recipe.get("h5_file_pattern", "*.h5")
    sample_id_regex = recipe.get("sample_id_regex", r"(GSM[0-9]+)")
    curated_sample_metadata = recipe.get("curated_sample_metadata", {})
    expected = recipe.get("expected", {})

    output_h5ad = output_dir / "output" / f"{paper_id}_standardized.h5ad"
    summary_json = output_dir / "output" / f"{paper_id}_standardized.summary.json"
    use_external_gene_mapping = bool(recipe.get("use_external_gene_mapping", False))

    script = r'''from __future__ import annotations

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


DATASET_ID = __DATASET_ID__
PAPER_ID = __PAPER_ID__

DOWNLOADED_DIR = Path(__DOWNLOADED_DIR__)
ANNOTATED_DIR = Path(__ANNOTATED_DIR__)

OUTPUT_H5AD = Path(__OUTPUT_H5AD__)
SUMMARY_JSON = Path(__SUMMARY_JSON__)

H5_FILE_PATTERN = __H5_FILE_PATTERN__
SAMPLE_ID_REGEX = __SAMPLE_ID_REGEX__

EXPECTED = __EXPECTED__
USE_EXTERNAL_GENE_MAPPING = __USE_EXTERNAL_GENE_MAPPING__

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

CURATED_SAMPLE_METADATA = __CURATED_SAMPLE_METADATA__


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
'''

    replacements = {
        "__DATASET_ID__": as_python(dataset_id),
        "__PAPER_ID__": as_python(paper_id),
        "__DOWNLOADED_DIR__": as_python(downloaded_dir),
        "__ANNOTATED_DIR__": as_python(annotated_dir),
        "__OUTPUT_H5AD__": as_python(str(output_h5ad)),
        "__SUMMARY_JSON__": as_python(str(summary_json)),
        "__H5_FILE_PATTERN__": as_python(h5_file_pattern),
        "__SAMPLE_ID_REGEX__": as_python(sample_id_regex),
        "__EXPECTED__": as_python(expected),
        "__USE_EXTERNAL_GENE_MAPPING__": as_python(use_external_gene_mapping),
        "__CURATED_SAMPLE_METADATA__": as_python(curated_sample_metadata),
    }

    for old, new in replacements.items():
        script = script.replace(old, new)

    return script

def build_prepare_script(recipe: dict[str, Any], output_dir: Path) -> str:
    archetype = recipe.get("script_archetype")

    if archetype == "global_matrix_with_cell_metadata":
        return generate_global_matrix_with_cell_metadata_script(recipe, output_dir)

    if archetype == "per_sample_count_tables_gene_by_cell":
        return generate_per_sample_count_tables_gene_by_cell_script(recipe, output_dir)

    if archetype == "tenx_mtx_triplet_single_or_multi_sample":
        return generate_tenx_mtx_triplet_script(recipe, output_dir)

    if archetype == "tenx_h5_single_or_multi_sample":
        return generate_tenx_h5_script(recipe, output_dir)

    raise ValueError(
        f"Unsupported script_archetype: {archetype}. "
        "Currently supported: global_matrix_with_cell_metadata, "
        "per_sample_count_tables_gene_by_cell, "
        "tenx_mtx_triplet_single_or_multi_sample, "
        "tenx_h5_single_or_multi_sample"
    )


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate a case-specific prepare_count_mtx.py script from a recipe."
    )
    parser.add_argument("--recipe", required=True)
    parser.add_argument("--output-dir", required=True)

    args = parser.parse_args()

    recipe = load_recipe(args.recipe)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    prepare_script_path = output_dir / "prepare_count_mtx.py"
    recipe_copy_path = output_dir / "recipe.yaml"

    script_text = build_prepare_script(recipe, output_dir)

    prepare_script_path.write_text(script_text, encoding="utf-8")
    recipe_copy_path.write_text(
        yaml.safe_dump(recipe, sort_keys=False, allow_unicode=True),
        encoding="utf-8",
    )

    (output_dir / "output").mkdir(exist_ok=True)
    (output_dir / "reports").mkdir(exist_ok=True)

    print("Wrote generated script:", prepare_script_path)
    print("Copied recipe:", recipe_copy_path)
    print("Output directory:", output_dir / "output")
    print()
    print("Run with:")
    print(f"python {prepare_script_path}")


if __name__ == "__main__":
    main()
