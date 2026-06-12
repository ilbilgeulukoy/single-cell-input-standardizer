from __future__ import annotations

import json
from pathlib import Path

import anndata as ad

from src.global_sparse_writer import write_global_gene_by_cell_csv_to_h5ad_sparse
from src.recipe_loader import validate_expected_obs_columns


# =============================================================================
# 1. Set script variables
# =============================================================================

DATASET_ID = 'GSE200997'
INPUT_FORMAT = 'single_global_compressed_csv_count_matrix'
READER = 'sparse_global_csv'

MATRIX_PATH = Path('data/independent_gse200997/raw/GSE200997_GEO_processed_CRC_10X_raw_UMI_count_matrix.csv.gz')
METADATA_PATH = Path('data/independent_gse200997/raw/GSE200997_GEO_processed_CRC_10X_cell_annotation.csv.gz') if 'data/independent_gse200997/raw/GSE200997_GEO_processed_CRC_10X_cell_annotation.csv.gz' else None

OUTPUT_H5AD = Path('generated_cases/GSE200997/output/GSE200997_standardized.h5ad')
SUMMARY_JSON = Path('generated_cases/GSE200997/output/GSE200997_standardized.summary.json')


# =============================================================================
# 2. Define input interpretation
# =============================================================================

CELL_ID_COLUMN = 'Unnamed: 0'
SAMPLE_ID_COLUMN = 'samples'
OBS_MAPPINGS = {'sample_id': 'samples'}

KEEP_METADATA_CELLS_ONLY = True
KEEP_EXTRA_OBS_COLUMNS = False
STANDARDIZE_OBS_SCHEMA = True


# =============================================================================
# 3. Expected output checks
# =============================================================================

EXPECTED = {'n_cells': 49859,
 'n_genes': 23828,
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


# =============================================================================
# 4. Curated sample metadata
#
# Human-curated fields may intentionally contain REVIEW_REQUIRED values.
# These values are preserved to avoid unsafe biological assumptions.
# =============================================================================

CURATED_SAMPLE_METADATA = {'B_cac10': {'patient_id': 'B_cac10',
             'dataset_id': 'GSE200997',
             'cancer_type': 'Colorectal',
             'tumor_site': 'REVIEW_REQUIRED_FROM_CONDITION',
             'metastasis_site': 'Unknown',
             'tumor_treatment': 'Unknown',
             'cancer_site_origin': 'Colon',
             'tumour_grade': 'Unknown',
             'tumour_stage': 'Unknown',
             'histological_subtype': 'Colorectal',
             'patient_ethnicity': 'Unknown'},
 'B_cac11': {'patient_id': 'B_cac11',
             'dataset_id': 'GSE200997',
             'cancer_type': 'Colorectal',
             'tumor_site': 'REVIEW_REQUIRED_FROM_CONDITION',
             'metastasis_site': 'Unknown',
             'tumor_treatment': 'Unknown',
             'cancer_site_origin': 'Colon',
             'tumour_grade': 'Unknown',
             'tumour_stage': 'Unknown',
             'histological_subtype': 'Colorectal',
             'patient_ethnicity': 'Unknown'},
 'B_cac14': {'patient_id': 'B_cac14',
             'dataset_id': 'GSE200997',
             'cancer_type': 'Colorectal',
             'tumor_site': 'REVIEW_REQUIRED_FROM_CONDITION',
             'metastasis_site': 'Unknown',
             'tumor_treatment': 'Unknown',
             'cancer_site_origin': 'Colon',
             'tumour_grade': 'Unknown',
             'tumour_stage': 'Unknown',
             'histological_subtype': 'Colorectal',
             'patient_ethnicity': 'Unknown'},
 'B_cac15': {'patient_id': 'B_cac15',
             'dataset_id': 'GSE200997',
             'cancer_type': 'Colorectal',
             'tumor_site': 'REVIEW_REQUIRED_FROM_CONDITION',
             'metastasis_site': 'Unknown',
             'tumor_treatment': 'Unknown',
             'cancer_site_origin': 'Colon',
             'tumour_grade': 'Unknown',
             'tumour_stage': 'Unknown',
             'histological_subtype': 'Colorectal',
             'patient_ethnicity': 'Unknown'},
 'B_cac4': {'patient_id': 'B_cac4',
            'dataset_id': 'GSE200997',
            'cancer_type': 'Colorectal',
            'tumor_site': 'REVIEW_REQUIRED_FROM_CONDITION',
            'metastasis_site': 'Unknown',
            'tumor_treatment': 'Unknown',
            'cancer_site_origin': 'Colon',
            'tumour_grade': 'Unknown',
            'tumour_stage': 'Unknown',
            'histological_subtype': 'Colorectal',
            'patient_ethnicity': 'Unknown'},
 'B_cac6': {'patient_id': 'B_cac6',
            'dataset_id': 'GSE200997',
            'cancer_type': 'Colorectal',
            'tumor_site': 'REVIEW_REQUIRED_FROM_CONDITION',
            'metastasis_site': 'Unknown',
            'tumor_treatment': 'Unknown',
            'cancer_site_origin': 'Colon',
            'tumour_grade': 'Unknown',
            'tumour_stage': 'Unknown',
            'histological_subtype': 'Colorectal',
            'patient_ethnicity': 'Unknown'},
 'B_cac7': {'patient_id': 'B_cac7',
            'dataset_id': 'GSE200997',
            'cancer_type': 'Colorectal',
            'tumor_site': 'REVIEW_REQUIRED_FROM_CONDITION',
            'metastasis_site': 'Unknown',
            'tumor_treatment': 'Unknown',
            'cancer_site_origin': 'Colon',
            'tumour_grade': 'Unknown',
            'tumour_stage': 'Unknown',
            'histological_subtype': 'Colorectal',
            'patient_ethnicity': 'Unknown'},
 'T_cac1': {'patient_id': 'T_cac1',
            'dataset_id': 'GSE200997',
            'cancer_type': 'Colorectal',
            'tumor_site': 'REVIEW_REQUIRED_FROM_CONDITION',
            'metastasis_site': 'Unknown',
            'tumor_treatment': 'Unknown',
            'cancer_site_origin': 'Colon',
            'tumour_grade': 'Unknown',
            'tumour_stage': 'Unknown',
            'histological_subtype': 'Colorectal',
            'patient_ethnicity': 'Unknown'},
 'T_cac10': {'patient_id': 'T_cac10',
             'dataset_id': 'GSE200997',
             'cancer_type': 'Colorectal',
             'tumor_site': 'REVIEW_REQUIRED_FROM_CONDITION',
             'metastasis_site': 'Unknown',
             'tumor_treatment': 'Unknown',
             'cancer_site_origin': 'Colon',
             'tumour_grade': 'Unknown',
             'tumour_stage': 'Unknown',
             'histological_subtype': 'Colorectal',
             'patient_ethnicity': 'Unknown'},
 'T_cac11': {'patient_id': 'T_cac11',
             'dataset_id': 'GSE200997',
             'cancer_type': 'Colorectal',
             'tumor_site': 'REVIEW_REQUIRED_FROM_CONDITION',
             'metastasis_site': 'Unknown',
             'tumor_treatment': 'Unknown',
             'cancer_site_origin': 'Colon',
             'tumour_grade': 'Unknown',
             'tumour_stage': 'Unknown',
             'histological_subtype': 'Colorectal',
             'patient_ethnicity': 'Unknown'},
 'T_cac12': {'patient_id': 'T_cac12',
             'dataset_id': 'GSE200997',
             'cancer_type': 'Colorectal',
             'tumor_site': 'REVIEW_REQUIRED_FROM_CONDITION',
             'metastasis_site': 'Unknown',
             'tumor_treatment': 'Unknown',
             'cancer_site_origin': 'Colon',
             'tumour_grade': 'Unknown',
             'tumour_stage': 'Unknown',
             'histological_subtype': 'Colorectal',
             'patient_ethnicity': 'Unknown'},
 'T_cac13': {'patient_id': 'T_cac13',
             'dataset_id': 'GSE200997',
             'cancer_type': 'Colorectal',
             'tumor_site': 'REVIEW_REQUIRED_FROM_CONDITION',
             'metastasis_site': 'Unknown',
             'tumor_treatment': 'Unknown',
             'cancer_site_origin': 'Colon',
             'tumour_grade': 'Unknown',
             'tumour_stage': 'Unknown',
             'histological_subtype': 'Colorectal',
             'patient_ethnicity': 'Unknown'},
 'T_cac14': {'patient_id': 'T_cac14',
             'dataset_id': 'GSE200997',
             'cancer_type': 'Colorectal',
             'tumor_site': 'REVIEW_REQUIRED_FROM_CONDITION',
             'metastasis_site': 'Unknown',
             'tumor_treatment': 'Unknown',
             'cancer_site_origin': 'Colon',
             'tumour_grade': 'Unknown',
             'tumour_stage': 'Unknown',
             'histological_subtype': 'Colorectal',
             'patient_ethnicity': 'Unknown'},
 'T_cac15': {'patient_id': 'T_cac15',
             'dataset_id': 'GSE200997',
             'cancer_type': 'Colorectal',
             'tumor_site': 'REVIEW_REQUIRED_FROM_CONDITION',
             'metastasis_site': 'Unknown',
             'tumor_treatment': 'Unknown',
             'cancer_site_origin': 'Colon',
             'tumour_grade': 'Unknown',
             'tumour_stage': 'Unknown',
             'histological_subtype': 'Colorectal',
             'patient_ethnicity': 'Unknown'},
 'T_cac16': {'patient_id': 'T_cac16',
             'dataset_id': 'GSE200997',
             'cancer_type': 'Colorectal',
             'tumor_site': 'REVIEW_REQUIRED_FROM_CONDITION',
             'metastasis_site': 'Unknown',
             'tumor_treatment': 'Unknown',
             'cancer_site_origin': 'Colon',
             'tumour_grade': 'Unknown',
             'tumour_stage': 'Unknown',
             'histological_subtype': 'Colorectal',
             'patient_ethnicity': 'Unknown'},
 'T_cac2': {'patient_id': 'T_cac2',
            'dataset_id': 'GSE200997',
            'cancer_type': 'Colorectal',
            'tumor_site': 'REVIEW_REQUIRED_FROM_CONDITION',
            'metastasis_site': 'Unknown',
            'tumor_treatment': 'Unknown',
            'cancer_site_origin': 'Colon',
            'tumour_grade': 'Unknown',
            'tumour_stage': 'Unknown',
            'histological_subtype': 'Colorectal',
            'patient_ethnicity': 'Unknown'},
 'T_cac3': {'patient_id': 'T_cac3',
            'dataset_id': 'GSE200997',
            'cancer_type': 'Colorectal',
            'tumor_site': 'REVIEW_REQUIRED_FROM_CONDITION',
            'metastasis_site': 'Unknown',
            'tumor_treatment': 'Unknown',
            'cancer_site_origin': 'Colon',
            'tumour_grade': 'Unknown',
            'tumour_stage': 'Unknown',
            'histological_subtype': 'Colorectal',
            'patient_ethnicity': 'Unknown'},
 'T_cac4': {'patient_id': 'T_cac4',
            'dataset_id': 'GSE200997',
            'cancer_type': 'Colorectal',
            'tumor_site': 'REVIEW_REQUIRED_FROM_CONDITION',
            'metastasis_site': 'Unknown',
            'tumor_treatment': 'Unknown',
            'cancer_site_origin': 'Colon',
            'tumour_grade': 'Unknown',
            'tumour_stage': 'Unknown',
            'histological_subtype': 'Colorectal',
            'patient_ethnicity': 'Unknown'},
 'T_cac5': {'patient_id': 'T_cac5',
            'dataset_id': 'GSE200997',
            'cancer_type': 'Colorectal',
            'tumor_site': 'REVIEW_REQUIRED_FROM_CONDITION',
            'metastasis_site': 'Unknown',
            'tumor_treatment': 'Unknown',
            'cancer_site_origin': 'Colon',
            'tumour_grade': 'Unknown',
            'tumour_stage': 'Unknown',
            'histological_subtype': 'Colorectal',
            'patient_ethnicity': 'Unknown'},
 'T_cac6': {'patient_id': 'T_cac6',
            'dataset_id': 'GSE200997',
            'cancer_type': 'Colorectal',
            'tumor_site': 'REVIEW_REQUIRED_FROM_CONDITION',
            'metastasis_site': 'Unknown',
            'tumor_treatment': 'Unknown',
            'cancer_site_origin': 'Colon',
            'tumour_grade': 'Unknown',
            'tumour_stage': 'Unknown',
            'histological_subtype': 'Colorectal',
            'patient_ethnicity': 'Unknown'},
 'T_cac7': {'patient_id': 'T_cac7',
            'dataset_id': 'GSE200997',
            'cancer_type': 'Colorectal',
            'tumor_site': 'REVIEW_REQUIRED_FROM_CONDITION',
            'metastasis_site': 'Unknown',
            'tumor_treatment': 'Unknown',
            'cancer_site_origin': 'Colon',
            'tumour_grade': 'Unknown',
            'tumour_stage': 'Unknown',
            'histological_subtype': 'Colorectal',
            'patient_ethnicity': 'Unknown'},
 'T_cac8': {'patient_id': 'T_cac8',
            'dataset_id': 'GSE200997',
            'cancer_type': 'Colorectal',
            'tumor_site': 'REVIEW_REQUIRED_FROM_CONDITION',
            'metastasis_site': 'Unknown',
            'tumor_treatment': 'Unknown',
            'cancer_site_origin': 'Colon',
            'tumour_grade': 'Unknown',
            'tumour_stage': 'Unknown',
            'histological_subtype': 'Colorectal',
            'patient_ethnicity': 'Unknown'},
 'T_cac9': {'patient_id': 'T_cac9',
            'dataset_id': 'GSE200997',
            'cancer_type': 'Colorectal',
            'tumor_site': 'REVIEW_REQUIRED_FROM_CONDITION',
            'metastasis_site': 'Unknown',
            'tumor_treatment': 'Unknown',
            'cancer_site_origin': 'Colon',
            'tumour_grade': 'Unknown',
            'tumour_stage': 'Unknown',
            'histological_subtype': 'Colorectal',
            'patient_ethnicity': 'Unknown'}}


# =============================================================================
# 5. Utility functions
# =============================================================================

def write_json(path: Path, data: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, indent=2), encoding="utf-8")


def check_input_files() -> None:
    if not MATRIX_PATH.exists():
        raise FileNotFoundError(f"Matrix file not found: {MATRIX_PATH}")

    if METADATA_PATH is not None and not METADATA_PATH.exists():
        raise FileNotFoundError(f"Metadata file not found: {METADATA_PATH}")


def validate_result(result) -> list[str]:
    validation_errors = list(result.validation_errors)

    expected_cells = EXPECTED.get("n_cells")
    expected_genes = EXPECTED.get("n_genes")
    expected_obs_columns = EXPECTED.get("obs_columns")

    if expected_cells is not None and expected_cells != "REVIEW_REQUIRED":
        if int(expected_cells) != result.n_cells:
            validation_errors.append(
                f"Expected {expected_cells} cells, found {result.n_cells}."
            )

    if expected_genes is not None and expected_genes != "REVIEW_REQUIRED":
        if int(expected_genes) != result.n_genes:
            validation_errors.append(
                f"Expected {expected_genes} genes, found {result.n_genes}."
            )

    missing_obs_columns = validate_expected_obs_columns(
        result.obs_columns,
        expected_obs_columns,
    )

    if missing_obs_columns:
        validation_errors.append(
            f"Expected obs columns missing from h5ad: {missing_obs_columns}"
        )

    return validation_errors


def write_summary(result, validation_errors: list[str]) -> None:
    summary = {
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
    }

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
