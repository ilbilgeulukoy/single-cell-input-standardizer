from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
import gzip
import csv

import numpy as np
import pandas as pd
import anndata as ad
from scipy import sparse

from src.metadata_checker import check_cell_metadata_coverage
from src.sample_metadata_standardizer import apply_curated_sample_metadata_to_obs
from src.validation import (
    check_obs_var_alignment,
    check_duplicate_names,
    check_nan_sparse_safe,
    merge_validation_results,
)


@dataclass
class SparseGlobalWriteResult:
    output_h5ad: str
    shape: tuple[int, int]
    n_cells: int
    n_genes: int
    nnz: int
    obs_columns: list[str]
    var_columns: list[str]
    validation_passed: bool
    validation_errors: list[str]
    validation_warnings: list[str]
    coverage_warnings: list[str] = field(default_factory=list)
    metadata_standardization_passed: bool | None = None
    metadata_standardization_errors: list[str] = field(default_factory=list)
    metadata_standardization_warnings: list[str] = field(default_factory=list)
    evidence: dict[str, Any] = field(default_factory=dict)


def _open_text(path: str | Path):
    path = Path(path)
    if path.name.endswith(".gz"):
        return gzip.open(path, "rt", newline="")
    return open(path, "rt", newline="")


def count_global_csv_shape(matrix_path: str | Path) -> tuple[int, int, list[str], list[str]]:
    matrix_path = Path(matrix_path)

    with _open_text(matrix_path) as handle:
        reader = csv.reader(handle)
        header = next(reader)

        cell_ids = [str(x).strip().strip('"').strip("'") for x in header[1:]]
        first_genes = []
        n_genes = 0

        for row in reader:
            if not row:
                continue
            n_genes += 1
            if len(first_genes) < 5:
                first_genes.append(row[0])

    return len(cell_ids), n_genes, cell_ids, first_genes


def read_global_gene_by_cell_csv_sparse(
    matrix_path: str | Path,
    dtype=np.float32,
    report_every: int = 5000,
) -> tuple[sparse.csr_matrix, list[str], list[str], dict[str, Any]]:
    matrix_path = Path(matrix_path)

    with _open_text(matrix_path) as handle:
        reader = csv.reader(handle)
        header = next(reader)

        cell_ids = [str(x).strip().strip('"').strip("'") for x in header[1:]]
        n_cells = len(cell_ids)

        data = []
        row_indices = []
        col_indices = []
        gene_names = []

        for gene_idx, row in enumerate(reader):
            if not row:
                continue

            gene_name = str(row[0]).strip().strip('"').strip("'")
            values = row[1:]

            if len(values) != n_cells:
                raise ValueError(
                    f"Row {gene_idx} has {len(values)} values, expected {n_cells}."
                )

            gene_names.append(gene_name)

            for cell_idx, raw_value in enumerate(values):
                if raw_value in {"", "0", "0.0"}:
                    continue

                value = float(raw_value)
                if value != 0:
                    row_indices.append(cell_idx)
                    col_indices.append(gene_idx)
                    data.append(value)

            if report_every and (gene_idx + 1) % report_every == 0:
                print(f"Processed genes: {gene_idx + 1:,}; nnz so far: {len(data):,}")

    matrix = sparse.csr_matrix(
        (
            np.asarray(data, dtype=dtype),
            (
                np.asarray(row_indices, dtype=np.int64),
                np.asarray(col_indices, dtype=np.int64),
            ),
        ),
        shape=(n_cells, len(gene_names)),
    )

    evidence = {
        "matrix_path": str(matrix_path),
        "n_cells": n_cells,
        "n_genes": len(gene_names),
        "nnz": int(matrix.nnz),
        "density": float(matrix.nnz / (matrix.shape[0] * matrix.shape[1])) if matrix.shape[0] and matrix.shape[1] else 0.0,
    }

    return matrix, cell_ids, gene_names, evidence


def attach_metadata_to_obs(
    obs: pd.DataFrame,
    metadata_path: str | Path,
    cell_id_column: str,
    sample_id_column: str | None = None,
    keep_metadata_cells_only: bool = False,
) -> tuple[pd.DataFrame, list[str], Any]:
    metadata_path = Path(metadata_path)
    compression = "gzip" if metadata_path.name.endswith(".gz") else None
    sep = "," if metadata_path.name.endswith((".csv", ".csv.gz")) else "\t"

    metadata = pd.read_csv(metadata_path, sep=sep, compression=compression)
    metadata[cell_id_column] = metadata[cell_id_column].astype(str)
    metadata = metadata.drop_duplicates(subset=[cell_id_column]).set_index(cell_id_column)

    coverage = check_cell_metadata_coverage(
        matrix_cell_ids=list(obs.index),
        metadata_cell_ids=list(metadata.index),
        require_all_matrix_cells_in_metadata=False,
    )

    if keep_metadata_cells_only:
        common_cells = obs.index.intersection(metadata.index)
        obs = obs.loc[common_cells].copy()

    common_cells = obs.index.intersection(metadata.index)

    for column in metadata.columns:
        obs[column] = pd.NA
        if len(common_cells) > 0:
            obs.loc[common_cells, column] = metadata.loc[common_cells, column].values

    if sample_id_column and sample_id_column in obs.columns:
        obs["sample_id"] = obs[sample_id_column].astype(str)

    return obs, list(obs.index), coverage


def write_global_gene_by_cell_csv_to_h5ad_sparse(
    matrix_path: str | Path,
    output_h5ad: str | Path,
    metadata_path: str | Path | None = None,
    cell_id_column: str | None = None,
    sample_id_column: str | None = None,
    keep_metadata_cells_only: bool = False,
    obs_mappings: dict[str, str] | None = None,
    curated_sample_metadata: dict[str, dict[str, Any]] | None = None,
    standardize_obs_schema: bool = False,
    keep_extra_obs_columns: bool = True,
    dtype=np.float32,
) -> SparseGlobalWriteResult:
    output_h5ad = Path(output_h5ad)

    matrix, cell_ids, gene_names, evidence = read_global_gene_by_cell_csv_sparse(
        matrix_path,
        dtype=dtype,
    )

    obs = pd.DataFrame(index=pd.Index(cell_ids, name="cell_id"))
    var = pd.DataFrame(index=pd.Index(gene_names, name="gene_symbol"))

    coverage = None
    metadata_standardization = None

    if metadata_path is not None:
        if cell_id_column is None:
            raise ValueError("cell_id_column is required when metadata_path is provided.")

        obs, kept_cells, coverage = attach_metadata_to_obs(
            obs=obs,
            metadata_path=metadata_path,
            cell_id_column=cell_id_column,
            sample_id_column=sample_id_column,
            keep_metadata_cells_only=keep_metadata_cells_only,
        )

        if keep_metadata_cells_only:
            cell_position = {cell: idx for idx, cell in enumerate(cell_ids)}
            selected_indices = [cell_position[cell] for cell in kept_cells]
            matrix = matrix[selected_indices, :]

    if standardize_obs_schema:
        if curated_sample_metadata is None:
            raise ValueError(
                "curated_sample_metadata is required when standardize_obs_schema=True."
            )

        if sample_id_column is None and "sample_id" not in obs.columns:
            raise ValueError(
                "sample_id_column or obs['sample_id'] is required for standard obs schema."
            )

        metadata_standardization_sample_column = "sample_id" if "sample_id" in obs.columns else sample_id_column

        obs, metadata_standardization = apply_curated_sample_metadata_to_obs(
            obs=obs,
            curated_sample_metadata=curated_sample_metadata,
            sample_id_column=metadata_standardization_sample_column,
            keep_extra_obs_columns=keep_extra_obs_columns,
        )

    adata = ad.AnnData(X=matrix, obs=obs, var=var)

    validation = merge_validation_results(
        [
            check_obs_var_alignment(
                adata.X,
                obs_names=list(adata.obs_names),
                var_names=list(adata.var_names),
            ),
            check_duplicate_names(list(adata.obs_names), axis_label="obs", allow_duplicates=False),
            check_duplicate_names(list(adata.var_names), axis_label="var", allow_duplicates=True),
            check_nan_sparse_safe(adata.X),
        ]
    )

    output_h5ad.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(output_h5ad)

    return SparseGlobalWriteResult(
        output_h5ad=str(output_h5ad),
        shape=tuple(int(x) for x in adata.shape),
        n_cells=int(adata.n_obs),
        n_genes=int(adata.n_vars),
        nnz=int(adata.X.nnz) if sparse.issparse(adata.X) else int(np.count_nonzero(adata.X)),
        obs_columns=list(adata.obs.columns),
        var_columns=list(adata.var.columns),
        validation_passed=validation.passed,
        validation_errors=validation.errors,
        validation_warnings=validation.warnings,
        coverage_warnings=[] if coverage is None else coverage.warnings,
        metadata_standardization_passed=None if metadata_standardization is None else metadata_standardization.passed,
        metadata_standardization_errors=[] if metadata_standardization is None else metadata_standardization.errors,
        metadata_standardization_warnings=[] if metadata_standardization is None else metadata_standardization.warnings,
        evidence={
            **evidence,
            "metadata_standardization": None if metadata_standardization is None else metadata_standardization.evidence,
        },
    )
