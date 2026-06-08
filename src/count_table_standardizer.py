from __future__ import annotations

import os
from typing import Any

import pandas as pd
import numpy as np
from anndata import AnnData
from scipy.sparse import issparse


def read_count_table(filepath: str, compression: str | None = "infer") -> pd.DataFrame:
    """
    Read a CSV/TSV count table.

    The separator is inferred from the file extension:
    - .csv or .csv.gz -> comma
    - .tsv, .txt, .tsv.gz, .txt.gz -> tab
    """
    basename = os.path.basename(filepath)

    if basename.endswith(".csv") or basename.endswith(".csv.gz"):
        sep = ","
    elif (
        basename.endswith(".tsv")
        or basename.endswith(".txt")
        or basename.endswith(".tsv.gz")
        or basename.endswith(".txt.gz")
    ):
        sep = "\t"
    else:
        sep = ","

    return pd.read_csv(filepath, sep=sep, compression=compression)


def infer_count_table_structure(df: pd.DataFrame) -> dict[str, Any]:
    """
    Infer basic structure of a count table.

    Current supported pattern:
    - first column contains gene symbols or gene IDs
    - remaining columns contain cell barcodes
    - raw matrix orientation is genes x cells
    """
    if df.empty:
        return {
            "is_valid": False,
            "reason": "empty_dataframe",
        }

    if df.shape[1] < 2:
        return {
            "is_valid": False,
            "reason": "less_than_two_columns",
        }

    first_column = df.columns[0]
    numeric_part = df.iloc[:, 1:]

    numeric_column_count = sum(
        pd.api.types.is_numeric_dtype(numeric_part[col])
        for col in numeric_part.columns
    )

    numeric_ratio = numeric_column_count / max(numeric_part.shape[1], 1)

    duplicated_feature_count = df.iloc[:, 0].duplicated().sum()

    likely_gene_by_cell = numeric_ratio > 0.95

    return {
        "is_valid": True,
        "first_column": str(first_column),
        "first_column_role": "gene_or_feature_id",
        "remaining_columns_role": "cells_or_barcodes",
        "raw_orientation": "genes_x_cells" if likely_gene_by_cell else "unknown",
        "requires_transpose_for_anndata": bool(likely_gene_by_cell),
        "numeric_ratio_after_first_column": float(numeric_ratio),
        "duplicated_feature_count": int(duplicated_feature_count),
        "n_rows": int(df.shape[0]),
        "n_columns": int(df.shape[1]),
    }


def validate_gene_by_cell_count_table(df: pd.DataFrame, filepath: str | None = None) -> None:
    """
    Validate a gene-by-cell count table before AnnData creation.
    """
    label = filepath or "input dataframe"

    if df.empty:
        raise ValueError(f"Input dataframe is empty: {label}")

    if df.shape[1] < 2:
        raise ValueError(
            f"Input dataframe should contain one gene column and at least one cell column: {label}"
        )

    numeric_part = df.iloc[:, 1:]
    non_numeric_columns = numeric_part.columns[
        ~numeric_part.apply(lambda col: pd.api.types.is_numeric_dtype(col))
    ]

    if len(non_numeric_columns) > 0:
        raise ValueError(
            f"Non-numeric count columns detected in {label}: {list(non_numeric_columns[:10])}"
        )


def build_anndata_from_gene_by_cell_table(
    df: pd.DataFrame,
    sample_id: str,
    sample_metadata: dict[str, Any] | None = None,
    source_file: str | None = None,
    gene_column_name: str = "gene",
) -> AnnData:
    """
    Build an AnnData object from a gene-by-cell count table.

    Input dataframe:
        rows    = genes
        columns = first column + cells

    AnnData:
        rows    = cells
        columns = genes

    Therefore:
        AnnData(df.T)
    """
    validate_gene_by_cell_count_table(df, filepath=source_file)

    original_first_column = df.columns[0]
    df = df.rename(columns={original_first_column: gene_column_name})
    df = df.set_index(gene_column_name)

    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)

    adata = AnnData(df.T)

    adata.obs_names = df.columns.astype(str)
    adata.var_names = df.index.astype(str)

    adata.var["gene_symbol"] = adata.var_names
    adata.obs["sample_id"] = sample_id

    if sample_metadata:
        for key, value in sample_metadata.items():
            adata.obs[key] = value

    if source_file:
        adata.obs["source_file"] = os.path.basename(source_file)

    return adata


def align_common_genes(list_adata: list[AnnData]) -> list[AnnData]:
    """
    Keep only genes shared across all AnnData objects.
    """
    if not list_adata:
        raise ValueError("No AnnData objects provided.")

    common_genes = list_adata[0].var_names

    for adata_sample in list_adata[1:]:
        common_genes = common_genes.intersection(adata_sample.var_names)

    return [adata_sample[:, common_genes].copy() for adata_sample in list_adata]


def check_nan_in_adata(adata: AnnData) -> dict[str, Any]:
    """
    Return NaN validation summary for X, obs, and var.
    """
    if issparse(adata.X):
        nan_count_x = int(np.isnan(adata.X.data).sum())
    else:
        nan_count_x = int(np.isnan(adata.X).sum())

    return {
        "nan_count_x": nan_count_x,
        "nan_exists_obs": bool(adata.obs.isna().values.any()),
        "nan_exists_var": bool(adata.var.isna().values.any()),
    }
