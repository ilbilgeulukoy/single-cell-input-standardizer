from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping
import re

import pandas as pd
from scipy import sparse

from src.sample_id_parser import summarize_barcode_suffixes
from src.metadata_checker import check_cell_metadata_coverage
from src.validation import (
    check_obs_var_alignment,
    check_duplicate_names,
    check_nan_sparse_safe,
    merge_validation_results,
)


@dataclass
class GlobalMatrixInspection:
    matrix_path: str
    separator: str
    compression: str | None
    preview_shape: tuple[int, int]
    first_column_name: str | None
    first_column_role: str
    orientation_label: str
    transpose_required_for_anndata: bool
    feature_id_type: str
    n_cells_estimate: int | None
    n_features_estimate: int | None
    barcode_suffix_summary: dict[str, int] = field(default_factory=dict)
    first_features: list[str] = field(default_factory=list)
    first_cells: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    evidence: dict[str, Any] = field(default_factory=dict)


@dataclass
class MetadataTableInspection:
    metadata_path: str
    separator: str
    compression: str | None
    preview_shape: tuple[int, int]
    candidate_cell_id_columns: list[str]
    candidate_sample_columns: list[str]
    first_columns: list[str]
    warnings: list[str] = field(default_factory=list)
    evidence: dict[str, Any] = field(default_factory=dict)


def infer_separator(path: str | Path) -> str:
    name = Path(path).name.lower()

    if name.endswith(".csv") or name.endswith(".csv.gz"):
        return ","

    if name.endswith(".tsv") or name.endswith(".tsv.gz") or name.endswith(".txt") or name.endswith(".txt.gz"):
        return "\t"

    return "\t"


def infer_compression(path: str | Path) -> str | None:
    return "gzip" if Path(path).name.endswith(".gz") else None


def _looks_like_cell_barcode(value: str) -> bool:
    value = str(value).strip().strip('"').strip("'")

    patterns = [
        r"^[ACGT]{8,}-\d+$",
        r"^[ACGT]{8,}\.\d+_\d+$",
        r"^[ACGT]{8,}_\d+$",
        r"^cell[\w.-]*$",
    ]

    return any(re.match(pattern, value, flags=re.IGNORECASE) for pattern in patterns)


def _looks_like_ensembl(value: str) -> bool:
    return str(value).startswith(("ENSG", "ENSMUSG"))


def _infer_feature_id_type(values: list[str]) -> str:
    if not values:
        return "unknown"

    ensembl_count = sum(_looks_like_ensembl(v) for v in values)
    if ensembl_count / len(values) >= 0.5:
        return "ensembl_id"

    return "gene_symbol_or_custom_feature"


def _infer_first_column_role(first_column_name: str | None, values: list[str]) -> str:
    name = str(first_column_name or "").lower()

    if name in {"gene", "genes", "gene_id", "gene_symbol", "features", "unnamed: 0"}:
        return "feature_id"

    if name in {"cell", "cell_id", "cellid", "barcode", "barcodes"}:
        return "cell_id"

    if values and sum(_looks_like_cell_barcode(v) for v in values) / len(values) >= 0.5:
        return "cell_id"

    if values and sum(_looks_like_ensembl(v) for v in values) / len(values) >= 0.5:
        return "feature_id"

    return "unknown"


def inspect_global_matrix(matrix_path: str | Path, preview_rows: int = 20) -> GlobalMatrixInspection:
    matrix_path = Path(matrix_path)
    sep = infer_separator(matrix_path)
    compression = infer_compression(matrix_path)

    preview = pd.read_csv(
        matrix_path,
        sep=sep,
        compression=compression,
        nrows=preview_rows,
    )

    warnings: list[str] = []
    first_column_name = str(preview.columns[0]) if len(preview.columns) else None
    first_column_values = preview.iloc[:, 0].astype(str).head(10).tolist() if preview.shape[1] else []
    column_names = [str(c).strip().strip('"').strip("'") for c in preview.columns]

    first_column_role = _infer_first_column_role(first_column_name, first_column_values)
    feature_id_type = _infer_feature_id_type(first_column_values)

    if first_column_role == "feature_id":
        orientation_label = "gene_by_cell_global_matrix"
        transpose_required = True
        first_features = first_column_values[:5]
        first_cells = column_names[1:6]
        n_cells_estimate = max(preview.shape[1] - 1, 0)
        n_features_estimate = None

        suffix_result = summarize_barcode_suffixes(column_names[1:])
        barcode_suffix_summary = suffix_result.suffix_counts

        if suffix_result.parser_label == "multiple_suffixes_detected":
            warnings.append("Multiple barcode suffixes detected; sample identity may be encoded in cell names.")
        elif suffix_result.parser_label == "no_suffix_detected":
            warnings.append("No barcode suffixes detected in preview columns.")

    elif first_column_role == "cell_id":
        orientation_label = "cell_by_gene_global_matrix"
        transpose_required = False
        first_cells = first_column_values[:5]
        first_features = column_names[1:6]
        n_cells_estimate = None
        n_features_estimate = max(preview.shape[1] - 1, 0)

        suffix_result = summarize_barcode_suffixes(first_column_values)
        barcode_suffix_summary = suffix_result.suffix_counts

    else:
        orientation_label = "unknown_global_matrix_orientation"
        transpose_required = True
        first_features = first_column_values[:5]
        first_cells = column_names[1:6]
        n_cells_estimate = None
        n_features_estimate = None
        barcode_suffix_summary = {}
        warnings.append("Could not confidently infer global matrix orientation.")

    return GlobalMatrixInspection(
        matrix_path=str(matrix_path),
        separator=sep,
        compression=compression,
        preview_shape=tuple(int(x) for x in preview.shape),
        first_column_name=first_column_name,
        first_column_role=first_column_role,
        orientation_label=orientation_label,
        transpose_required_for_anndata=transpose_required,
        feature_id_type=feature_id_type,
        n_cells_estimate=n_cells_estimate,
        n_features_estimate=n_features_estimate,
        barcode_suffix_summary=barcode_suffix_summary,
        first_features=first_features,
        first_cells=first_cells,
        warnings=warnings,
        evidence={
            "preview_rows": preview_rows,
            "columns_preview": column_names[:10],
            "first_column_values_preview": first_column_values[:10],
        },
    )


def inspect_metadata_table(metadata_path: str | Path, preview_rows: int = 50) -> MetadataTableInspection:
    metadata_path = Path(metadata_path)
    sep = infer_separator(metadata_path)
    compression = infer_compression(metadata_path)

    preview = pd.read_csv(
        metadata_path,
        sep=sep,
        compression=compression,
        nrows=preview_rows,
    )

    columns = [str(c) for c in preview.columns]
    lower_to_original = {c.lower(): c for c in columns}

    candidate_cell_keywords = ["cell", "cell_id", "cellid", "barcode", "barcodes"]
    candidate_sample_keywords = ["sample", "sample_id", "gsm", "patient", "patient_id", "library", "batch"]

    candidate_cell_id_columns = [
        col for col in columns
        if any(keyword in col.lower() for keyword in candidate_cell_keywords)
    ]

    candidate_sample_columns = [
        col for col in columns
        if any(keyword in col.lower() for keyword in candidate_sample_keywords)
    ]

    warnings: list[str] = []

    if not candidate_cell_id_columns:
        warnings.append("No obvious cell ID column detected in metadata preview.")

    if not candidate_sample_columns:
        warnings.append("No obvious sample/patient column detected in metadata preview.")

    return MetadataTableInspection(
        metadata_path=str(metadata_path),
        separator=sep,
        compression=compression,
        preview_shape=tuple(int(x) for x in preview.shape),
        candidate_cell_id_columns=candidate_cell_id_columns,
        candidate_sample_columns=candidate_sample_columns,
        first_columns=columns[:20],
        warnings=warnings,
        evidence={
            "preview_rows": preview_rows,
            "columns": columns,
            "lower_to_original": lower_to_original,
        },
    )


def read_global_matrix_with_metadata(
    matrix_path: str | Path,
    metadata_path: str | Path | None = None,
    cell_id_column: str | None = None,
    sample_id_column: str | None = None,
    suffix_to_sample: Mapping[str, str] | None = None,
    keep_metadata_cells_only: bool = False,
):
    try:
        import anndata as ad
    except ImportError as exc:
        raise ImportError("anndata is required to build AnnData objects from global matrices.") from exc

    inspection = inspect_global_matrix(matrix_path)

    sep = infer_separator(matrix_path)
    compression = infer_compression(matrix_path)
    df = pd.read_csv(matrix_path, sep=sep, compression=compression)

    if inspection.first_column_role == "feature_id":
        features = df.iloc[:, 0].astype(str).tolist()
        cells = [str(c).strip().strip('"').strip("'") for c in df.columns[1:]]
        matrix = sparse.csr_matrix(df.iloc[:, 1:].to_numpy()).T
        obs_names = cells
        var_names = features
    elif inspection.first_column_role == "cell_id":
        cells = df.iloc[:, 0].astype(str).tolist()
        features = [str(c) for c in df.columns[1:]]
        matrix = sparse.csr_matrix(df.iloc[:, 1:].to_numpy())
        obs_names = cells
        var_names = features
    else:
        raise ValueError("Could not infer global matrix orientation. Manual review required.")

    adata = ad.AnnData(X=matrix)
    adata.obs_names = obs_names
    adata.var_names = var_names

    if suffix_to_sample:
        suffix_result = summarize_barcode_suffixes(list(adata.obs_names))
        suffix_by_cell = {}
        for cell in adata.obs_names:
            parsed = summarize_barcode_suffixes([cell])
            suffix = next(iter(parsed.suffix_counts), None)
            suffix_by_cell[cell] = suffix_to_sample.get(suffix, suffix) if suffix is not None else None
        adata.obs["sample_id"] = [suffix_by_cell[cell] for cell in adata.obs_names]
        adata.uns["barcode_suffix_counts"] = suffix_result.suffix_counts

    coverage_result = None

    if metadata_path is not None:
        metadata_sep = infer_separator(metadata_path)
        metadata_compression = infer_compression(metadata_path)
        metadata = pd.read_csv(metadata_path, sep=metadata_sep, compression=metadata_compression)

        if cell_id_column is None:
            metadata_inspection = inspect_metadata_table(metadata_path)
            if not metadata_inspection.candidate_cell_id_columns:
                raise ValueError("No cell ID column was provided and none could be inferred from metadata.")
            cell_id_column = metadata_inspection.candidate_cell_id_columns[0]

        metadata[cell_id_column] = metadata[cell_id_column].astype(str)
        metadata = metadata.drop_duplicates(subset=[cell_id_column]).set_index(cell_id_column)

        coverage_result = check_cell_metadata_coverage(
            matrix_cell_ids=list(adata.obs_names),
            metadata_cell_ids=list(metadata.index),
            require_all_matrix_cells_in_metadata=False,
        )

        if keep_metadata_cells_only:
            common_cells = adata.obs_names.intersection(metadata.index)
            adata = adata[common_cells, :].copy()

        common_cells = adata.obs_names.intersection(metadata.index)
        if len(common_cells) > 0:
            for column in metadata.columns:
                adata.obs[column] = pd.NA
                adata.obs.loc[common_cells, column] = metadata.loc[common_cells, column].values

        if sample_id_column and sample_id_column in adata.obs.columns:
            adata.obs["sample_id"] = adata.obs[sample_id_column].astype(str)

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

    return adata, inspection, coverage_result, validation


def print_global_matrix_inspection(inspection: GlobalMatrixInspection) -> None:
    print("=== Global matrix inspection ===")
    print("Matrix path:", inspection.matrix_path)
    print("Separator:", repr(inspection.separator))
    print("Compression:", inspection.compression)
    print("Preview shape:", inspection.preview_shape)
    print("First column name:", inspection.first_column_name)
    print("First column role:", inspection.first_column_role)
    print("Orientation:", inspection.orientation_label)
    print("Transpose required:", inspection.transpose_required_for_anndata)
    print("Feature ID type:", inspection.feature_id_type)
    print("Cell estimate:", inspection.n_cells_estimate)
    print("Feature estimate:", inspection.n_features_estimate)
    print("Barcode suffix summary:", inspection.barcode_suffix_summary)

    print("\nFirst features:")
    for item in inspection.first_features:
        print("-", item)

    print("\nFirst cells:")
    for item in inspection.first_cells:
        print("-", item)

    if inspection.warnings:
        print("\nWarnings:")
        for warning in inspection.warnings:
            print("-", warning)

    print("\nEvidence:")
    for key, value in inspection.evidence.items():
        print(f"- {key}: {value}")


def print_metadata_table_inspection(inspection: MetadataTableInspection) -> None:
    print("=== Metadata table inspection ===")
    print("Metadata path:", inspection.metadata_path)
    print("Separator:", repr(inspection.separator))
    print("Compression:", inspection.compression)
    print("Preview shape:", inspection.preview_shape)
    print("Candidate cell ID columns:", inspection.candidate_cell_id_columns)
    print("Candidate sample columns:", inspection.candidate_sample_columns)
    print("First columns:", inspection.first_columns)

    if inspection.warnings:
        print("\nWarnings:")
        for warning in inspection.warnings:
            print("-", warning)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Inspect a global count matrix and optional metadata table.")
    parser.add_argument("matrix_path", help="Global matrix file.")
    parser.add_argument("--metadata-path", help="Optional metadata table.")

    args = parser.parse_args()

    matrix_inspection = inspect_global_matrix(args.matrix_path)
    print_global_matrix_inspection(matrix_inspection)

    if args.metadata_path:
        print()
        metadata_inspection = inspect_metadata_table(args.metadata_path)
        print_metadata_table_inspection(metadata_inspection)
