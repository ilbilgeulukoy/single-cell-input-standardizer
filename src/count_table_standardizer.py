from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping
import re

import pandas as pd
import numpy as np
from scipy import sparse

from src.sample_id_parser import extract_gsm_id
from src.validation import (
    check_obs_var_alignment,
    check_duplicate_names,
    check_nan_sparse_safe,
    merge_validation_results,
)


@dataclass
class CountTableInspection:
    path: str
    sample_id: str
    separator: str
    compression: str | None
    raw_shape: tuple[int, int]
    orientation_label: str
    transpose_required_for_anndata: bool
    first_column_name: str | None
    first_column_role: str
    feature_id_type: str
    n_cells_estimate: int | None
    n_features_estimate: int | None
    first_features: list[str] = field(default_factory=list)
    first_cells: list[str] = field(default_factory=list)
    duplicated_features: int | None = None
    warnings: list[str] = field(default_factory=list)
    evidence: dict[str, Any] = field(default_factory=dict)


def infer_separator(path: str | Path) -> str:
    name = Path(path).name.lower()

    if name.endswith(".csv") or name.endswith(".csv.gz"):
        return ","

    if name.endswith(".tsv") or name.endswith(".tsv.gz") or name.endswith(".txt") or name.endswith(".txt.gz"):
        return "\t"

    return ","


def infer_compression(path: str | Path) -> str | None:
    return "gzip" if Path(path).name.endswith(".gz") else None


def infer_sample_id_from_path(path: str | Path) -> str:
    path = Path(path)
    gsm = extract_gsm_id(path.name)
    if gsm:
        return gsm

    base = path.name
    for suffix in [".csv.gz", ".txt.gz", ".tsv.gz", ".csv", ".txt", ".tsv"]:
        if base.endswith(suffix):
            base = base[: -len(suffix)]
            break

    return base


def _looks_like_cell_barcode(value: str) -> bool:
    value = str(value)

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

    if name in {"cellid", "cell_id", "barcode", "barcodes"}:
        return "cell_id"

    if name in {"gene", "genes", "gene_id", "gene_symbol", "features", "unnamed: 0"}:
        return "feature_id"

    if values and sum(_looks_like_cell_barcode(v) for v in values) / len(values) >= 0.5:
        return "cell_id"

    if values and (
        sum(_looks_like_ensembl(v) for v in values) / len(values) >= 0.5
        or name.startswith("unnamed")
    ):
        return "feature_id"

    return "unknown"


def inspect_count_table(path: str | Path, preview_rows: int = 20) -> CountTableInspection:
    path = Path(path)
    sep = infer_separator(path)
    compression = infer_compression(path)
    sample_id = infer_sample_id_from_path(path)

    preview = pd.read_csv(
        path,
        sep=sep,
        compression=compression,
        nrows=preview_rows,
    )

    raw_shape_preview = tuple(int(x) for x in preview.shape)
    warnings: list[str] = []

    first_column_name = str(preview.columns[0]) if len(preview.columns) else None
    first_column_values = preview.iloc[:, 0].astype(str).head(10).tolist() if preview.shape[1] else []
    column_names = [str(c) for c in preview.columns]

    first_column_role = _infer_first_column_role(first_column_name, first_column_values)
    feature_id_type = _infer_feature_id_type(first_column_values)

    data_columns = column_names[1:] if preview.shape[1] > 1 else []
    barcode_like_columns = sum(_looks_like_cell_barcode(c) for c in data_columns[:20])

    if first_column_role == "feature_id":
        orientation_label = "gene_by_cell"
        transpose_required = True
        n_features_estimate = None
        n_cells_estimate = max(preview.shape[1] - 1, 0)
        first_features = first_column_values[:5]
        first_cells = data_columns[:5]
    elif first_column_role == "cell_id":
        orientation_label = "cell_by_gene"
        transpose_required = False
        n_cells_estimate = None
        n_features_estimate = max(preview.shape[1] - 1, 0)
        first_features = data_columns[:5]
        first_cells = first_column_values[:5]
    else:
        if barcode_like_columns >= 3:
            orientation_label = "gene_by_cell_likely"
            transpose_required = True
            n_features_estimate = None
            n_cells_estimate = max(preview.shape[1] - 1, 0)
            first_features = first_column_values[:5]
            first_cells = data_columns[:5]
            warnings.append("First column role was ambiguous, but columns look like cell barcodes.")
        else:
            orientation_label = "unknown"
            transpose_required = True
            n_features_estimate = None
            n_cells_estimate = None
            first_features = first_column_values[:5]
            first_cells = data_columns[:5]
            warnings.append("Could not confidently infer count table orientation.")

    duplicated_features = None
    if first_column_role == "feature_id":
        duplicated_features = len(first_column_values) - len(set(first_column_values))

    return CountTableInspection(
        path=str(path),
        sample_id=sample_id,
        separator=sep,
        compression=compression,
        raw_shape=raw_shape_preview,
        orientation_label=orientation_label,
        transpose_required_for_anndata=transpose_required,
        first_column_name=first_column_name,
        first_column_role=first_column_role,
        feature_id_type=feature_id_type,
        n_cells_estimate=n_cells_estimate,
        n_features_estimate=n_features_estimate,
        first_features=first_features,
        first_cells=first_cells,
        duplicated_features=duplicated_features,
        warnings=warnings,
        evidence={
            "preview_rows": preview_rows,
            "preview_shape": raw_shape_preview,
            "barcode_like_columns_in_first_20": barcode_like_columns,
            "columns_preview": column_names[:10],
            "first_column_values_preview": first_column_values[:10],
        },
    )


def inspect_many_count_tables(paths: list[str | Path]) -> list[CountTableInspection]:
    return [inspect_count_table(path) for path in paths]


def discover_count_tables(input_path: str | Path) -> list[Path]:
    root = Path(input_path)

    if root.is_file():
        files = [root]
    else:
        files = sorted([p for p in root.rglob("*") if p.is_file()])

    allowed_suffixes = (".csv.gz", ".csv", ".txt.gz", ".txt", ".tsv.gz", ".tsv")

    return [
        p for p in files
        if p.name.endswith(allowed_suffixes)
        and not p.name.endswith(("features.tsv.gz", "genes.tsv.gz", "barcodes.tsv.gz"))
    ]


def summarize_count_table_inspections(inspections: list[CountTableInspection]) -> pd.DataFrame:
    rows = []

    for item in inspections:
        rows.append(
            {
                "path": item.path,
                "sample_id": item.sample_id,
                "separator": item.separator,
                "compression": item.compression,
                "raw_shape_preview": item.raw_shape,
                "orientation_label": item.orientation_label,
                "transpose_required_for_anndata": item.transpose_required_for_anndata,
                "first_column_name": item.first_column_name,
                "first_column_role": item.first_column_role,
                "feature_id_type": item.feature_id_type,
                "n_cells_estimate": item.n_cells_estimate,
                "n_features_estimate": item.n_features_estimate,
                "duplicated_features_in_preview": item.duplicated_features,
                "warnings": " | ".join(item.warnings),
            }
        )

    return pd.DataFrame(rows)


def read_count_table_with_metadata(
    path: str | Path,
    sample_id: str | None = None,
    metadata: Mapping[str, Any] | None = None,
    prefix_obs_names: bool = True,
):
    try:
        import anndata as ad
    except ImportError as exc:
        raise ImportError("anndata is required to build AnnData objects from count tables.") from exc

    path = Path(path)
    inspection = inspect_count_table(path)

    sep = infer_separator(path)
    compression = infer_compression(path)

    df = pd.read_csv(path, sep=sep, compression=compression)
    resolved_sample_id = sample_id or inspection.sample_id

    if inspection.first_column_role == "feature_id":
        features = df.iloc[:, 0].astype(str).tolist()
        cells = [str(c) for c in df.columns[1:]]
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
        raise ValueError("Could not infer count table orientation. Manual review required.")

    if prefix_obs_names:
        obs_names = [f"{resolved_sample_id}_{name}" for name in obs_names]

    adata = ad.AnnData(X=matrix)
    adata.obs_names = obs_names
    adata.var_names = var_names
    adata.obs["sample_id"] = resolved_sample_id

    if metadata:
        for key, value in metadata.items():
            adata.obs[key] = value

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

    return adata, inspection, validation


def print_count_table_inspection(inspection: CountTableInspection) -> None:
    print("=== Count table inspection ===")
    print("Path:", inspection.path)
    print("Sample ID:", inspection.sample_id)
    print("Separator:", repr(inspection.separator))
    print("Compression:", inspection.compression)
    print("Raw preview shape:", inspection.raw_shape)
    print("Orientation:", inspection.orientation_label)
    print("Transpose required:", inspection.transpose_required_for_anndata)
    print("First column name:", inspection.first_column_name)
    print("First column role:", inspection.first_column_role)
    print("Feature ID type:", inspection.feature_id_type)
    print("Cells estimate:", inspection.n_cells_estimate)
    print("Features estimate:", inspection.n_features_estimate)
    print("Duplicated features in preview:", inspection.duplicated_features)

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


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Inspect compressed CSV/TXT count tables.")
    parser.add_argument("input_path", help="Input file or directory.")
    parser.add_argument("--summary-csv", help="Optional CSV path for inspection summary.")

    args = parser.parse_args()

    tables = discover_count_tables(args.input_path)
    inspections = inspect_many_count_tables(tables)

    print("Detected count tables:", len(tables))
    print()

    for inspection in inspections:
        print_count_table_inspection(inspection)
        print()

    if args.summary_csv:
        df = summarize_count_table_inspections(inspections)
        df.to_csv(args.summary_csv, index=False)
        print("Wrote:", args.summary_csv)
