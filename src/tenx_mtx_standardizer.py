from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping
import gzip

import pandas as pd
from scipy.io import mmread

from src.sample_id_parser import extract_gsm_id
from src.validation import (
    check_obs_var_alignment,
    check_duplicate_names,
    check_nan_sparse_safe,
    merge_validation_results,
)


@dataclass
class TenxMtxTriplet:
    sample_id: str
    matrix_path: str
    features_path: str
    barcodes_path: str
    feature_table_type: str
    warnings: list[str] = field(default_factory=list)


@dataclass
class TenxMtxInspection:
    sample_id: str
    matrix_path: str
    features_path: str
    barcodes_path: str
    raw_matrix_shape: tuple[int, int] | None
    expected_anndata_shape: tuple[int, int] | None
    n_features: int | None
    n_barcodes: int | None
    feature_table_type: str
    gene_id_column: int | None
    gene_symbol_column: int | None
    feature_type_column: int | None
    first_barcodes: list[str] = field(default_factory=list)
    first_gene_ids: list[str] = field(default_factory=list)
    first_gene_symbols: list[str] = field(default_factory=list)
    duplicated_gene_symbols: int | None = None
    transpose_required_for_anndata: bool = True
    warnings: list[str] = field(default_factory=list)
    evidence: dict[str, Any] = field(default_factory=dict)


def _read_first_lines(path: str | Path, n: int = 5) -> list[str]:
    path = Path(path)

    opener = gzip.open if path.name.endswith(".gz") else open
    lines = []

    with opener(path, "rt", encoding="utf-8") as handle:
        for _ in range(n):
            line = handle.readline()
            if not line:
                break
            lines.append(line.rstrip("\n"))

    return lines


def _read_tsv_no_header(path: str | Path) -> pd.DataFrame:
    return pd.read_csv(path, sep="\t", header=None, compression="infer")


def _count_duplicates(values: list[str]) -> int:
    return len(values) - len(set(values))


def _role_from_name(path: Path) -> str | None:
    name = path.name.lower()

    if name.endswith(".mtx.gz") or name.endswith(".mtx"):
        if "matrix" in name:
            return "matrix"

    if name.endswith(".tsv.gz") or name.endswith(".tsv"):
        if "barcodes" in name:
            return "barcodes"
        if "features" in name or "genes" in name:
            return "features"

    return None


def _sample_key_from_name(path: Path) -> str:
    gsm = extract_gsm_id(path.name)
    if gsm:
        return gsm

    name = path.name
    for token in ["_matrix", ".matrix", "_features", ".features", "_genes", ".genes", "_barcodes", ".barcodes"]:
        if token in name:
            return name.split(token)[0]

    return path.stem


def discover_10x_mtx_triplets(input_path: str | Path) -> list[TenxMtxTriplet]:
    root = Path(input_path)

    if root.is_file():
        files = [root]
    else:
        files = sorted([p for p in root.rglob("*") if p.is_file()])

    grouped: dict[str, dict[str, list[Path]]] = {}

    for path in files:
        role = _role_from_name(path)
        if role is None:
            continue

        sample_key = _sample_key_from_name(path)
        grouped.setdefault(sample_key, {"matrix": [], "features": [], "barcodes": []})
        grouped[sample_key][role].append(path)

    triplets: list[TenxMtxTriplet] = []

    for sample_id, roles in sorted(grouped.items()):
        warnings = []

        if len(roles["matrix"]) != 1:
            warnings.append(f"Expected 1 matrix file, found {len(roles['matrix'])}.")
        if len(roles["features"]) != 1:
            warnings.append(f"Expected 1 features/genes file, found {len(roles['features'])}.")
        if len(roles["barcodes"]) != 1:
            warnings.append(f"Expected 1 barcodes file, found {len(roles['barcodes'])}.")

        if warnings:
            continue

        features_path = roles["features"][0]
        feature_table_type = "features.tsv" if "features" in features_path.name.lower() else "genes.tsv"

        triplets.append(
            TenxMtxTriplet(
                sample_id=sample_id,
                matrix_path=str(roles["matrix"][0]),
                features_path=str(features_path),
                barcodes_path=str(roles["barcodes"][0]),
                feature_table_type=feature_table_type,
                warnings=warnings,
            )
        )

    return triplets


def inspect_10x_mtx_triplet(triplet: TenxMtxTriplet) -> TenxMtxInspection:
    warnings: list[str] = []

    matrix = mmread(triplet.matrix_path)
    raw_shape = tuple(int(x) for x in matrix.shape)

    features = _read_tsv_no_header(triplet.features_path)
    barcodes = _read_tsv_no_header(triplet.barcodes_path)

    n_features = int(features.shape[0])
    n_barcodes = int(barcodes.shape[0])

    if raw_shape[0] == n_features and raw_shape[1] == n_barcodes:
        transpose_required = True
        expected_shape = (n_barcodes, n_features)
    elif raw_shape[0] == n_barcodes and raw_shape[1] == n_features:
        transpose_required = False
        expected_shape = raw_shape
        warnings.append("Matrix appears already cell x gene; transpose may not be required.")
    else:
        transpose_required = True
        expected_shape = None
        warnings.append("Matrix shape does not match feature/barcode counts.")

    if features.shape[1] >= 3:
        gene_id_column = 0
        gene_symbol_column = 1
        feature_type_column = 2
    elif features.shape[1] == 2:
        gene_id_column = 0
        gene_symbol_column = 1
        feature_type_column = None
    elif features.shape[1] == 1:
        gene_id_column = None
        gene_symbol_column = 0
        feature_type_column = None
        warnings.append("Feature table has one column; treating it as gene symbol.")
    else:
        gene_id_column = None
        gene_symbol_column = None
        feature_type_column = None
        warnings.append("Feature table is empty or malformed.")

    first_barcodes = barcodes.iloc[:5, 0].astype(str).tolist() if n_barcodes else []

    first_gene_ids: list[str] = []
    first_gene_symbols: list[str] = []
    duplicated_gene_symbols = None

    if gene_id_column is not None:
        first_gene_ids = features.iloc[:5, gene_id_column].astype(str).tolist()

    if gene_symbol_column is not None:
        first_gene_symbols = features.iloc[:5, gene_symbol_column].astype(str).tolist()
        all_gene_symbols = features.iloc[:, gene_symbol_column].astype(str).tolist()
        duplicated_gene_symbols = _count_duplicates(all_gene_symbols)

    return TenxMtxInspection(
        sample_id=triplet.sample_id,
        matrix_path=triplet.matrix_path,
        features_path=triplet.features_path,
        barcodes_path=triplet.barcodes_path,
        raw_matrix_shape=raw_shape,
        expected_anndata_shape=expected_shape,
        n_features=n_features,
        n_barcodes=n_barcodes,
        feature_table_type=triplet.feature_table_type,
        gene_id_column=gene_id_column,
        gene_symbol_column=gene_symbol_column,
        feature_type_column=feature_type_column,
        first_barcodes=first_barcodes,
        first_gene_ids=first_gene_ids,
        first_gene_symbols=first_gene_symbols,
        duplicated_gene_symbols=duplicated_gene_symbols,
        transpose_required_for_anndata=transpose_required,
        warnings=warnings,
        evidence={
            "feature_table_columns": int(features.shape[1]),
            "matrix_nnz": int(matrix.nnz) if hasattr(matrix, "nnz") else None,
            "first_feature_lines": _read_first_lines(triplet.features_path, n=3),
            "first_barcode_lines": _read_first_lines(triplet.barcodes_path, n=3),
        },
    )


def inspect_many_10x_mtx_triplets(triplets: list[TenxMtxTriplet]) -> list[TenxMtxInspection]:
    return [inspect_10x_mtx_triplet(triplet) for triplet in triplets]


def summarize_mtx_inspections(inspections: list[TenxMtxInspection]) -> pd.DataFrame:
    rows = []

    for item in inspections:
        rows.append(
            {
                "sample_id": item.sample_id,
                "raw_matrix_shape": item.raw_matrix_shape,
                "expected_anndata_shape": item.expected_anndata_shape,
                "n_features": item.n_features,
                "n_barcodes": item.n_barcodes,
                "feature_table_type": item.feature_table_type,
                "gene_symbol_column": item.gene_symbol_column,
                "duplicated_gene_symbols": item.duplicated_gene_symbols,
                "transpose_required_for_anndata": item.transpose_required_for_anndata,
                "warnings": " | ".join(item.warnings),
            }
        )

    return pd.DataFrame(rows)


def read_10x_mtx_triplet_with_metadata(
    triplet: TenxMtxTriplet,
    metadata: Mapping[str, Any] | None = None,
    prefix_obs_names: bool = True,
):
    try:
        import anndata as ad
    except ImportError as exc:
        raise ImportError("anndata is required to build AnnData objects from 10x MTX triplets.") from exc

    inspection = inspect_10x_mtx_triplet(triplet)

    matrix = mmread(triplet.matrix_path).tocsr()
    if inspection.transpose_required_for_anndata:
        matrix = matrix.T.tocsr()

    features = _read_tsv_no_header(triplet.features_path)
    barcodes = _read_tsv_no_header(triplet.barcodes_path)

    gene_symbol_column = inspection.gene_symbol_column
    if gene_symbol_column is None:
        raise ValueError("Could not identify gene symbol column.")

    obs_names = barcodes.iloc[:, 0].astype(str).tolist()
    if prefix_obs_names:
        obs_names = [f"{triplet.sample_id}_{barcode}" for barcode in obs_names]

    var_names = features.iloc[:, gene_symbol_column].astype(str).tolist()

    adata = ad.AnnData(X=matrix)
    adata.obs_names = obs_names
    adata.var_names = var_names

    adata.obs["sample_id"] = triplet.sample_id

    if metadata:
        for key, value in metadata.items():
            adata.obs[key] = value

    if inspection.gene_id_column is not None:
        adata.var["gene_id"] = features.iloc[:, inspection.gene_id_column].astype(str).values

    if inspection.feature_type_column is not None:
        adata.var["feature_type"] = features.iloc[:, inspection.feature_type_column].astype(str).values

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


def print_mtx_inspection(inspection: TenxMtxInspection) -> None:
    print("=== 10x MTX triplet inspection ===")
    print("Sample ID:", inspection.sample_id)
    print("Matrix:", inspection.matrix_path)
    print("Features:", inspection.features_path)
    print("Barcodes:", inspection.barcodes_path)
    print("Raw matrix shape:", inspection.raw_matrix_shape)
    print("Expected AnnData shape:", inspection.expected_anndata_shape)
    print("Features:", inspection.n_features)
    print("Barcodes:", inspection.n_barcodes)
    print("Feature table type:", inspection.feature_table_type)
    print("Gene ID column:", inspection.gene_id_column)
    print("Gene symbol column:", inspection.gene_symbol_column)
    print("Feature type column:", inspection.feature_type_column)
    print("Transpose required:", inspection.transpose_required_for_anndata)
    print("Duplicated gene symbols:", inspection.duplicated_gene_symbols)

    print("\nFirst barcodes:")
    for item in inspection.first_barcodes:
        print("-", item)

    print("\nFirst gene IDs:")
    for item in inspection.first_gene_ids:
        print("-", item)

    print("\nFirst gene symbols:")
    for item in inspection.first_gene_symbols:
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

    parser = argparse.ArgumentParser(description="Inspect 10x Matrix Market triplets.")
    parser.add_argument("input_path", help="Input directory containing 10x Matrix Market files.")
    parser.add_argument("--summary-csv", help="Optional CSV path for inspection summary.")

    args = parser.parse_args()

    triplets = discover_10x_mtx_triplets(args.input_path)
    inspections = inspect_many_10x_mtx_triplets(triplets)

    print("Detected complete 10x MTX triplets:", len(triplets))
    print()

    for inspection in inspections:
        print_mtx_inspection(inspection)
        print()

    if args.summary_csv:
        df = summarize_mtx_inspections(inspections)
        df.to_csv(args.summary_csv, index=False)
        print("Wrote:", args.summary_csv)
