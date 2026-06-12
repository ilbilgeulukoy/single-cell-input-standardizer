from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Mapping

import h5py
import pandas as pd

from src.sample_id_parser import extract_gsm_id
from src.validation import (
    check_obs_var_alignment,
    check_duplicate_names,
    check_nan_sparse_safe,
    merge_validation_results,
)


@dataclass
class TenxH5Inspection:
    path: str
    layout_label: str
    genome_group: str | None
    n_cells: int | None
    n_genes: int | None
    gene_id_source: str | None
    gene_name_source: str | None
    barcode_source: str | None
    first_barcodes: list[str] = field(default_factory=list)
    first_gene_ids: list[str] = field(default_factory=list)
    first_gene_names: list[str] = field(default_factory=list)
    duplicated_gene_names: int | None = None
    warnings: list[str] = field(default_factory=list)
    evidence: dict[str, Any] = field(default_factory=dict)


def _decode_array(values: Any, limit: int | None = None) -> list[str]:
    values = values[:limit] if limit is not None else values[:]
    decoded = []

    for value in values:
        if isinstance(value, bytes):
            decoded.append(value.decode("utf-8"))
        else:
            decoded.append(str(value))

    return decoded


def _count_duplicates(values: list[str]) -> int:
    return len(values) - len(set(values))


def detect_10x_h5_layout(path: str | Path) -> TenxH5Inspection:
    path = Path(path)
    warnings: list[str] = []

    with h5py.File(path, "r") as handle:
        top_level_keys = list(handle.keys())

        if "matrix" in handle:
            matrix = handle["matrix"]
            matrix_subkeys = list(matrix.keys())

            shape = matrix["shape"][:] if "shape" in matrix else None
            n_genes = int(shape[0]) if shape is not None else None
            n_cells = int(shape[1]) if shape is not None else None

            first_barcodes = (
                _decode_array(matrix["barcodes"], limit=5)
                if "barcodes" in matrix
                else []
            )

            gene_id_source = None
            gene_name_source = None
            first_gene_ids = []
            first_gene_names = []
            duplicated_gene_names = None
            feature_keys = []

            if "features" in matrix:
                features = matrix["features"]
                feature_keys = list(features.keys())

                if "id" in features:
                    gene_id_source = "matrix/features/id"
                    first_gene_ids = _decode_array(features["id"], limit=5)

                if "name" in features:
                    gene_name_source = "matrix/features/name"
                    first_gene_names = _decode_array(features["name"], limit=5)
                    all_gene_names = _decode_array(features["name"])
                    duplicated_gene_names = _count_duplicates(all_gene_names)
            else:
                warnings.append("No matrix/features group found.")

            return TenxH5Inspection(
                path=str(path),
                layout_label="modern_10x_h5_v3_matrix_group",
                genome_group=None,
                n_cells=n_cells,
                n_genes=n_genes,
                gene_id_source=gene_id_source,
                gene_name_source=gene_name_source,
                barcode_source="matrix/barcodes" if "barcodes" in matrix else None,
                first_barcodes=first_barcodes,
                first_gene_ids=first_gene_ids,
                first_gene_names=first_gene_names,
                duplicated_gene_names=duplicated_gene_names,
                warnings=warnings,
                evidence={
                    "top_level_keys": top_level_keys,
                    "matrix_subkeys": matrix_subkeys,
                    "feature_keys": feature_keys,
                    "raw_shape_genes_x_cells": tuple(shape) if shape is not None else None,
                },
            )

        genome_groups = [
            key for key in top_level_keys
            if isinstance(handle[key], h5py.Group)
        ]

        if genome_groups:
            genome_group = genome_groups[0]
            matrix = handle[genome_group]
            subkeys = list(matrix.keys())

            shape = matrix["shape"][:] if "shape" in matrix else None
            n_genes = int(shape[0]) if shape is not None else None
            n_cells = int(shape[1]) if shape is not None else None

            first_barcodes = (
                _decode_array(matrix["barcodes"], limit=5)
                if "barcodes" in matrix
                else []
            )
            first_gene_ids = (
                _decode_array(matrix["genes"], limit=5)
                if "genes" in matrix
                else []
            )
            first_gene_names = (
                _decode_array(matrix["gene_names"], limit=5)
                if "gene_names" in matrix
                else []
            )

            duplicated_gene_names = None
            if "gene_names" in matrix:
                all_gene_names = _decode_array(matrix["gene_names"])
                duplicated_gene_names = _count_duplicates(all_gene_names)

            return TenxH5Inspection(
                path=str(path),
                layout_label="old_10x_h5_v2_genome_group",
                genome_group=genome_group,
                n_cells=n_cells,
                n_genes=n_genes,
                gene_id_source=f"{genome_group}/genes" if "genes" in matrix else None,
                gene_name_source=f"{genome_group}/gene_names" if "gene_names" in matrix else None,
                barcode_source=f"{genome_group}/barcodes" if "barcodes" in matrix else None,
                first_barcodes=first_barcodes,
                first_gene_ids=first_gene_ids,
                first_gene_names=first_gene_names,
                duplicated_gene_names=duplicated_gene_names,
                warnings=warnings,
                evidence={
                    "top_level_keys": top_level_keys,
                    "genome_group_subkeys": subkeys,
                    "raw_shape_genes_x_cells": tuple(shape) if shape is not None else None,
                },
            )

    return TenxH5Inspection(
        path=str(path),
        layout_label="unknown_h5_layout",
        genome_group=None,
        n_cells=None,
        n_genes=None,
        gene_id_source=None,
        gene_name_source=None,
        barcode_source=None,
        warnings=["Could not detect a supported 10x h5 layout."],
        evidence={"top_level_keys": top_level_keys},
    )


def inspect_many_10x_h5(paths: list[str | Path]) -> list[TenxH5Inspection]:
    return [detect_10x_h5_layout(path) for path in paths]


def summarize_inspections(inspections: list[TenxH5Inspection]) -> pd.DataFrame:
    rows = []

    for item in inspections:
        rows.append(
            {
                "path": item.path,
                "layout_label": item.layout_label,
                "genome_group": item.genome_group,
                "n_cells": item.n_cells,
                "n_genes": item.n_genes,
                "duplicated_gene_names": item.duplicated_gene_names,
                "gene_id_source": item.gene_id_source,
                "gene_name_source": item.gene_name_source,
                "barcode_source": item.barcode_source,
                "warnings": " | ".join(item.warnings),
            }
        )

    return pd.DataFrame(rows)


def read_10x_h5_with_metadata(
    path: str | Path,
    sample_id: str | None = None,
    metadata: Mapping[str, Any] | None = None,
    prefix_obs_names: bool = True,
):
    try:
        import scanpy as sc
    except ImportError as exc:
        raise ImportError(
            "scanpy is required to read 10x h5 into AnnData. "
            "Use detect_10x_h5_layout for h5py-only inspection."
        ) from exc

    path = Path(path)
    inspection = detect_10x_h5_layout(path)
    adata = sc.read_10x_h5(str(path))

    inferred_sample_id = sample_id or extract_gsm_id(path.name) or path.stem

    if prefix_obs_names:
        adata.obs_names = [f"{inferred_sample_id}_{barcode}" for barcode in adata.obs_names]

    adata.obs["sample_id"] = inferred_sample_id

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


def print_h5_inspection(inspection: TenxH5Inspection) -> None:
    print("=== 10x h5 inspection ===")
    print("Path:", inspection.path)
    print("Layout:", inspection.layout_label)
    print("Genome group:", inspection.genome_group)
    print("Cells:", inspection.n_cells)
    print("Genes:", inspection.n_genes)
    print("Gene ID source:", inspection.gene_id_source)
    print("Gene name source:", inspection.gene_name_source)
    print("Barcode source:", inspection.barcode_source)
    print("Duplicated gene names:", inspection.duplicated_gene_names)

    print("\nFirst barcodes:")
    for item in inspection.first_barcodes:
        print("-", item)

    print("\nFirst gene IDs:")
    for item in inspection.first_gene_ids:
        print("-", item)

    print("\nFirst gene names:")
    for item in inspection.first_gene_names:
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

    parser = argparse.ArgumentParser(description="Inspect 10x h5 files.")
    parser.add_argument("h5_files", nargs="+", help="One or more 10x h5 files.")
    parser.add_argument("--summary-csv", help="Optional CSV path for inspection summary.")

    args = parser.parse_args()

    inspections = inspect_many_10x_h5(args.h5_files)

    for inspection in inspections:
        print_h5_inspection(inspection)
        print()

    if args.summary_csv:
        df = summarize_inspections(inspections)
        df.to_csv(args.summary_csv, index=False)
        print("Wrote:", args.summary_csv)
