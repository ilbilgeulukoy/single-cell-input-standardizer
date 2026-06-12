from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any


@dataclass
class ModalityFilterResult:
    input_path: str
    selected_modality: str
    selected_files: list[str]
    excluded_files: list[str]
    modality_counts: dict[str, int]
    warnings: list[str] = field(default_factory=list)
    evidence: dict[str, Any] = field(default_factory=dict)


def classify_file_modality(path: str | Path) -> str:
    name = Path(path).name.lower()

    if "spatial" in name or name.endswith("spatial.zip") or "tissue_positions" in name or "scalefactors" in name:
        return "spatial_visium"

    if "fragments.tsv" in name or "atac" in name or "peaks" in name:
        return "single_cell_atac"

    if "bulk" in name or "star" in name or "readspergene" in name:
        return "bulk_rna_seq"

    if name.endswith((".mtx.gz", ".mtx", ".h5", ".csv.gz", ".txt.gz", ".tsv.gz")):
        return "single_cell_rna_seq"

    if name.endswith((".tar", ".tar.gz", ".tgz", ".zip")):
        return "archive"

    return "unknown"


def filter_by_modality(
    input_path: str | Path,
    selected_modality: str = "single_cell_rna_seq",
) -> ModalityFilterResult:
    root = Path(input_path)
    files = [root] if root.is_file() else sorted([p for p in root.rglob("*") if p.is_file()])

    selected = []
    excluded = []
    counts: dict[str, int] = {}

    for path in files:
        modality = classify_file_modality(path)
        counts[modality] = counts.get(modality, 0) + 1

        if modality == selected_modality:
            selected.append(str(path))
        else:
            excluded.append(str(path))

    warnings = []
    if len(counts) > 1:
        warnings.append("Multiple modalities or file categories detected. Review excluded files before standardization.")

    if not selected:
        warnings.append(f"No files selected for modality: {selected_modality}")

    return ModalityFilterResult(
        input_path=str(root),
        selected_modality=selected_modality,
        selected_files=selected,
        excluded_files=excluded,
        modality_counts=counts,
        warnings=warnings,
        evidence={
            "n_files": len(files),
            "n_selected": len(selected),
            "n_excluded": len(excluded),
        },
    )


def print_modality_filter_result(result: ModalityFilterResult) -> None:
    print("=== Modality filter result ===")
    print("Input:", result.input_path)
    print("Selected modality:", result.selected_modality)
    print("Selected files:", len(result.selected_files))
    print("Excluded files:", len(result.excluded_files))
    print("Modality counts:", result.modality_counts)

    if result.warnings:
        print("\nWarnings:")
        for warning in result.warnings:
            print("-", warning)


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Filter files by modality/category.")
    parser.add_argument("input_path")
    parser.add_argument("--selected-modality", default="single_cell_rna_seq")

    args = parser.parse_args()
    result = filter_by_modality(args.input_path, selected_modality=args.selected_modality)
    print_modality_filter_result(result)
