from __future__ import annotations

import argparse
import json
from dataclasses import asdict, is_dataclass
from pathlib import Path
from typing import Any

import pandas as pd

from src.file_classifier import classify_input, print_classification
from src.data_size_policy import evaluate_data_size, print_data_size_decision
from src.modality_filter import filter_by_modality
from src.tenx_h5_standardizer import inspect_many_10x_h5, summarize_inspections as summarize_h5_inspections
from src.tenx_mtx_standardizer import (
    discover_10x_mtx_triplets,
    inspect_many_10x_mtx_triplets,
    summarize_mtx_inspections,
)
from src.count_table_standardizer import (
    discover_count_tables,
    inspect_many_count_tables,
    summarize_count_table_inspections,
)
from src.global_matrix_metadata_standardizer import (
    inspect_global_matrix,
    inspect_metadata_table,
)


def make_json_safe(value: Any) -> Any:
    if is_dataclass(value):
        return make_json_safe(asdict(value))

    if isinstance(value, dict):
        return {str(k): make_json_safe(v) for k, v in value.items()}

    if isinstance(value, list):
        return [make_json_safe(v) for v in value]

    if isinstance(value, tuple):
        return [make_json_safe(v) for v in value]

    if isinstance(value, Path):
        return str(value)

    if hasattr(value, "item"):
        try:
            return value.item()
        except Exception:
            return str(value)

    return value


def write_json(path: Path, data: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w", encoding="utf-8") as handle:
        json.dump(make_json_safe(data), handle, indent=2)


def write_text(path: Path, text: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding="utf-8")


def write_dataframe(path: Path, df: pd.DataFrame) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def get_any(obj: Any, *names: str, default: Any = None) -> Any:
    for name in names:
        if hasattr(obj, name):
            return getattr(obj, name)

    if isinstance(obj, dict):
        for name in names:
            if name in obj:
                return obj[name]

    return default


def collect_h5_files(input_path: Path) -> list[Path]:
    if input_path.is_file() and input_path.name.endswith(".h5"):
        return [input_path]

    if input_path.is_dir():
        return sorted([p for p in input_path.rglob("*.h5") if p.is_file()])

    return []


def build_manual_review_report(
    input_path: Path,
    classification: Any,
    size_decision: Any,
    modality_result: Any,
    format_summary: dict[str, Any],
) -> str:
    detected_format = get_any(
        classification,
        "detected_format",
        "format_label",
        default="unknown",
    )
    recommended_module = get_any(
        classification,
        "recommended_module",
        default="unknown",
    )
    size_policy = get_any(
        size_decision,
        "decision_label",
        "policy_label",
        default="unknown",
    )
    total_size_bytes = get_any(
        size_decision,
        "total_size_bytes",
        default="not_available",
    )
    total_size_readable = get_any(
        size_decision,
        "total_size_readable",
        "total_size_mb",
        default="not_available",
    )

    lines = []

    lines.append("# Manual review report")
    lines.append("")
    lines.append(f"Input path: `{input_path}`")
    lines.append("")

    lines.append("## File classification")
    lines.append("")
    lines.append(f"- Detected format: `{detected_format}`")
    lines.append(f"- Recommended module: `{recommended_module}`")
    lines.append("")

    lines.append("## Data size policy")
    lines.append("")
    lines.append(f"- Decision: `{size_policy}`")
    lines.append(f"- Total size bytes: `{total_size_bytes}`")
    lines.append(f"- Total size readable: `{total_size_readable}`")
    lines.append("")

    if getattr(size_decision, "warnings", []):
        lines.append("Warnings:")
        for warning in size_decision.warnings:
            lines.append(f"- {warning}")
        lines.append("")

    lines.append("## Modality filter")
    lines.append("")
    lines.append(f"- Selected modality: `{modality_result.selected_modality}`")
    lines.append(f"- Selected files: `{len(modality_result.selected_files)}`")
    lines.append(f"- Excluded files: `{len(modality_result.excluded_files)}`")
    lines.append(f"- Modality counts: `{modality_result.modality_counts}`")
    lines.append("")

    if modality_result.warnings:
        lines.append("Warnings:")
        for warning in modality_result.warnings:
            lines.append(f"- {warning}")
        lines.append("")

    lines.append("## Format-specific inspection")
    lines.append("")

    if not format_summary:
        lines.append("- No format-specific inspection was run.")
    else:
        for key, value in format_summary.items():
            lines.append(f"- {key}: `{value}`")
    lines.append("")

    lines.append("## Manual review checklist")
    lines.append("")
    lines.append("- Confirm that detected sample IDs match the biological metadata.")
    lines.append("- Confirm whether excluded files are intentionally excluded.")
    lines.append("- Confirm whether duplicate gene symbols should be aggregated or preserved.")
    lines.append("- Confirm whether metadata is complete, partial, or authoritative subset only.")
    lines.append("- Confirm whether the selected standardizer matches the actual input format.")
    lines.append("")

    lines.append("## Scope note")
    lines.append("")
    lines.append(
        "This toolkit does not claim to fully automate biological metadata interpretation. "
        "It provides structured inspection, candidate identifier extraction, consistency checks, "
        "and manual review outputs before AnnData standardization."
    )
    lines.append("")

    return "\n".join(lines)


def run_inspection(input_path: str | Path, output_dir: str | Path, metadata_path: str | Path | None = None) -> None:
    input_path = Path(input_path)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    print("=== Unified single-cell input inspection ===")
    print("Input:", input_path)
    print("Output:", output_dir)
    print()

    classification = classify_input(input_path)
    size_decision = evaluate_data_size(input_path)
    modality_result = filter_by_modality(input_path)

    print_classification(classification)
    print()
    print_data_size_decision(size_decision)
    print()

    write_json(output_dir / "file_classification.json", classification)
    write_json(output_dir / "data_size_policy.json", size_decision)
    write_json(output_dir / "modality_filter.json", modality_result)

    format_summary: dict[str, Any] = {}

    detected_format = get_any(classification, "detected_format", "format_label", default="unknown")

    if detected_format == "per_sample_10x_h5_files":
        h5_files = collect_h5_files(input_path)
        inspections = inspect_many_10x_h5(h5_files)
        df = summarize_h5_inspections(inspections)

        write_dataframe(output_dir / "format_inspection_summary.csv", df)
        write_json(output_dir / "format_inspection_details.json", inspections)

        format_summary = {
            "format": "10x h5",
            "n_h5_files": len(h5_files),
            "summary_csv": "format_inspection_summary.csv",
        }

    elif detected_format == "per_sample_10x_matrix_market_triplets":
        triplets = discover_10x_mtx_triplets(input_path)
        inspections = inspect_many_10x_mtx_triplets(triplets)
        df = summarize_mtx_inspections(inspections)

        write_dataframe(output_dir / "format_inspection_summary.csv", df)
        write_json(output_dir / "format_inspection_details.json", inspections)

        format_summary = {
            "format": "10x Matrix Market",
            "n_triplets": len(triplets),
            "summary_csv": "format_inspection_summary.csv",
        }

    elif detected_format in {"compressed_txt_count_tables", "multiple_compressed_count_tables"}:
        tables = discover_count_tables(input_path)
        inspections = inspect_many_count_tables(tables)
        df = summarize_count_table_inspections(inspections)

        write_dataframe(output_dir / "format_inspection_summary.csv", df)
        write_json(output_dir / "format_inspection_details.json", inspections)

        format_summary = {
            "format": "compressed count tables",
            "n_tables": len(tables),
            "summary_csv": "format_inspection_summary.csv",
        }

    elif detected_format == "single_global_compressed_csv_count_matrix":
        matrix_inspection = inspect_global_matrix(input_path)
        write_json(output_dir / "format_inspection_details.json", matrix_inspection)

        df = pd.DataFrame(
            [
                {
                    "matrix_path": matrix_inspection.matrix_path,
                    "orientation_label": matrix_inspection.orientation_label,
                    "transpose_required_for_anndata": matrix_inspection.transpose_required_for_anndata,
                    "first_column_role": matrix_inspection.first_column_role,
                    "feature_id_type": matrix_inspection.feature_id_type,
                    "n_cells_estimate": matrix_inspection.n_cells_estimate,
                    "n_features_estimate": matrix_inspection.n_features_estimate,
                    "warnings": " | ".join(matrix_inspection.warnings),
                }
            ]
        )
        write_dataframe(output_dir / "format_inspection_summary.csv", df)

        format_summary = {
            "format": "single global matrix",
            "summary_csv": "format_inspection_summary.csv",
        }

    elif detected_format == "archive_or_nested_archive":
        format_summary = {
            "format": "archive or nested archive",
            "note": "Run nested_archive_extractor.py before format-specific inspection.",
        }

    else:
        format_summary = {
            "format": detected_format,
            "note": "No dedicated format-specific inspection implemented for this label yet.",
        }

    if metadata_path is not None:
        metadata_inspection = inspect_metadata_table(metadata_path)
        write_json(output_dir / "metadata_table_inspection.json", metadata_inspection)
        format_summary["metadata_table_inspection"] = "metadata_table_inspection.json"

    report = build_manual_review_report(
        input_path=input_path,
        classification=classification,
        size_decision=size_decision,
        modality_result=modality_result,
        format_summary=format_summary,
    )

    write_text(output_dir / "manual_review_report.md", report)

    manifest = {
        "input_path": str(input_path),
        "output_dir": str(output_dir),
        "files_written": sorted([str(p.relative_to(output_dir)) for p in output_dir.rglob("*") if p.is_file()]),
    }
    write_json(output_dir / "inspection_manifest.json", manifest)

    print("Inspection complete.")
    print("Files written:")
    for item in manifest["files_written"]:
        print("-", item)


def main() -> None:
    parser = argparse.ArgumentParser(description="Unified single-cell input inspection CLI.")
    parser.add_argument("input_path", help="Input file or directory to inspect.")
    parser.add_argument("--output", required=True, help="Output report directory.")
    parser.add_argument("--metadata-path", help="Optional metadata table to inspect.")

    args = parser.parse_args()

    run_inspection(
        input_path=args.input_path,
        output_dir=args.output,
        metadata_path=args.metadata_path,
    )


if __name__ == "__main__":
    main()
