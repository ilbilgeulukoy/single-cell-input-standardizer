from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from src.global_matrix_metadata_standardizer import inspect_metadata_table
from src.standard_schema import STANDARD_OBS_COLUMNS


REVIEW_REQUIRED_FIELDS = {
    "patient_id",
    "cancer_type",
    "tumor_site",
    "metastasis_site",
    "tumor_treatment",
    "cancer_site_origin",
    "tumour_grade",
    "tumour_stage",
    "histological_subtype",
    "patient_ethnicity",
}


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def infer_separator(path: Path) -> str:
    if path.name.endswith((".csv", ".csv.gz")):
        return ","
    return "\t"


def infer_compression(path: Path) -> str | None:
    return "gzip" if path.name.endswith(".gz") else None


def load_metadata_preview(metadata_path: Path) -> pd.DataFrame:
    return pd.read_csv(
        metadata_path,
        sep=infer_separator(metadata_path),
        compression=infer_compression(metadata_path),
    )


def summarize_metadata_samples(
    metadata: pd.DataFrame,
    sample_id_column: str,
) -> pd.DataFrame:
    if sample_id_column not in metadata.columns:
        raise ValueError(f"sample_id_column '{sample_id_column}' not found in metadata.")

    aggregation: dict[str, Any] = {
        "n_cells": (sample_id_column, "size"),
    }

    for column in ["Condition", "Location", "MSI_Status", "condition", "location", "msi_status"]:
        if column in metadata.columns:
            aggregation[column] = (column, lambda values: sorted(set(map(str, values))))

    return metadata.groupby(sample_id_column).agg(**aggregation).reset_index()


def make_curated_template(
    sample_ids: list[str],
    dataset_id: str,
    cancer_type_default: str,
) -> dict[str, dict[str, Any]]:
    curated = {}

    for sample_id in sample_ids:
        curated[str(sample_id)] = {
            "patient_id": "REVIEW_REQUIRED",
            "dataset_id": dataset_id,
            "cancer_type": cancer_type_default or "REVIEW_REQUIRED",
            "tumor_site": "REVIEW_REQUIRED",
            "metastasis_site": "Unknown",
            "tumor_treatment": "REVIEW_REQUIRED",
            "cancer_site_origin": "REVIEW_REQUIRED",
            "tumour_grade": "Unknown",
            "tumour_stage": "Unknown",
            "histological_subtype": "REVIEW_REQUIRED",
            "patient_ethnicity": "Unknown",
        }

    return curated


def build_recipe(
    dataset_id: str,
    dataset_name: str,
    inspection_dir: Path,
    matrix_path: Path,
    metadata_path: Path,
    output_h5ad: Path,
    cancer_type_default: str,
) -> dict[str, Any]:
    classification = read_json(inspection_dir / "file_classification.json")
    data_size = read_json(inspection_dir / "data_size_policy.json")

    format_label = classification.get("format_label") or classification.get("detected_format")
    recommended_module = classification.get("recommended_module")

    metadata_inspection = inspect_metadata_table(metadata_path)

    if not metadata_inspection.candidate_cell_id_columns:
        raise ValueError("No candidate cell ID column found in metadata.")

    if not metadata_inspection.candidate_sample_columns:
        raise ValueError("No candidate sample column found in metadata.")

    cell_id_column = metadata_inspection.candidate_cell_id_columns[0]
    sample_id_column = metadata_inspection.candidate_sample_columns[0]

    metadata = load_metadata_preview(metadata_path)
    sample_summary = summarize_metadata_samples(metadata, sample_id_column)
    sample_ids = sample_summary[sample_id_column].astype(str).tolist()

    recipe = {
        "dataset_id": dataset_id,
        "dataset_name": dataset_name,
        "input_format": format_label,
        "reader": "sparse_global_csv" if format_label == "single_global_compressed_csv_count_matrix" else "REVIEW_REQUIRED",
        "matrix_path": str(matrix_path),
        "metadata_path": str(metadata_path),
        "output_h5ad": str(output_h5ad),
        "orientation": "gene_by_cell",
        "cell_id_column": cell_id_column,
        "sample_id_column": sample_id_column,
        "keep_metadata_cells_only": True,
        "standardize_obs_schema": True,
        "obs_mappings": {
            "sample_id": sample_id_column,
        },
        "expected": {
            "n_cells": "REVIEW_REQUIRED",
            "n_genes": "REVIEW_REQUIRED",
            "obs_columns": STANDARD_OBS_COLUMNS,
        },
        "curated_sample_metadata": make_curated_template(
            sample_ids=sample_ids,
            dataset_id=dataset_id,
            cancer_type_default=cancer_type_default,
        ),
        "review_required": {
            "fields": sorted(REVIEW_REQUIRED_FIELDS),
            "reason": "Biological metadata interpretation must be manually reviewed before final h5ad generation.",
        },
        "evidence": {
            "inspection_dir": str(inspection_dir),
            "detected_format": format_label,
            "recommended_module": recommended_module,
            "data_size_policy": data_size.get("policy_label") or data_size.get("decision_label"),
            "metadata_cell_id_column_candidates": metadata_inspection.candidate_cell_id_columns,
            "metadata_sample_column_candidates": metadata_inspection.candidate_sample_columns,
            "n_detected_samples": len(sample_ids),
        },
        "metadata_sample_summary": sample_summary.to_dict(orient="records"),
        "notes": [
            "This recipe is semi-automatic and must be reviewed before final use.",
            "Technical fields were inferred from inspection outputs.",
            "Curated biological metadata fields marked REVIEW_REQUIRED must be completed by a human reviewer.",
        ],
    }

    return recipe


def main() -> None:
    parser = argparse.ArgumentParser(description="Generate a semi-automatic GEO-to-h5ad recipe template.")

    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--dataset-name", default="")
    parser.add_argument("--inspection-dir", required=True)
    parser.add_argument("--matrix-path", required=True)
    parser.add_argument("--metadata-path", required=True)
    parser.add_argument("--output-h5ad", required=True)
    parser.add_argument("--output-recipe", required=True)
    parser.add_argument("--cancer-type-default", default="REVIEW_REQUIRED")

    args = parser.parse_args()

    recipe = build_recipe(
        dataset_id=args.dataset_id,
        dataset_name=args.dataset_name,
        inspection_dir=Path(args.inspection_dir),
        matrix_path=Path(args.matrix_path),
        metadata_path=Path(args.metadata_path),
        output_h5ad=Path(args.output_h5ad),
        cancer_type_default=args.cancer_type_default,
    )

    output_recipe = Path(args.output_recipe)
    output_recipe.parent.mkdir(parents=True, exist_ok=True)

    with output_recipe.open("w", encoding="utf-8") as handle:
        yaml.safe_dump(recipe, handle, sort_keys=False, allow_unicode=True)

    print("Wrote recipe:", output_recipe)
    print("Detected format:", recipe["input_format"])
    print("Cell ID column:", recipe["cell_id_column"])
    print("Sample ID column:", recipe["sample_id_column"])
    print("Detected samples:", recipe["evidence"]["n_detected_samples"])
    print("Review fields:", ", ".join(recipe["review_required"]["fields"]))


if __name__ == "__main__":
    main()
