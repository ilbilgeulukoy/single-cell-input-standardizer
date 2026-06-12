from pathlib import Path
import yaml
import pandas as pd


CASES_DIR = Path("cases")
OUTPUT_DIR = Path("reports")
OUTPUT_DIR.mkdir(exist_ok=True)


def load_yaml(path: Path) -> dict:
    with path.open("r", encoding="utf-8") as handle:
        return yaml.safe_load(handle) or {}


def get_nested(data: dict, keys: list[str], default=None):
    current = data
    for key in keys:
        if not isinstance(current, dict):
            return default
        current = current.get(key)
        if current is None:
            return default
    return current


def main():
    rows = []

    for path in sorted(CASES_DIR.glob("case_*.yaml")):
        data = load_yaml(path)

        row = {
            "case_id": data.get("case_id"),
            "paper_id": data.get("paper_id"),
            "geo_accession": data.get("geo_accession"),
            "format_label": get_nested(data, ["detected_input_format", "format_label"]),
            "selected_format": get_nested(data, ["detected_input_format", "selected_format_for_standardization"]),
            "recommended_mode": get_nested(data, ["input_source", "recommended_processing_mode"]),
            "total_cells": get_nested(data, ["sample_summary", "total_cells"]),
            "selected_samples": get_nested(data, ["sample_summary", "selected_samples"]),
            "number_of_samples": get_nested(data, ["sample_summary", "number_of_samples"]),
            "unique_gene_counts": get_nested(data, ["sample_summary", "unique_gene_counts"]),
            "transpose_required": get_nested(data, ["detected_input_format", "transpose_required_for_anndata"]),
            "metadata_required": get_nested(data, ["metadata_strategy", "curated_sample_metadata_required"]),
            "gene_mapping_key": get_nested(data, ["gene_mapping_strategy", "mapping_key"]),
            "main_module": get_nested(data, ["tool_design_implication", "reusable_module_needed"]),
            "source_file": str(path),
        }

        rows.append(row)

    df = pd.DataFrame(rows)

    csv_path = OUTPUT_DIR / "case_atlas_summary.csv"
    md_path = OUTPUT_DIR / "case_atlas_summary.md"

    df.to_csv(csv_path, index=False)
    df.to_markdown(md_path, index=False)

    print("Cases summarized:", len(df))
    print("Wrote:", csv_path)
    print("Wrote:", md_path)

    print("\n=== Format counts ===")
    print(df["format_label"].value_counts(dropna=False))

    print("\n=== Main modules ===")
    print(df["main_module"].value_counts(dropna=False))


if __name__ == "__main__":
    main()
