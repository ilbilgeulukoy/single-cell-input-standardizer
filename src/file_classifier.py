from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


@dataclass
class FileClassification:
    input_path: str
    total_files: int
    total_directories: int
    format_label: str
    recommended_module: str
    recommended_mode: str
    selected_files: list[str] = field(default_factory=list)
    excluded_files: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    evidence: dict[str, object] = field(default_factory=dict)


def _list_files(input_path: Path) -> list[Path]:
    if input_path.is_file():
        return [input_path]

    return sorted([p for p in input_path.rglob("*") if p.is_file()])


def _list_directories(input_path: Path) -> list[Path]:
    if input_path.is_file():
        return []

    return sorted([p for p in input_path.rglob("*") if p.is_dir()])


def _rel(path: Path, root: Path) -> str:
    try:
        return str(path.relative_to(root))
    except ValueError:
        return str(path)


def _count_suffix(files: Iterable[Path], suffix: str) -> int:
    return sum(1 for p in files if p.name.endswith(suffix))


def _has_any(files: Iterable[Path], suffixes: tuple[str, ...]) -> bool:
    return any(p.name.endswith(suffixes) for p in files)


def _gsm_prefix(filename: str) -> str | None:
    if not filename.startswith("GSM"):
        return None
    return filename.split("_")[0]


def _classify_10x_mtx_triplets(files: list[Path]) -> dict[str, dict[str, list[Path]]]:
    triplets: dict[str, dict[str, list[Path]]] = {}

    for path in files:
        name = path.name
        gsm = _gsm_prefix(name)
        if gsm is None:
            continue

        triplets.setdefault(gsm, {"matrix": [], "features": [], "barcodes": []})

        if name.endswith(".mtx.gz") and ("matrix" in name):
            triplets[gsm]["matrix"].append(path)
        elif name.endswith(".tsv.gz") and ("features" in name or "genes" in name):
            triplets[gsm]["features"].append(path)
        elif name.endswith(".tsv.gz") and "barcodes" in name:
            triplets[gsm]["barcodes"].append(path)

    return triplets


def _complete_triplets(triplets: dict[str, dict[str, list[Path]]]) -> dict[str, dict[str, list[Path]]]:
    complete = {}

    for gsm, roles in triplets.items():
        if len(roles["matrix"]) == 1 and len(roles["features"]) == 1 and len(roles["barcodes"]) == 1:
            complete[gsm] = roles

    return complete


def classify_input(input_path: str | Path) -> FileClassification:
    root = Path(input_path)
    files = _list_files(root)
    directories = _list_directories(root)

    names = [p.name for p in files]

    h5_files = [p for p in files if p.name.endswith(".h5")]
    h5ad_files = [p for p in files if p.name.endswith((".h5ad", ".h5ad.gz"))]
    csv_gz_files = [p for p in files if p.name.endswith(".csv.gz")]
    csv_files = [p for p in files if p.name.endswith(".csv")]
    txt_gz_files = [p for p in files if p.name.endswith(".txt.gz")]
    txt_files = [p for p in files if p.name.endswith(".txt")]
    tar_files = [p for p in files if p.name.endswith((".tar", ".tar.gz", ".tgz"))]
    spatial_zip_files = [p for p in files if p.name.endswith("spatial.zip")]
    bulk_like_files = [p for p in files if "bulk" in p.name.lower()]

    triplets = _classify_10x_mtx_triplets(files)
    complete_triplets = _complete_triplets(triplets)

    warnings: list[str] = []
    excluded_files: list[str] = []

    evidence = {
        "h5_files": len(h5_files),
        "h5ad_files": len(h5ad_files),
        "csv_gz_files": len(csv_gz_files),
        "csv_files": len(csv_files),
        "txt_gz_files": len(txt_gz_files),
        "txt_files": len(txt_files),
        "tar_files": len(tar_files),
        "spatial_zip_files": len(spatial_zip_files),
        "bulk_like_files": len(bulk_like_files),
        "complete_mtx_triplets": len(complete_triplets),
        "detected_gsms_with_triplet_roles": sorted(triplets),
    }

    if h5ad_files:
        return FileClassification(
            input_path=str(root),
            total_files=len(files),
            total_directories=len(directories),
            format_label="h5ad_or_h5ad_gz_matrix",
            recommended_module="src/h5ad_inspector.py",
            recommended_mode="server_or_local_depending_on_size",
            selected_files=[_rel(p, root) for p in h5ad_files],
            warnings=["Existing h5ad-like object detected. Inspect obs/var before standardization."],
            evidence=evidence,
        )

    if h5_files:
        return FileClassification(
            input_path=str(root),
            total_files=len(files),
            total_directories=len(directories),
            format_label="per_sample_10x_h5_files",
            recommended_module="src/tenx_h5_standardizer.py",
            recommended_mode="local_temporary_download",
            selected_files=[_rel(p, root) for p in h5_files],
            warnings=[
                "10x h5 files can use either modern /matrix layout or old genome-group layout such as /hg19.",
                "Prefix barcodes with GSM/sample ID before concatenation.",
            ],
            evidence=evidence,
        )

    if complete_triplets:
        selected = []
        for roles in complete_triplets.values():
            selected.extend(roles["matrix"])
            selected.extend(roles["features"])
            selected.extend(roles["barcodes"])

        if spatial_zip_files:
            warnings.append("Spatial Visium files detected. Filter modality before scRNA-seq standardization.")
            excluded_files.extend([_rel(p, root) for p in spatial_zip_files])

        if bulk_like_files:
            warnings.append("Bulk-like files detected. Exclude bulk files from single-cell triplet processing.")
            excluded_files.extend([_rel(p, root) for p in bulk_like_files])

        return FileClassification(
            input_path=str(root),
            total_files=len(files),
            total_directories=len(directories),
            format_label="per_sample_10x_matrix_market_triplets",
            recommended_module="src/tenx_mtx_standardizer.py",
            recommended_mode="local_temporary_download",
            selected_files=[_rel(p, root) for p in selected],
            excluded_files=excluded_files,
            warnings=warnings,
            evidence=evidence,
        )

    if txt_gz_files:
        return FileClassification(
            input_path=str(root),
            total_files=len(files),
            total_directories=len(directories),
            format_label="compressed_txt_count_tables",
            recommended_module="src/count_table_standardizer.py",
            recommended_mode="local_temporary_download",
            selected_files=[_rel(p, root) for p in txt_gz_files],
            excluded_files=[_rel(p, root) for p in spatial_zip_files],
            warnings=[
                "Read .txt.gz directly when possible. Avoid unnecessary decompression.",
                "Orientation must be inferred from row index, first column and cell barcode-like columns.",
            ],
            evidence=evidence,
        )

    if csv_gz_files:
        if len(csv_gz_files) == 1:
            module = "src/global_matrix_metadata_standardizer.py"
            format_label = "single_global_compressed_csv_count_matrix"
            warnings.append("Single global CSV detected. Sample IDs may be encoded in cell barcodes or external metadata.")
        else:
            module = "src/count_table_standardizer.py"
            format_label = "compressed_csv_count_tables"
            warnings.append("Multiple CSV count tables detected. Validate orientation per file.")

        return FileClassification(
            input_path=str(root),
            total_files=len(files),
            total_directories=len(directories),
            format_label=format_label,
            recommended_module=module,
            recommended_mode="local_temporary_download",
            selected_files=[_rel(p, root) for p in csv_gz_files],
            warnings=warnings,
            evidence=evidence,
        )

    if csv_files:
        return FileClassification(
            input_path=str(root),
            total_files=len(files),
            total_directories=len(directories),
            format_label="csv_count_tables",
            recommended_module="src/count_table_standardizer.py",
            recommended_mode="local_temporary_download",
            selected_files=[_rel(p, root) for p in csv_files],
            warnings=["Validate whether rows are cells or genes before AnnData construction."],
            evidence=evidence,
        )

    if tar_files:
        return FileClassification(
            input_path=str(root),
            total_files=len(files),
            total_directories=len(directories),
            format_label="archive_or_nested_archive",
            recommended_module="src/nested_archive_extractor.py",
            recommended_mode="local_temporary_download",
            selected_files=[_rel(p, root) for p in tar_files],
            warnings=[
                "Archive detected. Extract to temporary directory and classify extracted files.",
                "Nested GSM-level tar archives may hide 10x triplets.",
            ],
            evidence=evidence,
        )

    return FileClassification(
        input_path=str(root),
        total_files=len(files),
        total_directories=len(directories),
        format_label="unknown",
        recommended_module="manual_review_required",
        recommended_mode="manual_review",
        selected_files=[_rel(p, root) for p in files[:20]],
        warnings=["No known single-cell input pattern detected."],
        evidence=evidence,
    )


def print_classification(result: FileClassification) -> None:
    print("=== File classification ===")
    print("Input path:", result.input_path)
    print("Total files:", result.total_files)
    print("Total directories:", result.total_directories)
    print("Format label:", result.format_label)
    print("Recommended module:", result.recommended_module)
    print("Recommended mode:", result.recommended_mode)

    print("\nSelected files:")
    for item in result.selected_files[:30]:
        print("-", item)

    if len(result.selected_files) > 30:
        print(f"... {len(result.selected_files) - 30} more selected files")

    print("\nExcluded files:")
    for item in result.excluded_files[:30]:
        print("-", item)

    if result.warnings:
        print("\nWarnings:")
        for warning in result.warnings:
            print("-", warning)

    print("\nEvidence:")
    for key, value in result.evidence.items():
        print(f"- {key}: {value}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Classify single-cell input files.")
    parser.add_argument("input_path", help="Input file or directory to classify.")

    args = parser.parse_args()
    classification = classify_input(args.input_path)
    print_classification(classification)
