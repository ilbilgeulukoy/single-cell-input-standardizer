from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any
import tarfile
import zipfile
import shutil


@dataclass
class ArchiveExtractionResult:
    archive_path: str
    output_dir: str
    extracted_files: list[str]
    nested_archives: list[str]
    warnings: list[str] = field(default_factory=list)
    evidence: dict[str, Any] = field(default_factory=dict)


def is_archive(path: str | Path) -> bool:
    name = Path(path).name.lower()
    return name.endswith((".tar", ".tar.gz", ".tgz", ".zip"))


def find_archives(input_path: str | Path) -> list[Path]:
    root = Path(input_path)
    if root.is_file():
        return [root] if is_archive(root) else []

    return sorted([p for p in root.rglob("*") if p.is_file() and is_archive(p)])


def safe_extract_archive(
    archive_path: str | Path,
    output_dir: str | Path,
    overwrite: bool = False,
) -> ArchiveExtractionResult:
    archive_path = Path(archive_path)
    output_dir = Path(output_dir)

    warnings: list[str] = []

    if output_dir.exists() and overwrite:
        shutil.rmtree(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)

    if not archive_path.exists():
        raise FileNotFoundError(f"Archive not found: {archive_path}")

    if archive_path.name.endswith((".tar", ".tar.gz", ".tgz")):
        with tarfile.open(archive_path, "r:*") as tar:
            members = tar.getmembers()

            for member in members:
                target = output_dir / member.name
                if not str(target.resolve()).startswith(str(output_dir.resolve())):
                    raise RuntimeError(f"Unsafe archive member path: {member.name}")

            tar.extractall(output_dir)

    elif archive_path.name.endswith(".zip"):
        with zipfile.ZipFile(archive_path, "r") as zf:
            for member in zf.namelist():
                target = output_dir / member
                if not str(target.resolve()).startswith(str(output_dir.resolve())):
                    raise RuntimeError(f"Unsafe archive member path: {member}")
            zf.extractall(output_dir)

    else:
        raise ValueError(f"Unsupported archive format: {archive_path}")

    extracted_files = sorted([str(p) for p in output_dir.rglob("*") if p.is_file()])
    nested_archives = sorted([str(p) for p in output_dir.rglob("*") if p.is_file() and is_archive(p)])

    if nested_archives:
        warnings.append("Nested archives detected. A second extraction/classification step is required.")

    return ArchiveExtractionResult(
        archive_path=str(archive_path),
        output_dir=str(output_dir),
        extracted_files=extracted_files,
        nested_archives=nested_archives,
        warnings=warnings,
        evidence={
            "n_extracted_files": len(extracted_files),
            "n_nested_archives": len(nested_archives),
        },
    )


def extract_nested_archives_once(
    input_dir: str | Path,
    output_base_dir: str | Path,
) -> list[ArchiveExtractionResult]:
    input_dir = Path(input_dir)
    output_base_dir = Path(output_base_dir)
    archives = find_archives(input_dir)

    results = []

    for archive in archives:
        clean_name = archive.name.replace(".tar.gz", "").replace(".tgz", "").replace(".tar", "").replace(".zip", "")
        outdir = output_base_dir / clean_name
        results.append(safe_extract_archive(archive, outdir, overwrite=True))

    return results


def print_archive_result(result: ArchiveExtractionResult) -> None:
    print("=== Archive extraction result ===")
    print("Archive:", result.archive_path)
    print("Output dir:", result.output_dir)
    print("Extracted files:", len(result.extracted_files))
    print("Nested archives:", len(result.nested_archives))

    if result.nested_archives:
        print("\nNested archives:")
        for item in result.nested_archives[:20]:
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

    parser = argparse.ArgumentParser(description="Safely extract archives and detect nested archives.")
    parser.add_argument("archive_path")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--overwrite", action="store_true")

    args = parser.parse_args()

    result = safe_extract_archive(args.archive_path, args.output_dir, overwrite=args.overwrite)
    print_archive_result(result)
