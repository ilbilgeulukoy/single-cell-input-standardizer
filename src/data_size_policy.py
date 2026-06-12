from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path


MB = 1024 * 1024
GB = 1024 * MB


@dataclass
class DataSizeDecision:
    input_path: str
    total_size_bytes: int
    total_size_mb: float
    total_size_gb: float
    policy_label: str
    recommended_mode: str
    allow_local_processing: bool
    allow_decompression: bool
    warnings: list[str] = field(default_factory=list)
    evidence: dict[str, object] = field(default_factory=dict)


def _iter_files(path: Path) -> list[Path]:
    if path.is_file():
        return [path]

    return sorted([p for p in path.rglob("*") if p.is_file()])


def _file_size(path: Path) -> int:
    try:
        return path.stat().st_size
    except FileNotFoundError:
        return 0


def evaluate_data_size(
    input_path: str | Path,
    large_file_mb: float = 500,
    server_recommended_gb: float = 5,
    server_required_gb: float = 20,
) -> DataSizeDecision:
    path = Path(input_path)
    files = _iter_files(path)

    total_size_bytes = sum(_file_size(p) for p in files)
    total_size_mb = round(total_size_bytes / MB, 2)
    total_size_gb = round(total_size_bytes / GB, 2)

    largest_files = sorted(files, key=_file_size, reverse=True)[:10]

    compressed_files = [
        p for p in files
        if p.name.endswith((".gz", ".zip", ".tar", ".tar.gz", ".tgz"))
    ]

    h5ad_files = [
        p for p in files
        if p.name.endswith((".h5ad", ".h5ad.gz"))
    ]

    warnings: list[str] = []

    allow_local_processing = True
    allow_decompression = True
    recommended_mode = "local_temporary_processing"
    policy_label = "small_or_moderate_input"

    if total_size_gb >= server_required_gb:
        policy_label = "server_required"
        recommended_mode = "server_or_hpc_required"
        allow_local_processing = False
        allow_decompression = False
        warnings.append("Input is very large. Use server/HPC processing.")
        warnings.append("Avoid full local download, full extraction, and dense loading.")
    elif total_size_gb >= server_recommended_gb:
        policy_label = "server_recommended"
        recommended_mode = "server_or_hpc_recommended"
        allow_local_processing = False
        allow_decompression = False
        warnings.append("Input is large. Server/HPC processing is recommended.")
        warnings.append("Avoid unnecessary decompression and inspect files incrementally.")
    elif total_size_mb >= large_file_mb:
        policy_label = "large_local_with_caution"
        recommended_mode = "large_local_temporary_processing"
        allow_local_processing = True
        allow_decompression = False
        warnings.append("Input is large but may be processed locally with caution.")
        warnings.append("Avoid unnecessary decompression.")
        warnings.append("Prefer header inspection, streaming, chunking, or sparse-aware loading.")

    if compressed_files and not allow_decompression:
        warnings.append("Compressed files detected. Read compressed files directly when possible.")

    if h5ad_files and total_size_gb >= server_recommended_gb:
        warnings.append("Large h5ad-like file detected. Prefer backed/server inspection before loading into memory.")

    evidence = {
        "n_files": len(files),
        "n_compressed_files": len(compressed_files),
        "n_h5ad_files": len(h5ad_files),
        "largest_files": [
            {
                "path": str(p),
                "size_mb": round(_file_size(p) / MB, 2),
            }
            for p in largest_files
        ],
    }

    return DataSizeDecision(
        input_path=str(path),
        total_size_bytes=total_size_bytes,
        total_size_mb=total_size_mb,
        total_size_gb=total_size_gb,
        policy_label=policy_label,
        recommended_mode=recommended_mode,
        allow_local_processing=allow_local_processing,
        allow_decompression=allow_decompression,
        warnings=warnings,
        evidence=evidence,
    )


def print_data_size_decision(decision: DataSizeDecision) -> None:
    print("=== Data size policy ===")
    print("Input path:", decision.input_path)
    print("Total size MB:", decision.total_size_mb)
    print("Total size GB:", decision.total_size_gb)
    print("Policy label:", decision.policy_label)
    print("Recommended mode:", decision.recommended_mode)
    print("Allow local processing:", decision.allow_local_processing)
    print("Allow decompression:", decision.allow_decompression)

    if decision.warnings:
        print("\nWarnings:")
        for warning in decision.warnings:
            print("-", warning)

    print("\nEvidence:")
    for key, value in decision.evidence.items():
        print(f"- {key}: {value}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Evaluate data size policy for input files.")
    parser.add_argument("input_path", help="Input file or directory.")
    parser.add_argument("--large-file-mb", type=float, default=500)
    parser.add_argument("--server-recommended-gb", type=float, default=5)
    parser.add_argument("--server-required-gb", type=float, default=20)

    args = parser.parse_args()

    decision = evaluate_data_size(
        args.input_path,
        large_file_mb=args.large_file_mb,
        server_recommended_gb=args.server_recommended_gb,
        server_required_gb=args.server_required_gb,
    )
    print_data_size_decision(decision)
