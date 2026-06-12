from pathlib import Path
import shutil

from src.file_classifier import classify_input


TMP = Path("data/test_file_classifier")


def reset_tmp():
    if TMP.exists():
        shutil.rmtree(TMP)
    TMP.mkdir(parents=True)


def touch(relative_path: str):
    path = TMP / relative_path
    path.parent.mkdir(parents=True, exist_ok=True)
    path.touch()


def run_case(name: str, files: list[str], expected_format: str):
    reset_tmp()

    for file in files:
        touch(file)

    result = classify_input(TMP)

    print("\n===", name, "===")
    print("Detected:", result.format_label)
    print("Expected:", expected_format)
    print("Module:", result.recommended_module)

    assert result.format_label == expected_format, (
        f"{name}: expected {expected_format}, got {result.format_label}"
    )


def main():
    run_case(
        "per-sample 10x MTX triplet",
        [
            "GSM1_sample_matrix.mtx.gz",
            "GSM1_sample_features.tsv.gz",
            "GSM1_sample_barcodes.tsv.gz",
        ],
        "per_sample_10x_matrix_market_triplets",
    )

    run_case(
        "10x h5",
        [
            "GSM1_filtered_feature_bc_matrix.h5",
            "GSM2_filtered_feature_bc_matrix.h5",
        ],
        "per_sample_10x_h5_files",
    )

    run_case(
        "compressed txt",
        [
            "GSM1_counts_A.txt.gz",
            "GSM2_counts_B.txt.gz",
        ],
        "compressed_txt_count_tables",
    )

    run_case(
        "single global csv gz",
        [
            "GSE281120_counts.csv.gz",
        ],
        "single_global_compressed_csv_count_matrix",
    )

    run_case(
        "archive",
        [
            "GSE168652_RAW.tar",
        ],
        "archive_or_nested_archive",
    )

    print("\nAll file classifier tests passed.")


if __name__ == "__main__":
    main()
