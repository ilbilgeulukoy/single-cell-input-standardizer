from pathlib import Path
import shutil
import subprocess
import sys

import numpy as np
from scipy import sparse
from scipy.io import mmwrite
import gzip


TMP = Path("data/test_inspect_input_cli")
OUT = Path("reports/test_inspect_input_cli")


def reset_tmp():
    if TMP.exists():
        shutil.rmtree(TMP)
    if OUT.exists():
        shutil.rmtree(OUT)

    TMP.mkdir(parents=True)
    OUT.mkdir(parents=True)


def write_gzip_text(path: Path, text: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text)


def write_gzip_mtx(path: Path, matrix):
    plain_path = path.with_suffix("")
    mmwrite(plain_path, matrix)

    with open(plain_path, "rb") as source, gzip.open(path, "wb") as target:
        target.write(source.read())

    plain_path.unlink()


def write_10x_mtx_example():
    matrix = sparse.coo_matrix(
        np.array(
            [
                [1, 0, 2],
                [0, 3, 0],
                [4, 0, 5],
            ]
        )
    )

    write_gzip_mtx(TMP / "GSM1_matrix.mtx.gz", matrix)

    features_text = "\n".join(
        [
            "ENSG1\tGeneA\tGene Expression",
            "ENSG2\tGeneB\tGene Expression",
            "ENSG3\tGeneA\tGene Expression",
        ]
    ) + "\n"
    write_gzip_text(TMP / "GSM1_features.tsv.gz", features_text)

    barcodes_text = "\n".join(["cell1-1", "cell2-1", "cell3-1"]) + "\n"
    write_gzip_text(TMP / "GSM1_barcodes.tsv.gz", barcodes_text)


def test_inspect_input_cli_on_10x_mtx():
    reset_tmp()
    write_10x_mtx_example()

    result = subprocess.run(
        [
            sys.executable,
            "scripts/inspect_input.py",
            str(TMP),
            "--output",
            str(OUT),
        ],
        capture_output=True,
        text=True,
    )

    print("\n=== CLI stdout ===")
    print(result.stdout)

    print("\n=== CLI stderr ===")
    print(result.stderr)

    assert result.returncode == 0
    assert (OUT / "file_classification.json").exists()
    assert (OUT / "data_size_policy.json").exists()
    assert (OUT / "modality_filter.json").exists()
    assert (OUT / "format_inspection_summary.csv").exists()
    assert (OUT / "format_inspection_details.json").exists()
    assert (OUT / "manual_review_report.md").exists()
    assert (OUT / "inspection_manifest.json").exists()


def main():
    test_inspect_input_cli_on_10x_mtx()
    print("\nAll inspect_input CLI tests passed.")


if __name__ == "__main__":
    main()
