from pathlib import Path
import shutil

import h5py
import numpy as np

from src.tenx_h5_standardizer import (
    detect_10x_h5_layout,
    inspect_many_10x_h5,
    summarize_inspections,
)


TMP = Path("data/test_tenx_h5_standardizer")


def reset_tmp():
    if TMP.exists():
        shutil.rmtree(TMP)
    TMP.mkdir(parents=True)


def write_modern_h5(path: Path):
    with h5py.File(path, "w") as handle:
        matrix = handle.create_group("matrix")
        matrix.create_dataset("barcodes", data=np.array([b"cell1-1", b"cell2-1"]))
        matrix.create_dataset("data", data=np.array([1, 2, 3]))
        matrix.create_dataset("indices", data=np.array([0, 1, 0]))
        matrix.create_dataset("indptr", data=np.array([0, 1, 3]))
        matrix.create_dataset("shape", data=np.array([3, 2]))

        features = matrix.create_group("features")
        features.create_dataset("id", data=np.array([b"ENSG1", b"ENSG2", b"ENSG3"]))
        features.create_dataset("name", data=np.array([b"GeneA", b"GeneB", b"GeneA"]))
        features.create_dataset(
            "feature_type",
            data=np.array([b"Gene Expression", b"Gene Expression", b"Gene Expression"]),
        )


def write_old_h5(path: Path):
    with h5py.File(path, "w") as handle:
        matrix = handle.create_group("hg19")
        matrix.create_dataset("barcodes", data=np.array([b"cell1-1", b"cell2-1", b"cell3-1"]))
        matrix.create_dataset("data", data=np.array([1, 2, 3, 4]))
        matrix.create_dataset("indices", data=np.array([0, 1, 0, 2]))
        matrix.create_dataset("indptr", data=np.array([0, 1, 3, 4]))
        matrix.create_dataset("shape", data=np.array([4, 3]))
        matrix.create_dataset("genes", data=np.array([b"ENSG1", b"ENSG2", b"ENSG3", b"ENSG4"]))
        matrix.create_dataset("gene_names", data=np.array([b"GeneA", b"GeneB", b"GeneB", b"GeneD"]))


def test_modern_layout():
    reset_tmp()
    path = TMP / "GSM1_filtered_feature_bc_matrix.h5"
    write_modern_h5(path)

    result = detect_10x_h5_layout(path)

    print("\n=== Modern 10x h5 layout ===")
    print(result)

    assert result.layout_label == "modern_10x_h5_v3_matrix_group"
    assert result.n_cells == 2
    assert result.n_genes == 3
    assert result.duplicated_gene_names == 1
    assert result.gene_name_source == "matrix/features/name"


def test_old_layout():
    reset_tmp()
    path = TMP / "GSM2_old_filtered_gene_bc_matrices_h5.h5"
    write_old_h5(path)

    result = detect_10x_h5_layout(path)

    print("\n=== Old 10x h5 layout ===")
    print(result)

    assert result.layout_label == "old_10x_h5_v2_genome_group"
    assert result.genome_group == "hg19"
    assert result.n_cells == 3
    assert result.n_genes == 4
    assert result.duplicated_gene_names == 1
    assert result.gene_name_source == "hg19/gene_names"


def test_many_summary():
    reset_tmp()
    modern = TMP / "GSM1_filtered_feature_bc_matrix.h5"
    old = TMP / "GSM2_old_filtered_gene_bc_matrices_h5.h5"

    write_modern_h5(modern)
    write_old_h5(old)

    inspections = inspect_many_10x_h5([modern, old])
    df = summarize_inspections(inspections)

    print("\n=== 10x h5 summary ===")
    print(df)

    assert len(df) == 2
    assert set(df["layout_label"]) == {
        "modern_10x_h5_v3_matrix_group",
        "old_10x_h5_v2_genome_group",
    }


def main():
    test_modern_layout()
    test_old_layout()
    test_many_summary()
    print("\nAll 10x h5 standardizer tests passed.")


if __name__ == "__main__":
    main()
