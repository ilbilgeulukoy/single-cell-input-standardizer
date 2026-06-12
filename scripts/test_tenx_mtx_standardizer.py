from pathlib import Path
import gzip
import shutil

import numpy as np
from scipy import sparse
from scipy.io import mmwrite

from src.tenx_mtx_standardizer import (
    discover_10x_mtx_triplets,
    inspect_10x_mtx_triplet,
    inspect_many_10x_mtx_triplets,
    summarize_mtx_inspections,
)


TMP = Path("data/test_tenx_mtx_standardizer")


def reset_tmp():
    if TMP.exists():
        shutil.rmtree(TMP)
    TMP.mkdir(parents=True)


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


def write_triplet(sample: str, use_features: bool = True):
    matrix = sparse.coo_matrix(
        np.array(
            [
                [1, 0, 2],
                [0, 3, 0],
                [4, 0, 5],
                [0, 6, 0],
            ]
        )
    )

    write_gzip_mtx(TMP / f"{sample}_matrix.mtx.gz", matrix)

    if use_features:
        features_text = "\n".join(
            [
                "ENSG1\tGeneA\tGene Expression",
                "ENSG2\tGeneB\tGene Expression",
                "ENSG3\tGeneA\tGene Expression",
                "ENSG4\tGeneD\tGene Expression",
            ]
        ) + "\n"
        write_gzip_text(TMP / f"{sample}_features.tsv.gz", features_text)
    else:
        genes_text = "\n".join(
            [
                "ENSG1\tGeneA",
                "ENSG2\tGeneB",
                "ENSG3\tGeneA",
                "ENSG4\tGeneD",
            ]
        ) + "\n"
        write_gzip_text(TMP / f"{sample}_genes.tsv.gz", genes_text)

    barcodes_text = "\n".join(["cell1-1", "cell2-1", "cell3-1"]) + "\n"
    write_gzip_text(TMP / f"{sample}_barcodes.tsv.gz", barcodes_text)


def test_discover_triplets_features():
    reset_tmp()
    write_triplet("GSM1_sample", use_features=True)

    triplets = discover_10x_mtx_triplets(TMP)

    print("\n=== Discover triplets with features.tsv ===")
    print(triplets)

    assert len(triplets) == 1
    assert triplets[0].sample_id == "GSM1"
    assert triplets[0].feature_table_type == "features.tsv"


def test_inspect_triplet_features():
    reset_tmp()
    write_triplet("GSM1_sample", use_features=True)

    triplet = discover_10x_mtx_triplets(TMP)[0]
    result = inspect_10x_mtx_triplet(triplet)

    print("\n=== Inspect triplet with features.tsv ===")
    print(result)

    assert result.raw_matrix_shape == (4, 3)
    assert result.expected_anndata_shape == (3, 4)
    assert result.n_features == 4
    assert result.n_barcodes == 3
    assert result.gene_symbol_column == 1
    assert result.feature_type_column == 2
    assert result.duplicated_gene_symbols == 1
    assert result.transpose_required_for_anndata is True


def test_discover_triplets_genes():
    reset_tmp()
    write_triplet("GSM2_sample", use_features=False)

    triplets = discover_10x_mtx_triplets(TMP)

    print("\n=== Discover triplets with genes.tsv ===")
    print(triplets)

    assert len(triplets) == 1
    assert triplets[0].sample_id == "GSM2"
    assert triplets[0].feature_table_type == "genes.tsv"


def test_many_summary():
    reset_tmp()
    write_triplet("GSM1_sample", use_features=True)
    write_triplet("GSM2_sample", use_features=False)

    triplets = discover_10x_mtx_triplets(TMP)
    inspections = inspect_many_10x_mtx_triplets(triplets)
    df = summarize_mtx_inspections(inspections)

    print("\n=== MTX summary ===")
    print(df)

    assert len(df) == 2
    assert set(df["sample_id"]) == {"GSM1", "GSM2"}


def main():
    test_discover_triplets_features()
    test_inspect_triplet_features()
    test_discover_triplets_genes()
    test_many_summary()
    print("\nAll 10x MTX standardizer tests passed.")


if __name__ == "__main__":
    main()
