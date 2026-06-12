from pathlib import Path
import gzip
import shutil

from src.global_matrix_metadata_standardizer import (
    inspect_global_matrix,
    inspect_metadata_table,
    read_global_matrix_with_metadata,
)


TMP = Path("data/test_global_matrix_metadata_standardizer")


def reset_tmp():
    if TMP.exists():
        shutil.rmtree(TMP)
    TMP.mkdir(parents=True)


def write_gzip_text(path: Path, text: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text)


def write_gene_by_cell_global_csv(path: Path):
    text = "\n".join(
        [
            "Unnamed: 0,AAACCCAAGACGAGCT.1_1,AAACGAAAGAGAGAAC.1_2,AAACGAAAGCTTCTAG.1_10",
            "GeneA,1,0,2",
            "GeneB,0,3,1",
            "GeneC,2,1,0",
        ]
    ) + "\n"
    write_gzip_text(path, text)


def write_cell_metadata_tsv(path: Path):
    text = "\n".join(
        [
            "cell_id\tsample\tcell_type",
            "AAACCCAAGACGAGCT.1_1\tPt1\tTumor",
            "AAACGAAAGAGAGAAC.1_2\tPt2\tImmune",
        ]
    ) + "\n"
    write_gzip_text(path, text)


def write_cell_metadata_with_sample_prefixed_barcodes_csv(path: Path):
    text = "\n".join(
        [
            "Unnamed: 0,samples,Condition,Location",
            "B_cac10_AAACCTGAGTCAATAG,B_cac10,Normal,Left",
            "B_cac10_AAACCTGCACAGCCCA,B_cac10,Normal,Left",
            "B_cac11_AAACCTGCACTTCGAA,B_cac11,Tumor,Right",
        ]
    ) + "\n"
    write_gzip_text(path, text)


def write_cell_metadata_with_unnamed_index_csv(path: Path):
    text = "\n".join(
        [
            "Unnamed: 0,samples,Condition,Location",
            "AAACCCAAGACGAGCT.1_1,Pt1,Tumor,Primary",
            "AAACGAAAGAGAGAAC.1_2,Pt2,Normal,Adjacent",
        ]
    ) + "\n"
    write_gzip_text(path, text)


def test_inspect_global_matrix_suffixes():
    reset_tmp()
    matrix_path = TMP / "GSE281120_counts.csv.gz"
    write_gene_by_cell_global_csv(matrix_path)

    result = inspect_global_matrix(matrix_path)

    print("\n=== Global matrix suffix inspection ===")
    print(result)

    assert result.first_column_role == "feature_id"
    assert result.orientation_label == "gene_by_cell_global_matrix"
    assert result.transpose_required_for_anndata is True
    assert result.barcode_suffix_summary == {"1": 1, "10": 1, "2": 1}


def test_inspect_metadata_table():
    reset_tmp()
    metadata_path = TMP / "cell_metadata.tsv.gz"
    write_cell_metadata_tsv(metadata_path)

    result = inspect_metadata_table(metadata_path)

    print("\n=== Metadata table inspection ===")
    print(result)

    assert "cell_id" in result.candidate_cell_id_columns
    assert "sample" in result.candidate_sample_columns




def test_metadata_unnamed_index_detected_as_cell_id():
    reset_tmp()
    metadata_path = TMP / "cell_metadata_index.csv.gz"
    write_cell_metadata_with_unnamed_index_csv(metadata_path)

    result = inspect_metadata_table(metadata_path)

    print("\n=== Metadata table unnamed index inspection ===")
    print(result)

    assert "Unnamed: 0" in result.candidate_cell_id_columns
    assert "samples" in result.candidate_sample_columns
    assert not result.warnings




def test_metadata_sample_prefixed_barcodes_detected_as_cell_id():
    reset_tmp()
    metadata_path = TMP / "cell_metadata_sample_prefixed.csv.gz"
    write_cell_metadata_with_sample_prefixed_barcodes_csv(metadata_path)

    result = inspect_metadata_table(metadata_path)

    print("\n=== Metadata sample-prefixed barcode inspection ===")
    print(result)

    assert "Unnamed: 0" in result.candidate_cell_id_columns
    assert "samples" in result.candidate_sample_columns
    assert not result.warnings


def test_read_global_matrix_with_metadata():
    reset_tmp()
    matrix_path = TMP / "GSE281120_counts.csv.gz"
    metadata_path = TMP / "cell_metadata.tsv.gz"

    write_gene_by_cell_global_csv(matrix_path)
    write_cell_metadata_tsv(metadata_path)

    adata, inspection, coverage, validation = read_global_matrix_with_metadata(
        matrix_path,
        metadata_path=metadata_path,
        cell_id_column="cell_id",
        sample_id_column="sample",
        keep_metadata_cells_only=True,
    )

    print("\n=== Read global matrix with metadata ===")
    print("adata shape:", adata.shape)
    print("obs columns:", list(adata.obs.columns))
    print("coverage:", coverage)
    print("validation:", validation)

    assert adata.shape == (2, 3)
    assert "sample_id" in adata.obs.columns
    assert validation.passed is True
    assert coverage.evidence["n_common_cells"] == 2


def main():
    test_inspect_global_matrix_suffixes()
    test_inspect_metadata_table()
    test_metadata_unnamed_index_detected_as_cell_id()
    test_metadata_sample_prefixed_barcodes_detected_as_cell_id()
    test_read_global_matrix_with_metadata()
    print("\nAll global matrix metadata standardizer tests passed.")


if __name__ == "__main__":
    main()
