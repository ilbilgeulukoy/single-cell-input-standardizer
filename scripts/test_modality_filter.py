from pathlib import Path
import shutil

from src.modality_filter import classify_file_modality, filter_by_modality


TMP = Path("data/test_modality_filter")


def reset_tmp():
    if TMP.exists():
        shutil.rmtree(TMP)
    TMP.mkdir(parents=True)


def touch(name: str):
    path = TMP / name
    path.parent.mkdir(parents=True, exist_ok=True)
    path.touch()


def test_classify_file_modality():
    assert classify_file_modality("GSM1_matrix.mtx.gz") == "single_cell_rna_seq"
    assert classify_file_modality("GSM1_spatial.zip") == "spatial_visium"
    assert classify_file_modality("GSM1_fragments.tsv.gz") == "single_cell_atac"
    assert classify_file_modality("bulk_ReadsPerGene.out.tab") == "bulk_rna_seq"


def test_filter_by_modality():
    reset_tmp()
    touch("GSM1_matrix.mtx.gz")
    touch("GSM1_features.tsv.gz")
    touch("GSM1_barcodes.tsv.gz")
    touch("GSM1_spatial.zip")
    touch("bulk_ReadsPerGene.out.tab")

    result = filter_by_modality(TMP)

    print("\n=== Modality filter ===")
    print(result)

    assert result.modality_counts["single_cell_rna_seq"] == 3
    assert result.modality_counts["spatial_visium"] == 1
    assert result.modality_counts["bulk_rna_seq"] == 1
    assert len(result.selected_files) == 3
    assert result.warnings


def main():
    test_classify_file_modality()
    test_filter_by_modality()
    print("\nAll modality filter tests passed.")


if __name__ == "__main__":
    main()
