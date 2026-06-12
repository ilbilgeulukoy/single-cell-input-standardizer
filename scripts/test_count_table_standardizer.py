from pathlib import Path
import gzip
import shutil

import pandas as pd

from src.count_table_standardizer import (
    discover_count_tables,
    inspect_count_table,
    inspect_many_count_tables,
    summarize_count_table_inspections,
)


TMP = Path("data/test_count_table_standardizer")


def reset_tmp():
    if TMP.exists():
        shutil.rmtree(TMP)
    TMP.mkdir(parents=True)


def write_gzip_text(path: Path, text: str):
    path.parent.mkdir(parents=True, exist_ok=True)
    with gzip.open(path, "wt", encoding="utf-8") as handle:
        handle.write(text)


def write_gene_by_cell_csv(path: Path):
    text = "\n".join(
        [
            "Unnamed: 0,AAACCTGAGAGTAAGG-1,AAACCTGAGATCCTGT-1",
            "GeneA,1,0",
            "GeneB,0,3",
            "GeneC,2,1",
        ]
    ) + "\n"
    write_gzip_text(path, text)


def write_cell_by_gene_csv(path: Path):
    text = "\n".join(
        [
            "CellId,GeneA,GeneB,GeneC",
            "AAACCTGAGAGTAAGG-1,1,0,2",
            "AAACCTGAGATCCTGT-1,0,3,1",
        ]
    ) + "\n"
    write_gzip_text(path, text)


def write_ensembl_txt(path: Path):
    text = "\n".join(
        [
            "GENE\tAAACCTGAGAGTAAGG-1\tAAACCTGAGATCCTGT-1",
            "ENSG000001\t1\t0",
            "ENSG000002\t0\t3",
            "ENSG000003\t2\t1",
        ]
    ) + "\n"
    write_gzip_text(path, text)


def test_gene_by_cell_csv():
    reset_tmp()
    path = TMP / "GSM1_counts.csv.gz"
    write_gene_by_cell_csv(path)

    result = inspect_count_table(path)

    print("\n=== Gene by cell CSV ===")
    print(result)

    assert result.separator == ","
    assert result.compression == "gzip"
    assert result.first_column_role == "feature_id"
    assert result.orientation_label == "gene_by_cell"
    assert result.transpose_required_for_anndata is True
    assert result.feature_id_type == "gene_symbol_or_custom_feature"


def test_cell_by_gene_csv():
    reset_tmp()
    path = TMP / "GSM2_counts.csv.gz"
    write_cell_by_gene_csv(path)

    result = inspect_count_table(path)

    print("\n=== Cell by gene CSV ===")
    print(result)

    assert result.first_column_role == "cell_id"
    assert result.orientation_label == "cell_by_gene"
    assert result.transpose_required_for_anndata is False


def test_ensembl_txt():
    reset_tmp()
    path = TMP / "GSM3_counts.txt.gz"
    write_ensembl_txt(path)

    result = inspect_count_table(path)

    print("\n=== Ensembl TXT ===")
    print(result)

    assert result.separator == "\t"
    assert result.first_column_role == "feature_id"
    assert result.orientation_label == "gene_by_cell"
    assert result.feature_id_type == "ensembl_id"


def test_discover_and_summary():
    reset_tmp()
    write_gene_by_cell_csv(TMP / "GSM1_counts.csv.gz")
    write_cell_by_gene_csv(TMP / "GSM2_counts.csv.gz")
    write_ensembl_txt(TMP / "GSM3_counts.txt.gz")

    tables = discover_count_tables(TMP)
    inspections = inspect_many_count_tables(tables)
    df = summarize_count_table_inspections(inspections)

    print("\n=== Count table summary ===")
    print(df)

    assert len(tables) == 3
    assert len(df) == 3
    assert set(df["sample_id"]) == {"GSM1", "GSM2", "GSM3"}


def main():
    test_gene_by_cell_csv()
    test_cell_by_gene_csv()
    test_ensembl_txt()
    test_discover_and_summary()
    print("\nAll count table standardizer tests passed.")


if __name__ == "__main__":
    main()
