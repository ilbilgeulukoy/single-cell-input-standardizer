from pathlib import Path
import gzip
import shutil

import anndata as ad

from src.global_sparse_writer import (
    count_global_csv_shape,
    write_global_gene_by_cell_csv_to_h5ad_sparse,
)


TMP = Path("data/test_global_sparse_writer")
OUT = Path("outputs/test_global_sparse_writer")


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


def write_global_matrix(path: Path):
    text = "\n".join(
        [
            "Unnamed: 0,S1_AAACCTGAGTCAATAG,S1_AAACCTGCACAGCCCA,S2_AAACCTGCACTTCGAA",
            "GeneA,1,0,2",
            "GeneB,0,3,1",
            "GeneC,2,1,0",
            "GeneD,0,0,0",
        ]
    ) + "\n"
    write_gzip_text(path, text)


def write_metadata(path: Path):
    text = "\n".join(
        [
            "Unnamed: 0,samples,condition",
            "S1_AAACCTGAGTCAATAG,S1,Normal",
            "S1_AAACCTGCACAGCCCA,S1,Normal",
            "S2_AAACCTGCACTTCGAA,S2,Tumor",
        ]
    ) + "\n"
    write_gzip_text(path, text)


def test_count_global_csv_shape():
    reset_tmp()
    matrix_path = TMP / "global_counts.csv.gz"
    write_global_matrix(matrix_path)

    n_cells, n_genes, cell_ids, first_genes = count_global_csv_shape(matrix_path)

    print("\n=== Global CSV shape ===")
    print(n_cells, n_genes)
    print(cell_ids)
    print(first_genes)

    assert n_cells == 3
    assert n_genes == 4
    assert first_genes[:2] == ["GeneA", "GeneB"]


def test_write_sparse_h5ad_with_metadata():
    reset_tmp()
    matrix_path = TMP / "global_counts.csv.gz"
    metadata_path = TMP / "metadata.csv.gz"
    output_h5ad = OUT / "global_sparse.h5ad"

    write_global_matrix(matrix_path)
    write_metadata(metadata_path)

    result = write_global_gene_by_cell_csv_to_h5ad_sparse(
        matrix_path=matrix_path,
        metadata_path=metadata_path,
        cell_id_column="Unnamed: 0",
        sample_id_column="samples",
        keep_metadata_cells_only=True,
        output_h5ad=output_h5ad,
    )

    adata = ad.read_h5ad(output_h5ad)

    print("\n=== Sparse writer result ===")
    print(result)
    print("\n=== h5ad ===")
    print(adata)
    print(adata.obs)

    assert adata.shape == (3, 4)
    assert "sample_id" in adata.obs.columns
    assert "condition" in adata.obs.columns
    assert adata.X.nnz == 6
    assert result.validation_passed is True


def test_write_sparse_h5ad_with_curated_standard_obs():
    reset_tmp()
    matrix_path = TMP / "global_counts.csv.gz"
    metadata_path = TMP / "metadata.csv.gz"
    output_h5ad = OUT / "global_sparse_standard_obs.h5ad"

    write_global_matrix(matrix_path)
    write_metadata(metadata_path)

    curated = {
        "S1": {
            "patient_id": "PT1",
            "dataset_id": "GSE_TEST",
            "cancer_type": "Colorectal",
            "tumor_site": "Normal",
            "metastasis_site": "Unknown",
            "tumor_treatment": "Unknown",
            "cancer_site_origin": "Colon",
            "tumour_grade": "Unknown",
            "tumour_stage": "Unknown",
            "histological_subtype": "Unknown",
            "patient_ethnicity": "Unknown",
        },
        "S2": {
            "patient_id": "PT2",
            "dataset_id": "GSE_TEST",
            "cancer_type": "Colorectal",
            "tumor_site": "Tumor",
            "metastasis_site": "Unknown",
            "tumor_treatment": "Unknown",
            "cancer_site_origin": "Colon",
            "tumour_grade": "Unknown",
            "tumour_stage": "Unknown",
            "histological_subtype": "Unknown",
            "patient_ethnicity": "Unknown",
        },
    }

    result = write_global_gene_by_cell_csv_to_h5ad_sparse(
        matrix_path=matrix_path,
        metadata_path=metadata_path,
        cell_id_column="Unnamed: 0",
        sample_id_column="samples",
        keep_metadata_cells_only=True,
        output_h5ad=output_h5ad,
        curated_sample_metadata=curated,
        standardize_obs_schema=True,
    )

    adata = ad.read_h5ad(output_h5ad)

    print("\n=== Sparse writer standardized obs result ===")
    print(result)
    print(adata)
    print(adata.obs)

    expected_columns = [
        "sample_id",
        "patient_id",
        "dataset_id",
        "cancer_type",
        "tumor_site",
        "metastasis_site",
        "tumor_treatment",
        "cancer_site_origin",
        "tumour_grade",
        "tumour_stage",
        "histological_subtype",
        "patient_ethnicity",
    ]

    assert adata.shape == (3, 4)
    assert list(adata.obs.columns[: len(expected_columns)]) == expected_columns
    assert adata.obs.loc["S1_AAACCTGAGTCAATAG", "patient_id"] == "PT1"
    assert adata.obs.loc["S2_AAACCTGCACTTCGAA", "tumor_site"] == "Tumor"
    assert result.metadata_standardization_passed is True


def main():
    test_count_global_csv_shape()
    test_write_sparse_h5ad_with_metadata()
    test_write_sparse_h5ad_with_curated_standard_obs()
    print("\nAll global sparse writer tests passed.")


if __name__ == "__main__":
    main()
