import os
import argparse
from urllib.request import urlretrieve

import pandas as pd
import anndata as ad


PAPER_ID = "2022_Ignacio"
GEO_ACCESSION = "GSE180661"

BASE_DIR = os.path.join("data", "case_010_2022_ignacio")
RAW_DIR = os.path.join(BASE_DIR, "raw")
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")

METADATA_FILENAME = "GSE180661_GEO_cells.tsv.gz"
MATRIX_FILENAME = "GSE180661_matrix.h5"

METADATA_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180661/suppl/GSE180661%5FGEO%5Fcells%2Etsv%2Egz"
MATRIX_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE180nnn/GSE180661/suppl/GSE180661%5Fmatrix%2Eh5"

METADATA_PATH = os.path.join(RAW_DIR, METADATA_FILENAME)
MATRIX_PATH = os.path.join(RAW_DIR, MATRIX_FILENAME)


def ensure_directories() -> None:
    os.makedirs(RAW_DIR, exist_ok=True)
    os.makedirs(OUTPUT_DIR, exist_ok=True)


def download_file(url: str, path: str) -> None:
    if os.path.exists(path):
        print(f"File already exists: {path}")
        return

    print("Downloading:", url)
    print("Output:", path)

    urlretrieve(url, path)

    print("Download completed.")


def download_inputs() -> None:
    download_file(METADATA_URL, METADATA_PATH)
    download_file(MATRIX_URL, MATRIX_PATH)


def inspect_metadata() -> pd.DataFrame:
    print("\n=== Inspecting rich cell metadata ===")
    print("File:", METADATA_PATH)

    df_meta_preview = pd.read_csv(METADATA_PATH, sep="\t", compression="gzip", nrows=10)

    print("\nMetadata preview shape:")
    print(df_meta_preview.shape)

    print("\nMetadata columns:")
    for col in df_meta_preview.columns:
        print("-", col)

    print("\nMetadata preview rows:")
    print(df_meta_preview.head())

    required_columns = [
        "cell_id",
        "sample",
        "patient_id",
        "tumor_supersite",
        "tumor_subsite",
    ]

    print("\nRequired columns check:")
    for col in required_columns:
        print(f"- {col}:", col in df_meta_preview.columns)

    return df_meta_preview


def inspect_full_metadata_summary() -> pd.DataFrame:
    print("\n=== Inspecting full metadata summary ===")

    df_meta = pd.read_csv(METADATA_PATH, sep="\t", compression="gzip")

    print("Full metadata shape:")
    print(df_meta.shape)

    if "cell_id" in df_meta.columns:
        print("\nNumber of unique cell_id values:")
        print(df_meta["cell_id"].astype(str).nunique())

        print("\nFirst cell_id values:")
        print(df_meta["cell_id"].astype(str).head())

    for col in ["sample", "patient_id", "tumor_supersite", "tumor_subsite", "cell_type", "cell_type_super"]:
        if col in df_meta.columns:
            print(f"\nValue counts for {col}:")
            print(df_meta[col].astype(str).value_counts().head(30))

    return df_meta


def inspect_h5_matrix(df_meta: pd.DataFrame) -> None:
    print("\n=== Inspecting matrix.h5 ===")
    print("File:", MATRIX_PATH)

    try:
        adata_backed = ad.read_h5ad(MATRIX_PATH, backed="r")
        print("Read mode: anndata.read_h5ad(..., backed='r')")
        print("AnnData shape:")
        print(adata_backed.shape)

        print("\nFirst obs_names:")
        print(list(adata_backed.obs_names[:10]))

        print("\nFirst var_names:")
        print(list(adata_backed.var_names[:10]))

        print("\nExisting obs columns in matrix:")
        print(list(adata_backed.obs.columns))

        print("\nExisting var columns in matrix:")
        print(list(adata_backed.var.columns))

        n_obs = adata_backed.n_obs
        n_meta = df_meta.shape[0]

        print("\nMetadata row alignment check:")
        print("AnnData n_obs:", n_obs)
        print("Metadata rows:", n_meta)
        print("Same number of cells:", n_obs == n_meta)

        if "cell_id" in df_meta.columns:
            meta_cell_ids = df_meta["cell_id"].astype(str)

            exact_match_count = sum(
                obs_name == cell_id
                for obs_name, cell_id in zip(adata_backed.obs_names.astype(str), meta_cell_ids)
            )

            obs_name_set = set(adata_backed.obs_names.astype(str))
            meta_cell_set = set(meta_cell_ids)

            print("\nCell ID matching check:")
            print("Position-wise obs_names == metadata cell_id:", exact_match_count)
            print("Metadata cell IDs found in AnnData obs_names:", len(meta_cell_set & obs_name_set))
            print("Metadata cell IDs missing from AnnData obs_names:", len(meta_cell_set - obs_name_set))

        adata_backed.file.close()

    except Exception as exc:
        print("Could not read matrix.h5 with anndata.read_h5ad backed mode.")
        print("Error:", repr(exc))
        print("\nThis may be a non-h5ad HDF5 object. The standardizer should inspect HDF5 keys before choosing reader.")


def explain_metadata_harmonization() -> None:
    print("\n=== Clinical metadata harmonization rule ===")

    print("tumor_supersite -> tumor_site:")
    print("- Adnexa -> Primary")
    print("- Ascites -> Ascites")
    print("- everything else -> Metastasis")

    print("\ntumor_supersite -> metastasis_site:")
    print("- if tumor_site == Metastasis, keep tumor_supersite")
    print("- otherwise set empty string")

    print("\ntumor_supersite -> cancer_site_origin:")
    print("- if tumor_site == Primary, keep tumor_supersite")
    print("- otherwise set empty string")

    print("\nConstant harmonized fields:")
    print("- dataset_id: GSE180661")
    print("- cancer_type: Ovarian")
    print("- tumor_treatment: No")
    print("- tumour_grade: HGSC")
    print("- tumour_stage: Unknown")
    print("- histological_subtype: Serous")
    print("- patient_ethnicity: Unknown")


def cleanup_raw_files() -> None:
    print("\n=== Cleaning raw downloaded files ===")

    removed = 0

    for filename in os.listdir(RAW_DIR):
        filepath = os.path.join(RAW_DIR, filename)

        if os.path.isfile(filepath):
            os.remove(filepath)
            removed += 1

    print(f"Removed {removed} raw files from: {RAW_DIR}")


def run_case_010_inspection(keep_raw: bool = False) -> None:
    print("Case:", PAPER_ID)
    print("GEO accession:", GEO_ACCESSION)
    print("Raw directory:", RAW_DIR)
    print("Keep raw files:", keep_raw)

    ensure_directories()
    download_inputs()

    inspect_metadata()
    df_meta = inspect_full_metadata_summary()
    inspect_h5_matrix(df_meta)
    explain_metadata_harmonization()

    print("\n=== Tool module implications ===")
    print("- src/h5ad_inspector.py should detect whether .h5 is h5ad-like or 10x-like.")
    print("- src/external_metadata_joiner.py should attach rich cell metadata to AnnData obs.")
    print("- src/clinical_metadata_harmonizer.py should derive standardized clinical columns from tumor_supersite.")
    print("- src/metadata_column_selector.py should preserve useful metadata and drop analysis-specific columns when needed.")
    print("- src/gene_mapping_and_deduplication.py should use dataset_gene_symbol mapping for var_names.")

    if not keep_raw:
        cleanup_raw_files()

    print("\nInspection completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Inspect Case 010: 2022_Ignacio / GSE180661"
    )
    parser.add_argument(
        "--keep-raw",
        action="store_true",
        help="Keep downloaded raw files after inspection."
    )

    args = parser.parse_args()
    run_case_010_inspection(keep_raw=args.keep_raw)
