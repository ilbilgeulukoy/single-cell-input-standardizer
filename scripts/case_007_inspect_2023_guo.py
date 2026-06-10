import os
import tarfile
import argparse
from urllib.request import urlretrieve

import pandas as pd


PAPER_ID = "2023_Guo"
GEO_ACCESSION = "GSE181955"

BASE_DIR = os.path.join("data", "case_007_2023_guo")
RAW_DIR = os.path.join(BASE_DIR, "raw")
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")

RAW_TAR_FILENAME = "GSE181955_RAW.tar"
RAW_TAR_PATH = os.path.join(RAW_DIR, RAW_TAR_FILENAME)

GEO_DOWNLOAD_URL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE181955&format=file"


SAMPLE_METADATA = {
    "GSM5514788": {
        "patient_id": "OMT-1_CD45_POS",
        "dataset_id": "GSE181955",
        "cancer_type": "Ovarian",
        "tumor_site": "Metastasis",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Omentum",
        "tumour_grade": "HGSC",
        "tumour_stage": "Unknown",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    },
    "GSM5514789": {
        "patient_id": "OMT-1_CD45_NEG",
        "dataset_id": "GSE181955",
        "cancer_type": "Ovarian",
        "tumor_site": "Metastasis",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Omentum",
        "tumour_grade": "HGSC",
        "tumour_stage": "Unknown",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    },
    "GSM5514790": {
        "patient_id": "OMT-3_CD45_POS",
        "dataset_id": "GSE181955",
        "cancer_type": "Ovarian",
        "tumor_site": "Metastasis",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Omentum",
        "tumour_grade": "HGSC",
        "tumour_stage": "Unknown",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    },
    "GSM5514791": {
        "patient_id": "OMT-3_CD45_NEG",
        "dataset_id": "GSE181955",
        "cancer_type": "Ovarian",
        "tumor_site": "Metastasis",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Omentum",
        "tumour_grade": "HGSC",
        "tumour_stage": "Unknown",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    },
    "GSM5514792": {
        "patient_id": "T1",
        "dataset_id": "GSE181955",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Ovary",
        "tumour_grade": "HGSC",
        "tumour_stage": "Unknown",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    },
    "GSM5514794": {
        "patient_id": "T6_CD45_NEG",
        "dataset_id": "GSE181955",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Ovary",
        "tumour_grade": "HGSC",
        "tumour_stage": "Unknown",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    },
}


def ensure_directories() -> None:
    os.makedirs(RAW_DIR, exist_ok=True)
    os.makedirs(OUTPUT_DIR, exist_ok=True)


def download_raw_tar() -> None:
    if os.path.exists(RAW_TAR_PATH):
        print(f"RAW tar already exists: {RAW_TAR_PATH}")
        return

    print("Downloading GEO supplementary archive...")
    print("URL:", GEO_DOWNLOAD_URL)
    print("Output:", RAW_TAR_PATH)

    urlretrieve(GEO_DOWNLOAD_URL, RAW_TAR_PATH)

    print("Download completed.")


def extract_raw_tar() -> None:
    extracted_files = [
        f for f in os.listdir(RAW_DIR)
        if f.startswith("GSM") and f.endswith(".csv.gz")
    ]

    if extracted_files:
        print(f"RAW archive already appears extracted. Found {len(extracted_files)} GSM csv.gz files.")
        return

    print("Extracting RAW tar archive...")

    with tarfile.open(RAW_TAR_PATH, "r") as tar:
        tar.extractall(RAW_DIR)

    print("Extraction completed.")


def list_csv_files() -> list[str]:
    files = sorted([
        f for f in os.listdir(RAW_DIR)
        if f.startswith("GSM")
        and f.endswith(".csv.gz")
        and os.path.isfile(os.path.join(RAW_DIR, f))
    ])

    selected_files = [
        f for f in files
        if f.split("_")[0] in SAMPLE_METADATA
    ]

    skipped_files = [
        f for f in files
        if f.split("_")[0] not in SAMPLE_METADATA
    ]

    print("\n=== Detected GSM CSV files ===")
    print("Number of all GSM csv.gz files:", len(files))

    for f in files:
        print("-", f)

    print("\n=== Selected CSV files with curated metadata ===")
    print("Number of selected files:", len(selected_files))

    for f in selected_files:
        print("-", f)

    print("\n=== Skipped CSV files without curated metadata ===")
    print("Number of skipped files:", len(skipped_files))

    for f in skipped_files:
        print("-", f)

    return [os.path.join(RAW_DIR, f) for f in selected_files]


def infer_feature_id_type(values: pd.Series) -> str:
    as_str = values.dropna().astype(str)

    if as_str.empty:
        return "unknown"

    ensembl_like = as_str.str.startswith("ENSG").mean()
    symbol_like = as_str.str.match(r"^[A-Za-z][A-Za-z0-9_.-]*$").mean()

    if ensembl_like > 0.8:
        return "ensembl_gene_id"
    if symbol_like > 0.8:
        return "gene_symbol"
    return "ambiguous"


def inspect_first_csv(selected_paths: list[str]) -> None:
    if not selected_paths:
        raise FileNotFoundError("No selected GSM csv.gz files found.")

    first_file = selected_paths[0]

    print("\n=== Inspecting first selected CSV file ===")
    print("File:", first_file)

    df_preview = pd.read_csv(first_file, compression="gzip", nrows=5)

    print("\nPreview shape:")
    print(df_preview.shape)

    print("\nFirst columns:")
    print(list(df_preview.columns[:10]))

    print("\nPreview rows:")
    print(df_preview.head())

    first_col = df_preview.columns[0]
    feature_id_type = infer_feature_id_type(df_preview[first_col])

    print("\n=== Interpretation ===")
    print("Detected format: compressed CSV count table")
    print("First column name:", first_col)
    print("Inferred feature identifier type from preview:", feature_id_type)
    print("Remaining columns are treated as cell barcodes or cell IDs.")
    print("Raw orientation: genes/features x cells")
    print("Needs transpose for AnnData: True")

    print("\nCorrect AnnData construction rule:")
    print("- Read csv.gz with pandas.")
    print("- Rename first column to gene or feature_id.")
    print("- Set first column as index.")
    print("- Build AnnData with AnnData(df.T).")
    print("- Assign obs_names from original cell columns.")
    print("- Assign var_names from first column.")
    print("- Detect whether var_names are gene symbols or Ensembl IDs.")
    print("- If Ensembl-like, gene mapping should use dataset_gene_ensembl_id.")
    print("- If symbol-like, gene mapping should use dataset_gene_symbol.")


def cleanup_raw_files() -> None:
    print("\n=== Cleaning raw downloaded files ===")

    removed = 0

    for filename in os.listdir(RAW_DIR):
        filepath = os.path.join(RAW_DIR, filename)

        if os.path.isfile(filepath):
            os.remove(filepath)
            removed += 1

    print(f"Removed {removed} raw files from: {RAW_DIR}")


def run_case_007_inspection(keep_raw: bool = False) -> None:
    print("Case:", PAPER_ID)
    print("GEO accession:", GEO_ACCESSION)
    print("Raw directory:", RAW_DIR)
    print("Keep raw files:", keep_raw)

    ensure_directories()
    download_raw_tar()
    extract_raw_tar()

    selected_paths = list_csv_files()
    inspect_first_csv(selected_paths)

    print("\n=== Metadata selected samples ===")
    print("Curated sample metadata entries:", len(SAMPLE_METADATA))
    for gsm, meta in SAMPLE_METADATA.items():
        print(f"- {gsm}: {meta['patient_id']} / {meta['tumor_site']} / {meta['cancer_site_origin']} / {meta['tumour_grade']}")

    print("\n=== Tool module implications ===")
    print("- src/count_table_standardizer.py should support CSV count tables with Ensembl-like feature IDs.")
    print("- src/feature_identifier_type_detector.py is needed to distinguish gene symbols from Ensembl IDs.")
    print("- src/gene_mapping_and_deduplication.py should choose mapping key based on detected feature ID type.")
    print("- src/metadata_checker.py should report skipped GSM files without curated metadata.")

    if not keep_raw:
        cleanup_raw_files()

    print("\nInspection completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Inspect Case 007: 2023_Guo / GSE181955"
    )
    parser.add_argument(
        "--keep-raw",
        action="store_true",
        help="Keep downloaded raw files after inspection."
    )

    args = parser.parse_args()
    run_case_007_inspection(keep_raw=args.keep_raw)
