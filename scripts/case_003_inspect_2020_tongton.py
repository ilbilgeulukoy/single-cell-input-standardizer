import os
import tarfile
import argparse
from urllib.request import urlretrieve

import pandas as pd


PAPER_ID = "2020_Tongton"
GEO_ACCESSION = "GSE130000"

BASE_DIR = os.path.join("data", "case_003_2020_tongton")
RAW_DIR = os.path.join(BASE_DIR, "raw")
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")

RAW_TAR_FILENAME = "GSE130000_RAW.tar"
RAW_TAR_PATH = os.path.join(RAW_DIR, RAW_TAR_FILENAME)

GEO_DOWNLOAD_URL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE130000&format=file"


SAMPLE_METADATA = {
    "GSM3729170": {
        "patient_id": "P1",
        "dataset_id": "GSE130000",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "HGSC",
        "tumour_stage": "Unknown",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    },
    "GSM3729171": {
        "patient_id": "P2",
        "dataset_id": "GSE130000",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "LGSC",
        "tumour_stage": "Unknown",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    },
    "GSM3729172": {
        "patient_id": "P3",
        "dataset_id": "GSE130000",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "LGSC",
        "tumour_stage": "Unknown",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    },
    "GSM3729173": {
        "patient_id": "P4",
        "dataset_id": "GSE130000",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "Unknown",
        "tumour_stage": "Unknown",
        "histological_subtype": "Endometrioid",
        "patient_ethnicity": "Unknown",
    },
    "GSM3729174": {
        "patient_id": "M1",
        "dataset_id": "GSE130000",
        "cancer_type": "Ovarian",
        "tumor_site": "Metastasis",
        "metastasis_site": "Peritoneal",
        "tumor_treatment": "No",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "HGSC",
        "tumour_stage": "Unknown",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    },
    "GSM3729175": {
        "patient_id": "M2",
        "dataset_id": "GSE130000",
        "cancer_type": "Ovarian",
        "tumor_site": "Metastasis",
        "metastasis_site": "Peritoneal",
        "tumor_treatment": "No",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "LGSC",
        "tumour_stage": "Unknown",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    },
    "GSM3729176": {
        "patient_id": "R1",
        "dataset_id": "GSE130000",
        "cancer_type": "Ovarian",
        "tumor_site": "Relapse",
        "metastasis_site": "Unknown",
        "tumor_treatment": "Yes",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "Unknown",
        "tumour_stage": "Unknown",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    },
    "GSM3729177": {
        "patient_id": "R2",
        "dataset_id": "GSE130000",
        "cancer_type": "Ovarian",
        "tumor_site": "Relapse",
        "metastasis_site": "Unknown",
        "tumor_treatment": "Yes",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "LGSC",
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
        if f.startswith("GSM") and f.endswith(".txt.gz")
    ]

    if extracted_files:
        print(f"RAW archive already appears extracted. Found {len(extracted_files)} GSM txt.gz files.")
        return

    print("Extracting RAW tar archive...")

    with tarfile.open(RAW_TAR_PATH, "r") as tar:
        tar.extractall(RAW_DIR)

    print("Extraction completed.")


def list_gsm_txt_files() -> list[str]:
    sample_files = sorted([
        os.path.join(RAW_DIR, f)
        for f in os.listdir(RAW_DIR)
        if f.startswith("GSM") and f.endswith(".txt.gz")
    ])

    print("\n=== Detected GSM TXT files ===")
    print("Number of files:", len(sample_files))

    for fp in sample_files:
        print("-", os.path.basename(fp))

    return sample_files


def inspect_first_txt(sample_files: list[str]) -> None:
    if not sample_files:
        raise FileNotFoundError("No GSM txt.gz files found.")

    first_file = sample_files[0]

    print("\n=== Inspecting first TXT file ===")
    print("File:", first_file)

    df_preview = pd.read_csv(first_file, sep="\t", compression="gzip", nrows=5)

    print("\nPreview shape:")
    print(df_preview.shape)

    print("\nFirst columns:")
    print(list(df_preview.columns[:10]))

    print("\nPreview rows:")
    print(df_preview.head())

    print("\nInterpretation:")
    print("- This case contains compressed TXT count tables.")
    print("- The previous working workflow used sc.read_text(...).T.")
    print("- Therefore the raw table is treated as genes x cells.")
    print("- AnnData should be cells x genes after transposition.")
    print("- This case additionally requires gene mapping, gene cleanup, and duplicated gene aggregation.")


def cleanup_raw_files() -> None:
    print("\n=== Cleaning raw downloaded files ===")

    removed = 0

    for filename in os.listdir(RAW_DIR):
        filepath = os.path.join(RAW_DIR, filename)

        if os.path.isfile(filepath):
            os.remove(filepath)
            removed += 1

    print(f"Removed {removed} raw files from: {RAW_DIR}")


def run_case_003_inspection(keep_raw: bool = False) -> None:
    print("Case:", PAPER_ID)
    print("GEO accession:", GEO_ACCESSION)
    print("Raw directory:", RAW_DIR)
    print("Keep raw files:", keep_raw)

    ensure_directories()
    download_raw_tar()
    extract_raw_tar()

    sample_files = list_gsm_txt_files()
    inspect_first_txt(sample_files)

    print("\n=== Metadata samples ===")
    print("Curated sample metadata entries:", len(SAMPLE_METADATA))
    for gsm, meta in SAMPLE_METADATA.items():
        print(f"- {gsm}: {meta['patient_id']} / {meta['tumor_site']} / {meta['tumour_grade']} / {meta['histological_subtype']}")

    print("\n=== Tool module implications ===")
    print("- src/count_table_standardizer.py should support TXT/TSV gene-by-cell tables.")
    print("- src/gene_mapping_and_deduplication.py is needed for gene ID cleanup and duplicated symbol aggregation.")

    if not keep_raw:
        cleanup_raw_files()

    print("\nInspection completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Inspect Case 003: 2020_Tongton / GSE130000"
    )
    parser.add_argument(
        "--keep-raw",
        action="store_true",
        help="Keep downloaded raw files after inspection."
    )

    args = parser.parse_args()
    run_case_003_inspection(keep_raw=args.keep_raw)
