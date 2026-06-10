import os
import tarfile
import gzip
import argparse
from urllib.request import urlretrieve

import pandas as pd
from scipy.io import mmread


PAPER_ID = "2021_Hippen"
GEO_ACCESSION = "GSE158937"

BASE_DIR = os.path.join("data", "case_004_2021_hippen")
RAW_DIR = os.path.join(BASE_DIR, "raw")
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")

RAW_TAR_FILENAME = "GSE158937_RAW.tar"
RAW_TAR_PATH = os.path.join(RAW_DIR, RAW_TAR_FILENAME)

GEO_DOWNLOAD_URL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE158937&format=file"


SAMPLE_METADATA = {
    "GSM4816045": {
        "patient_id": "16030X2",
        "dataset_id": "GSE158937",
        "cancer_type": "Ovarian",
        "tumor_site": "Unknown",
        "metastasis_site": "Unknown",
        "tumor_treatment": "Unknown",
        "cancer_site_origin": "Ovary",
        "tumour_grade": "HGSC",
        "tumour_stage": "Unknown",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    },
    "GSM4816046": {
        "patient_id": "16030X3",
        "dataset_id": "GSE158937",
        "cancer_type": "Ovarian",
        "tumor_site": "Unknown",
        "metastasis_site": "Unknown",
        "tumor_treatment": "Unknown",
        "cancer_site_origin": "Ovary",
        "tumour_grade": "HGSC",
        "tumour_stage": "Unknown",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    },
    "GSM4816047": {
        "patient_id": "16030X4",
        "dataset_id": "GSE158937",
        "cancer_type": "Ovarian",
        "tumor_site": "Unknown",
        "metastasis_site": "Unknown",
        "tumor_treatment": "Unknown",
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
        if f.startswith("GSM") and f.endswith(".gz")
    ]

    if extracted_files:
        print(f"RAW archive already appears extracted. Found {len(extracted_files)} GSM gz files.")
        return

    print("Extracting RAW tar archive...")

    with tarfile.open(RAW_TAR_PATH, "r") as tar:
        tar.extractall(RAW_DIR)

    print("Extraction completed.")


def list_files() -> list[str]:
    files = sorted([
        f for f in os.listdir(RAW_DIR)
        if f.startswith("GSM") and os.path.isfile(os.path.join(RAW_DIR, f))
    ])

    print("\n=== Detected GSM files ===")
    print("Number of files:", len(files))

    for f in files:
        print("-", f)

    return files


def group_files_by_gsm(files: list[str]) -> dict[str, dict[str, str]]:
    grouped = {}

    for f in files:
        gsm = f.split("_")[0]
        grouped.setdefault(gsm, {})

        if "_matrix_" in f and f.endswith(".mtx.gz"):
            grouped[gsm]["matrix"] = f
        elif "_features_" in f and f.endswith(".tsv.gz"):
            grouped[gsm]["features"] = f
        elif "_barcodes_" in f and f.endswith(".tsv.gz"):
            grouped[gsm]["barcodes"] = f
        else:
            grouped[gsm].setdefault("other", [])
            grouped[gsm]["other"].append(f)

    print("\n=== Grouped files by GSM ===")

    for gsm, roles in sorted(grouped.items()):
        print(f"\n{gsm}")
        print("  matrix:", roles.get("matrix"))
        print("  features:", roles.get("features"))
        print("  barcodes:", roles.get("barcodes"))
        if roles.get("other"):
            print("  other:", roles.get("other"))

    return grouped


def inspect_first_triplet(grouped: dict[str, dict[str, str]]) -> None:
    complete_gsms = [
        gsm for gsm, roles in grouped.items()
        if {"matrix", "features", "barcodes"}.issubset(set(roles))
    ]

    if not complete_gsms:
        raise FileNotFoundError("No complete matrix/features/barcodes triplet found.")

    gsm = sorted(complete_gsms)[0]
    roles = grouped[gsm]

    matrix_path = os.path.join(RAW_DIR, roles["matrix"])
    features_path = os.path.join(RAW_DIR, roles["features"])
    barcodes_path = os.path.join(RAW_DIR, roles["barcodes"])

    print("\n=== Inspecting first 10x-like triplet ===")
    print("GSM:", gsm)
    print("Matrix:", roles["matrix"])
    print("Features:", roles["features"])
    print("Barcodes:", roles["barcodes"])

    with gzip.open(matrix_path, "rt") as f:
        matrix = mmread(f).tocsr()

    features = pd.read_csv(features_path, sep="\t", header=None, compression="gzip")
    barcodes = pd.read_csv(barcodes_path, sep="\t", header=None, compression="gzip")

    print("\nMatrix shape from Matrix Market:")
    print(matrix.shape)

    print("\nFeatures shape:")
    print(features.shape)

    print("\nBarcodes shape:")
    print(barcodes.shape)

    print("\nFirst features rows:")
    print(features.head())

    print("\nFirst barcodes rows:")
    print(barcodes.head())

    if matrix.shape[0] == features.shape[0] and matrix.shape[1] == barcodes.shape[0]:
        orientation = "genes_x_cells_before_transpose"
        needs_transpose = True
    elif matrix.shape[0] == barcodes.shape[0] and matrix.shape[1] == features.shape[0]:
        orientation = "cells_x_genes_already"
        needs_transpose = False
    else:
        orientation = "ambiguous"
        needs_transpose = None

    duplicated_gene_symbols = features.iloc[:, 1].astype(str).duplicated().sum() if features.shape[1] > 1 else None

    print("\n=== Interpretation ===")
    print("Detected format: 10x-like Matrix Market triplet")
    print("Matrix orientation:", orientation)
    print("Needs transpose for AnnData:", needs_transpose)
    print("Feature column used as gene symbols: column 1")
    print("Duplicated gene symbols in first sample:", duplicated_gene_symbols)

    print("\nCorrect AnnData construction rule:")
    print("- Read matrix.mtx.gz as sparse Matrix Market.")
    print("- Read barcodes.tsv.gz and assign real cell barcodes to obs_names.")
    print("- Prefix obs_names with GSM sample ID.")
    print("- Read features.tsv.gz and assign gene symbols to var_names.")
    print("- Transpose matrix if Matrix Market is genes x cells.")
    print("- Merge duplicated gene symbols using sparse aggregation, not dense pandas conversion.")


def cleanup_raw_files() -> None:
    print("\n=== Cleaning raw downloaded files ===")

    removed = 0

    for filename in os.listdir(RAW_DIR):
        filepath = os.path.join(RAW_DIR, filename)

        if os.path.isfile(filepath):
            os.remove(filepath)
            removed += 1

    print(f"Removed {removed} raw files from: {RAW_DIR}")


def run_case_004_inspection(keep_raw: bool = False) -> None:
    print("Case:", PAPER_ID)
    print("GEO accession:", GEO_ACCESSION)
    print("Raw directory:", RAW_DIR)
    print("Keep raw files:", keep_raw)

    ensure_directories()
    download_raw_tar()
    extract_raw_tar()

    files = list_files()
    grouped = group_files_by_gsm(files)
    inspect_first_triplet(grouped)

    print("\n=== Metadata samples ===")
    print("Curated sample metadata entries:", len(SAMPLE_METADATA))
    for gsm, meta in SAMPLE_METADATA.items():
        print(f"- {gsm}: {meta['patient_id']} / {meta['cancer_site_origin']} / {meta['tumour_grade']} / {meta['histological_subtype']}")

    print("\n=== Tool module implications ===")
    print("- src/tenx_mtx_standardizer.py should support local 10x-like triplet folders.")
    print("- src/gene_mapping_and_deduplication.py should support sparse duplicate-gene aggregation.")
    print("- src/metadata_checker.py should detect manual metadata dictionaries and missing GSM metadata.")

    if not keep_raw:
        cleanup_raw_files()

    print("\nInspection completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Inspect Case 004: 2021_Hippen / GSE158937"
    )
    parser.add_argument(
        "--keep-raw",
        action="store_true",
        help="Keep downloaded raw files after inspection."
    )

    args = parser.parse_args()
    run_case_004_inspection(keep_raw=args.keep_raw)
