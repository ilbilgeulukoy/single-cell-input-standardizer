import os
import tarfile
import gzip
import argparse
from urllib.request import urlretrieve

import pandas as pd
from scipy.io import mmread


PAPER_ID = "2023_Hippen"
GEO_ACCESSION = "GSE217517"

BASE_DIR = os.path.join("data", "case_011_2023_hippen")
RAW_DIR = os.path.join(BASE_DIR, "raw")
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")

RAW_TAR_FILENAME = "GSE217517_RAW.tar"
RAW_TAR_PATH = os.path.join(RAW_DIR, RAW_TAR_FILENAME)

GEO_DOWNLOAD_URL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE217517&format=file"


SAMPLE_METADATA = {
    "GSM6720925": {"patient_id": "2251", "dataset_id": "GSE217517"},
    "GSM6720926": {"patient_id": "2267", "dataset_id": "GSE217517"},
    "GSM6720927": {"patient_id": "2283", "dataset_id": "GSE217517"},
    "GSM6720928": {"patient_id": "2293", "dataset_id": "GSE217517"},
    "GSM6720929": {"patient_id": "2380", "dataset_id": "GSE217517"},
    "GSM6720930": {"patient_id": "2428", "dataset_id": "GSE217517"},
    "GSM6720931": {"patient_id": "2467", "dataset_id": "GSE217517"},
    "GSM6720932": {"patient_id": "2497", "dataset_id": "GSE217517"},
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
        if f.startswith("GSM") and f.endswith((".gz", ".mtx.gz", ".tsv.gz"))
    ]

    if extracted_files:
        print(f"RAW archive already appears extracted. Found {len(extracted_files)} GSM files.")
        return

    print("Extracting RAW tar archive...")

    with tarfile.open(RAW_TAR_PATH, "r") as tar:
        tar.extractall(RAW_DIR)

    print("Extraction completed.")


def list_files() -> list[str]:
    files = sorted([
        f for f in os.listdir(RAW_DIR)
        if f.startswith("GSM")
        and os.path.isfile(os.path.join(RAW_DIR, f))
    ])

    print("\n=== Detected GSM files ===")
    print("Number of GSM files:", len(files))

    for f in files:
        print("-", f)

    return files


def classify_files(files: list[str]) -> dict[str, list[str]]:
    matrix_files = [
        f for f in files
        if "single_cell" in f and "matrix" in f and f.endswith(".mtx.gz")
    ]
    feature_files = [
        f for f in files
        if "single_cell" in f and "features" in f and f.endswith(".tsv.gz")
    ]
    barcode_files = [
        f for f in files
        if "single_cell" in f and "barcodes" in f and f.endswith(".tsv.gz")
    ]

    bulk_or_other_files = [
        f for f in files
        if f not in matrix_files + feature_files + barcode_files
    ]

    print("\n=== Classified files ===")
    print("Single-cell matrix files:", len(matrix_files))
    for f in matrix_files:
        print("-", f)

    print("\nSingle-cell feature files:", len(feature_files))
    for f in feature_files:
        print("-", f)

    print("\nSingle-cell barcode files:", len(barcode_files))
    for f in barcode_files:
        print("-", f)

    print("\nOther / possible bulk files:", len(bulk_or_other_files))
    for f in bulk_or_other_files[:50]:
        print("-", f)

    return {
        "matrix": matrix_files,
        "features": feature_files,
        "barcodes": barcode_files,
        "other": bulk_or_other_files,
    }


def detect_triplet_mode(classified: dict[str, list[str]]) -> str:
    n_matrix = len(classified["matrix"])
    n_features = len(classified["features"])
    n_barcodes = len(classified["barcodes"])

    if n_matrix == 1 and n_features == 1 and n_barcodes == 1:
        return "single_pooled_10x_triplet"
    if n_matrix > 1 and n_features > 1 and n_barcodes > 1:
        return "per_sample_10x_triplets"
    return "ambiguous_or_incomplete"


def inspect_first_triplet(classified: dict[str, list[str]]) -> None:
    if not classified["matrix"] or not classified["features"] or not classified["barcodes"]:
        print("\nNo complete single-cell matrix/features/barcodes triplet detected.")
        return

    matrix_file = classified["matrix"][0]
    features_file = classified["features"][0]
    barcodes_file = classified["barcodes"][0]

    matrix_path = os.path.join(RAW_DIR, matrix_file)
    features_path = os.path.join(RAW_DIR, features_file)
    barcodes_path = os.path.join(RAW_DIR, barcodes_file)

    print("\n=== Inspecting first / pooled single-cell triplet ===")
    print("Matrix:", matrix_file)
    print("Features:", features_file)
    print("Barcodes:", barcodes_file)

    with gzip.open(matrix_path, "rt") as f:
        matrix = mmread(f).tocsr()

    features = pd.read_csv(features_path, sep="\t", header=None, compression="gzip")
    barcodes = pd.read_csv(barcodes_path, sep="\t", header=None, compression="gzip")

    print("\nMatrix shape:")
    print(matrix.shape)

    print("\nFeatures shape:")
    print(features.shape)

    print("\nBarcodes shape:")
    print(barcodes.shape)

    print("\nFirst feature rows:")
    print(features.head())

    print("\nFirst barcode rows:")
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

    gene_symbol_col = 1 if features.shape[1] > 1 else 0
    duplicated_gene_symbols = features.iloc[:, gene_symbol_col].astype(str).duplicated().sum()

    print("\n=== Interpretation ===")
    print("Detected format: 10x-like Matrix Market single-cell triplet")
    print("Triplet mode:", detect_triplet_mode(classified))
    print("Matrix orientation:", orientation)
    print("Needs transpose for AnnData:", needs_transpose)
    print("Gene symbol column:", gene_symbol_col)
    print("Duplicated gene symbols:", duplicated_gene_symbols)

    print("\nImportant script correction:")
    print("- Do not loop over GSM IDs if only one pooled single-cell triplet exists.")
    print("- Do not assign the same pooled triplet repeatedly to each GSM.")
    print("- If barcodes encode sample/patient identity, parse metadata from barcode or external annotation.")
    print("- If barcodes do not encode sample identity, this case may require pooled-level metadata only.")


def inspect_gsm_metadata_vs_files(files: list[str]) -> None:
    print("\n=== GSM metadata vs detected files ===")

    detected_gsms = sorted({f.split("_")[0] for f in files if f.startswith("GSM")})
    metadata_gsms = sorted(SAMPLE_METADATA)

    print("GSMs detected in file names:")
    print(detected_gsms)

    print("\nGSMs in curated metadata dictionary:")
    print(metadata_gsms)

    print("\nDetected GSMs missing from metadata:")
    print(sorted(set(detected_gsms) - set(metadata_gsms)))

    print("\nMetadata GSMs missing from detected files:")
    print(sorted(set(metadata_gsms) - set(detected_gsms)))


def cleanup_raw_files() -> None:
    print("\n=== Cleaning raw downloaded files ===")

    removed = 0

    for filename in os.listdir(RAW_DIR):
        filepath = os.path.join(RAW_DIR, filename)

        if os.path.isfile(filepath):
            os.remove(filepath)
            removed += 1

    print(f"Removed {removed} raw files from: {RAW_DIR}")


def run_case_011_inspection(keep_raw: bool = False) -> None:
    print("Case:", PAPER_ID)
    print("GEO accession:", GEO_ACCESSION)
    print("Raw directory:", RAW_DIR)
    print("Keep raw files:", keep_raw)

    ensure_directories()
    download_raw_tar()
    extract_raw_tar()

    files = list_files()
    classified = classify_files(files)
    inspect_gsm_metadata_vs_files(files)

    mode = detect_triplet_mode(classified)

    print("\n=== Triplet mode decision ===")
    print("Triplet mode:", mode)

    inspect_first_triplet(classified)

    print("\n=== Tool module implications ===")
    print("- src/tenx_mtx_standardizer.py should distinguish pooled triplet vs per-sample triplets.")
    print("- src/pooled_sample_resolver.py may be needed if barcodes encode sample identity.")
    print("- src/metadata_checker.py should detect when curated GSM metadata cannot be assigned cell-by-cell.")
    print("- src/gene_mapping_and_deduplication.py should be applied after AnnData construction.")

    if not keep_raw:
        cleanup_raw_files()

    print("\nInspection completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Inspect Case 011: 2023_Hippen / GSE217517"
    )
    parser.add_argument(
        "--keep-raw",
        action="store_true",
        help="Keep downloaded raw files after inspection."
    )

    args = parser.parse_args()
    run_case_011_inspection(keep_raw=args.keep_raw)
