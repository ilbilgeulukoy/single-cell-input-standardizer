import os
import tarfile
import gzip
import argparse
from urllib.request import urlretrieve

import pandas as pd
from scipy.io import mmread


PAPER_ID = "2021_Li"
GEO_ACCESSION = "GSE168652"

BASE_DIR = os.path.join("data", "case_013_2021_li")
RAW_DIR = os.path.join(BASE_DIR, "raw")

RAW_TAR_FILENAME = "GSE168652_RAW.tar"
RAW_TAR_PATH = os.path.join(RAW_DIR, RAW_TAR_FILENAME)

GEO_DOWNLOAD_URL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE168652&format=file"


CURATED_METADATA_FROM_ORIGINAL_SCRIPT = {
    "GSM5155196": {
        "patient_id": "N_OT_PT1",
        "dataset_id": "GSE168652",
        "cancer_type": "Ovarian",
        "tumor_site": "Metastasis",
        "metastasis_site": "Omentum",
        "tumor_treatment": "No",
        "cancer_site_origin": "Omentum",
        "tumour_grade": "SCC",
        "tumour_stage": "IIIc",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    }
}


EXPECTED_GEO_CONTEXT = {
    "disease": "cervical cancer",
    "sample": "cervical cancer tissue",
    "normal_sample_present": True,
    "expected_histology": "cervical squamous cell carcinoma",
    "hpv_status": "HPV18-positive",
}


def ensure_directories():
    os.makedirs(RAW_DIR, exist_ok=True)


def download_raw_tar():
    if os.path.exists(RAW_TAR_PATH):
        print(f"RAW tar already exists: {RAW_TAR_PATH}")
        return

    print("Downloading GEO supplementary archive...")
    print("URL:", GEO_DOWNLOAD_URL)
    print("Output:", RAW_TAR_PATH)

    urlretrieve(GEO_DOWNLOAD_URL, RAW_TAR_PATH)

    print("Download completed.")


def extract_raw_tar():
    extracted_files = [
        f for f in os.listdir(RAW_DIR)
        if f.startswith("GSM")
    ]

    if extracted_files:
        print(f"RAW archive already appears extracted. Found {len(extracted_files)} GSM files.")
        return

    print("Extracting GEO RAW tar archive...")

    with tarfile.open(RAW_TAR_PATH, "r") as tar:
        tar.extractall(RAW_DIR)

    print("Extraction completed.")


def list_raw_files():
    files = sorted(os.listdir(RAW_DIR))

    print("\n=== Detected top-level raw files ===")
    print("Number of files/directories:", len(files))

    for f in files:
        print("-", f)

    return files


def inspect_nested_tar(top_files):
    nested_tars = [
        f for f in top_files
        if f.startswith("GSM") and f.endswith(".tar.gz")
    ]

    print("\n=== Nested tar files ===")
    print("Nested tar.gz count:", len(nested_tars))

    for f in nested_tars:
        print("-", f)

    if not nested_tars:
        print("No nested tar.gz files detected.")
        return None

    tumor_tar = None
    for f in nested_tars:
        if "Tumor" in f:
            tumor_tar = f
            break

    if tumor_tar is None:
        tumor_tar = nested_tars[0]

    tumor_tar_path = os.path.join(RAW_DIR, tumor_tar)

    print("\nSelected nested tar for inspection:", tumor_tar)

    with tarfile.open(tumor_tar_path, "r:gz") as tar:
        members = tar.getnames()

        print("\nNested tar members:")
        for member in members:
            print("-", member)

        tar.extractall(RAW_DIR)

    print("\nNested tar extracted.")

    return tumor_tar


def find_tumor_triplet():
    tumor_dir = os.path.join(RAW_DIR, "Tumor")

    print("\n=== Tumor directory check ===")
    print("Tumor dir exists:", os.path.isdir(tumor_dir))

    if not os.path.isdir(tumor_dir):
        return None

    tumor_files = sorted(os.listdir(tumor_dir))

    print("Tumor directory files:")
    for f in tumor_files:
        print("-", f)

    matrix_files = [f for f in tumor_files if "matrix" in f and f.endswith(".mtx.gz")]
    feature_files = [f for f in tumor_files if "features" in f and f.endswith(".tsv.gz")]
    barcode_files = [f for f in tumor_files if "barcodes" in f and f.endswith(".tsv.gz")]

    print("\nMatrix files:", matrix_files)
    print("Feature files:", feature_files)
    print("Barcode files:", barcode_files)

    if len(matrix_files) != 1 or len(feature_files) != 1 or len(barcode_files) != 1:
        print("WARNING: expected exactly one matrix/features/barcodes triplet.")
        return None

    return {
        "tumor_dir": tumor_dir,
        "matrix": os.path.join(tumor_dir, matrix_files[0]),
        "features": os.path.join(tumor_dir, feature_files[0]),
        "barcodes": os.path.join(tumor_dir, barcode_files[0]),
    }


def inspect_triplet(paths):
    if not paths:
        print("\nNo complete triplet to inspect.")
        return

    print("\n=== Inspecting Tumor 10x Matrix Market triplet ===")

    print("Matrix:", paths["matrix"])
    print("Features:", paths["features"])
    print("Barcodes:", paths["barcodes"])

    with gzip.open(paths["matrix"], "rt") as f:
        matrix = mmread(f).tocsr()

    features = pd.read_csv(paths["features"], sep="\t", header=None, compression="gzip")
    barcodes = pd.read_csv(paths["barcodes"], sep="\t", header=None, compression="gzip")

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
    print("Detected format: nested tar containing single-sample 10x Matrix Market triplet")
    print("Selected sample: GSM5155196 Tumor")
    print("Matrix orientation:", orientation)
    print("Needs transpose for AnnData:", needs_transpose)
    print("Gene symbol column:", gene_symbol_col)
    print("Duplicated gene symbols:", duplicated_gene_symbols)


def inspect_metadata_sanity():
    print("\n=== Metadata sanity check ===")

    print("GEO context expected from accession:")
    for key, value in EXPECTED_GEO_CONTEXT.items():
        print(f"- {key}: {value}")

    print("\nOriginal script curated metadata:")
    for gsm, meta in CURATED_METADATA_FROM_ORIGINAL_SCRIPT.items():
        print("GSM:", gsm)
        for key, value in meta.items():
            print(f"  {key}: {value}")

    meta = CURATED_METADATA_FROM_ORIGINAL_SCRIPT["GSM5155196"]

    warnings = []

    if meta["cancer_type"].lower() != "cervical cancer" and meta["cancer_type"].lower() != "cervical":
        warnings.append("cancer_type in script is not cervical cancer")

    if meta["histological_subtype"].lower() == "serous":
        warnings.append("histological_subtype 'Serous' conflicts with cervical squamous cell carcinoma context")

    if meta["tumor_site"].lower() == "metastasis" or meta["metastasis_site"].lower() == "omentum":
        warnings.append("tumor/metastasis site looks ovarian-style and conflicts with GEO cervical tumor context")

    print("\nSanity warnings:")
    if warnings:
        for warning in warnings:
            print("-", warning)
    else:
        print("No obvious metadata conflicts detected.")

    print("\nRecommended policy:")
    print("- Keep this case as a metadata sanity-check example.")
    print("- Do not silently reuse ovarian-specific metadata defaults for cervical cancer datasets.")
    print("- If project scope is ovarian-only, flag this dataset as out-of-scope or non-ovarian.")
    print("- If kept for format testing, mark biological metadata as intentionally invalid / needs correction.")


def cleanup_raw_files():
    print("\n=== Cleaning raw downloaded files ===")

    removed_files = 0
    removed_dirs = 0

    for filename in os.listdir(RAW_DIR):
        filepath = os.path.join(RAW_DIR, filename)

        if os.path.isfile(filepath):
            os.remove(filepath)
            removed_files += 1
        elif os.path.isdir(filepath):
            for root, dirs, files in os.walk(filepath, topdown=False):
                for file in files:
                    os.remove(os.path.join(root, file))
                    removed_files += 1
                for directory in dirs:
                    os.rmdir(os.path.join(root, directory))
            os.rmdir(filepath)
            removed_dirs += 1

    print(f"Removed {removed_files} raw files and {removed_dirs} directories from: {RAW_DIR}")


def run_case_013_inspection(keep_raw=False):
    print("Case:", PAPER_ID)
    print("GEO accession:", GEO_ACCESSION)
    print("Raw directory:", RAW_DIR)
    print("Keep raw files:", keep_raw)
    print("Mode: scanpy-free nested tar + Matrix Market inspection")

    ensure_directories()
    download_raw_tar()
    extract_raw_tar()

    top_files = list_raw_files()
    inspect_nested_tar(top_files)
    triplet_paths = find_tumor_triplet()
    inspect_triplet(triplet_paths)
    inspect_metadata_sanity()

    print("\n=== Tool module implications ===")
    print("- src/nested_archive_extractor.py should support GSM-level nested tar.gz files.")
    print("- src/tenx_mtx_standardizer.py should process single-sample Matrix Market triplets.")
    print("- src/metadata_sanity_checker.py should compare curated metadata against GEO/project context.")
    print("- src/project_scope_filter.py may flag this dataset as non-ovarian if the project is ovarian-only.")

    if not keep_raw:
        cleanup_raw_files()

    print("\nInspection completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Inspect Case 013: 2021_Li / GSE168652"
    )
    parser.add_argument(
        "--keep-raw",
        action="store_true",
        help="Keep downloaded raw files after inspection."
    )

    args = parser.parse_args()
    run_case_013_inspection(keep_raw=args.keep_raw)
