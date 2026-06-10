import os
import tarfile
import argparse
from urllib.request import urlretrieve

import h5py
import pandas as pd


PAPER_ID = "2024_Brand"
GEO_ACCESSION = "GSE233615"

BASE_DIR = os.path.join("data", "case_012_2024_brand")
RAW_DIR = os.path.join(BASE_DIR, "raw")

RAW_TAR_FILENAME = "GSE233615_RAW.tar"
RAW_TAR_PATH = os.path.join(RAW_DIR, RAW_TAR_FILENAME)

GEO_DOWNLOAD_URL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE233615&format=file"


SAMPLE_METADATA = {
    "GSM7431434": {"patient_id": "E1", "tumor_site": "Metastasis", "metastasis_site": "Omentum", "tumour_grade": "LGSC", "patient_ethnicity": "Non-White"},
    "GSM7431435": {"patient_id": "E2", "tumor_site": "Metastasis", "metastasis_site": "Omentum", "tumour_grade": "LGSC", "patient_ethnicity": "Non-White"},
    "GSM7431436": {"patient_id": "E3", "tumor_site": "Unknown", "metastasis_site": "Left ovary", "tumour_grade": "LGSC", "patient_ethnicity": "Non-White"},
    "GSM7431437": {"patient_id": "E4", "tumor_site": "Unknown", "metastasis_site": "Right ovary", "tumour_grade": "LGSC", "patient_ethnicity": "Non-White"},
    "GSM7431438": {"patient_id": "E7", "tumor_site": "Metastasis", "metastasis_site": "Omentum", "tumour_grade": "HGSC", "patient_ethnicity": "White"},
    "GSM7431439": {"patient_id": "E9", "tumor_site": "Metastasis", "metastasis_site": "Omentum", "tumour_grade": "HGSC", "patient_ethnicity": "White"},
    "GSM7431440": {"patient_id": "E10", "tumor_site": "Metastasis", "metastasis_site": "Omentum", "tumour_grade": "HGSC", "patient_ethnicity": "White"},
    "GSM7431441": {"patient_id": "E11", "tumor_site": "Metastasis", "metastasis_site": "Abdominal Wall", "tumour_grade": "HGSC", "patient_ethnicity": "White"},
    "GSM7431442": {"patient_id": "E12", "tumor_site": "Metastasis", "metastasis_site": "Omentum", "tumour_grade": "HGSC", "patient_ethnicity": "White"},
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

    print("Extracting RAW tar archive...")

    with tarfile.open(RAW_TAR_PATH, "r") as tar:
        tar.extractall(RAW_DIR)

    print("Extraction completed.")


def decode_array(arr, n=5):
    values = arr[:n]
    decoded = []
    for x in values:
        if isinstance(x, bytes):
            decoded.append(x.decode("utf-8"))
        else:
            decoded.append(str(x))
    return decoded


def list_gsm_files():
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


def classify_files(files):
    h5_files = [f for f in files if f.endswith(".h5")]
    other_files = [f for f in files if not f.endswith(".h5")]

    print("\n=== Classified files ===")
    print("10x h5 files:", len(h5_files))
    for f in h5_files:
        print("-", f)

    print("\nOther files:", len(other_files))
    for f in other_files:
        print("-", f)

    return h5_files, other_files


def inspect_metadata_vs_files(files):
    detected_gsms = sorted({f.split("_")[0].strip() for f in files if f.startswith("GSM")})
    metadata_gsms = sorted(SAMPLE_METADATA)

    print("\n=== GSM metadata vs detected files ===")
    print("GSMs detected in file names:")
    print(detected_gsms)

    print("\nGSMs in curated metadata dictionary:")
    print(metadata_gsms)

    print("\nDetected GSMs missing from metadata:")
    print(sorted(set(detected_gsms) - set(metadata_gsms)))

    print("\nMetadata GSMs missing from detected files:")
    print(sorted(set(metadata_gsms) - set(detected_gsms)))


def inspect_10x_h5(filepath):
    filename = os.path.basename(filepath)
    sample_id = filename.split("_")[0].strip()

    print("\n=== Inspecting 10x h5 file ===")
    print("File:", filename)
    print("Sample ID:", sample_id)

    with h5py.File(filepath, "r") as f:
        print("Top-level keys:", list(f.keys()))

        if "matrix" not in f:
            print("WARNING: /matrix group not found. This may not be standard 10x h5.")
            return None

        matrix = f["matrix"]
        print("matrix subkeys:", list(matrix.keys()))

        shape = tuple(matrix["shape"][:])
        print("10x matrix shape from /matrix/shape:", shape)

        n_genes = int(shape[0])
        n_cells = int(shape[1])

        print("Interpreted raw orientation: genes x cells")
        print("AnnData expected shape after read_10x_h5:", (n_cells, n_genes))

        print("data shape:", matrix["data"].shape)
        print("indices shape:", matrix["indices"].shape)
        print("indptr shape:", matrix["indptr"].shape)
        print("barcodes shape:", matrix["barcodes"].shape)

        print("First barcodes:", decode_array(matrix["barcodes"], 5))

        if "features" in matrix:
            features = matrix["features"]
            print("features subkeys:", list(features.keys()))

            if "id" in features:
                print("First feature IDs:", decode_array(features["id"], 5))

            if "name" in features:
                print("First feature names:", decode_array(features["name"], 5))

            if "feature_type" in features:
                print("First feature types:", decode_array(features["feature_type"], 5))

            duplicated_names = None
            if "name" in features:
                names = [
                    x.decode("utf-8") if isinstance(x, bytes) else str(x)
                    for x in features["name"][:]
                ]
                duplicated_names = pd.Index(names).duplicated().sum()
                print("Duplicated feature names:", duplicated_names)

        patient_from_filename = None
        parts = filename.split("_")
        if len(parts) > 1:
            patient_from_filename = parts[1].replace(".h5", "").strip()

        print("\nFile-derived patient token:", patient_from_filename)

        if sample_id in SAMPLE_METADATA:
            print("Curated patient_id:", SAMPLE_METADATA[sample_id]["patient_id"])
            print("Patient token matches curated:", patient_from_filename == SAMPLE_METADATA[sample_id]["patient_id"])
        else:
            print("Sample ID not in curated metadata.")

        return {
            "sample_id": sample_id,
            "filename": filename,
            "n_cells": n_cells,
            "n_genes": n_genes,
            "patient_token": patient_from_filename,
            "duplicated_names": duplicated_names,
        }


def inspect_all_h5_files(h5_files):
    print("\n=== All h5 sample summaries ===")

    summaries = []
    for filename in h5_files:
        filepath = os.path.join(RAW_DIR, filename)
        summary = inspect_10x_h5(filepath)
        if summary:
            summaries.append(summary)

    print("\n=== Compact sample shape table ===")
    total_cells = 0

    for s in summaries:
        total_cells += s["n_cells"]
        print(
            s["sample_id"],
            s["filename"],
            "cells=", s["n_cells"],
            "genes=", s["n_genes"],
            "patient_token=", s["patient_token"],
            "dup_gene_names=", s["duplicated_names"],
        )

    print("\nTotal cells across h5 files:", total_cells)
    print("Unique gene counts:", sorted(set(s["n_genes"] for s in summaries)))
    print("Number of h5 files inspected:", len(summaries))


def cleanup_raw_files():
    print("\n=== Cleaning raw downloaded files ===")

    removed = 0

    for filename in os.listdir(RAW_DIR):
        filepath = os.path.join(RAW_DIR, filename)

        if os.path.isfile(filepath):
            os.remove(filepath)
            removed += 1

    print(f"Removed {removed} raw files from: {RAW_DIR}")


def run_case_012_inspection(keep_raw=False):
    print("Case:", PAPER_ID)
    print("GEO accession:", GEO_ACCESSION)
    print("Raw directory:", RAW_DIR)
    print("Keep raw files:", keep_raw)
    print("Mode: h5py-only inspection, no scanpy required")

    ensure_directories()
    download_raw_tar()
    extract_raw_tar()

    files = list_gsm_files()
    h5_files, other_files = classify_files(files)
    inspect_metadata_vs_files(files)
    inspect_all_h5_files(h5_files)

    print("\n=== Tool module implications ===")
    print("- src/tenx_h5_standardizer.py should support per-sample 10x h5 files.")
    print("- src/metadata_checker.py should validate GSM and patient IDs from filenames.")
    print("- src/gene_mapping_and_deduplication.py should run after concatenation.")
    print("- scanpy is useful for final processing, but h5py is enough for read-only inspection.")

    if not keep_raw:
        cleanup_raw_files()

    print("\nInspection completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Inspect Case 012: 2024_Brand / GSE233615"
    )
    parser.add_argument(
        "--keep-raw",
        action="store_true",
        help="Keep downloaded raw files after inspection."
    )

    args = parser.parse_args()
    run_case_012_inspection(keep_raw=args.keep_raw)
