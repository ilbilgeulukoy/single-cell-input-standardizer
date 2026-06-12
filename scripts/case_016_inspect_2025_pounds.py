import os
import gzip
import argparse
from urllib.request import urlretrieve
from collections import Counter

import pandas as pd


PAPER_ID = "2025_Pounds"
GEO_ACCESSION = "GSE281120"

BASE_DIR = os.path.join("data", "case_016_2025_pounds")
RAW_DIR = os.path.join(BASE_DIR, "raw")

COUNT_FILENAME = "GSE281120_counts.csv.gz"
COUNT_PATH = os.path.join(RAW_DIR, COUNT_FILENAME)

COUNT_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE281nnn/GSE281120/suppl/GSE281120%5Fcounts%2Ecsv%2Egz"


SAMPLE_METADATA = {
    "GSM8611246": {
        "patient_id": "Pt1",
        "sample_suffix": "1",
        "tumor_treatment": "Yes",
    },
    "GSM8611247": {
        "patient_id": "Pt2",
        "sample_suffix": "2",
        "tumor_treatment": "No",
    },
    "GSM8611248": {
        "patient_id": "Pt3",
        "sample_suffix": "3",
        "tumor_treatment": "Yes",
    },
    "GSM8611251": {
        "patient_id": "Pt6",
        "sample_suffix": "6",
        "tumor_treatment": "Yes",
    },
    "GSM8611252": {
        "patient_id": "Pt7",
        "sample_suffix": "7",
        "tumor_treatment": "No",
    },
    "GSM8611254": {
        "patient_id": "Pt8",
        "sample_suffix": "8",
        "tumor_treatment": "Yes",
    },
    "GSM8611255": {
        "patient_id": "Pt9",
        "sample_suffix": "9",
        "tumor_treatment": "No",
    },
    "GSM8611257": {
        "patient_id": "Pt10",
        "sample_suffix": "10",
        "tumor_treatment": "No",
    },
}


SUFFIX_TO_GSM = {
    "1": "GSM8611246",
    "2": "GSM8611247",
    "3": "GSM8611248",
    "6": "GSM8611251",
    "7": "GSM8611252",
    "8": "GSM8611254",
    "9": "GSM8611255",
    "10": "GSM8611257",
}


EXCLUDED_SUFFIXES = {
    "5": {
        "reason": "normal / non-metastasis sample excluded by original script",
        "geo_sample_hint": "s271",
        "geo_accession_hint": "GSM8611249",
    }
}


def ensure_directories():
    os.makedirs(RAW_DIR, exist_ok=True)


def download_count_matrix():
    if os.path.exists(COUNT_PATH):
        print(f"Count matrix already exists: {COUNT_PATH}")
        return

    print("Downloading GEO supplementary count matrix...")
    print("URL:", COUNT_URL)
    print("Output:", COUNT_PATH)

    urlretrieve(COUNT_URL, COUNT_PATH)

    print("Download completed.")


def read_csv_header(filepath):
    with gzip.open(filepath, "rt") as handle:
        header = handle.readline().rstrip("\n").split(",")

    return header


def count_gzip_lines(filepath):
    line_count = 0
    with gzip.open(filepath, "rt") as handle:
        for _ in handle:
            line_count += 1
    return line_count


def get_suffix_from_cell(cell_name):
    parts = str(cell_name).split("_")
    if len(parts) < 2:
        return "UNKNOWN"
    return parts[-1]


def inspect_count_matrix():
    print("\n=== File check ===")
    print("File exists:", os.path.exists(COUNT_PATH))

    if not os.path.exists(COUNT_PATH):
        raise FileNotFoundError(COUNT_PATH)

    print("File size MB:", round(os.path.getsize(COUNT_PATH) / 1024**2, 2))

    print("\n=== Header inspection ===")
    header = read_csv_header(COUNT_PATH)

    first_col = header[0]
    cell_columns = header[1:]

    print("First column name:", first_col)
    print("Number of columns including gene column:", len(header))
    print("Number of cell columns:", len(cell_columns))
    print("First 10 cell columns:")
    for col in cell_columns[:10]:
        print("-", col)

    suffix_counts = Counter(get_suffix_from_cell(cell) for cell in cell_columns)

    print("\n=== Cell suffix distribution from header ===")
    for suffix, count in sorted(suffix_counts.items(), key=lambda x: (len(x[0]), x[0])):
        print(suffix, count)

    mapped_suffixes = set(SUFFIX_TO_GSM)
    detected_suffixes = set(suffix_counts)

    print("\nMapped tumor suffixes:")
    print(sorted(mapped_suffixes, key=lambda x: (len(x), x)))

    print("\nDetected suffixes missing from sample_to_gsm:")
    print(sorted(detected_suffixes - mapped_suffixes, key=lambda x: (len(x), x)))

    print("\nMapped suffixes missing from count matrix:")
    print(sorted(mapped_suffixes - detected_suffixes, key=lambda x: (len(x), x)))

    print("\nExcluded suffix policy:")
    for suffix, info in EXCLUDED_SUFFIXES.items():
        print(f"- suffix {suffix}: {info}")

    selected_cell_count = sum(suffix_counts.get(suffix, 0) for suffix in SUFFIX_TO_GSM)
    excluded_cell_count = sum(suffix_counts.get(suffix, 0) for suffix in EXCLUDED_SUFFIXES)

    print("\nSelected cells from mapped tumor suffixes:", selected_cell_count)
    print("Excluded cells from excluded suffixes:", excluded_cell_count)

    print("\n=== Preview matrix ===")
    preview = pd.read_csv(COUNT_PATH, compression="gzip", nrows=5)

    print("Preview shape:", preview.shape)
    print("Preview first 10 columns:", preview.columns[:10].tolist())
    print(preview.head())

    genes_preview = preview.iloc[:, 0].astype(str).tolist()

    print("\nFirst five gene values:")
    for gene in genes_preview:
        print("-", gene)

    print("\n=== Gene row count ===")
    line_count = count_gzip_lines(COUNT_PATH)
    n_genes = max(line_count - 1, 0)

    print("Line count including header:", line_count)
    print("Gene rows:", n_genes)

    print("\n=== Interpretation ===")
    print("Detected format: global compressed CSV count matrix")
    print("Raw orientation: genes x cells")
    print("Transpose required for AnnData: True")
    print("Sample identity source: cell barcode suffix")
    print("Normal/excluded sample rule: remove suffix 5")
    print("Expected AnnData shape before exclusion:", (len(cell_columns), n_genes))
    print("Expected AnnData shape after exclusion:", (selected_cell_count, n_genes))

    print("\n=== Per-sample selected cell counts ===")
    for suffix, gsm in sorted(SUFFIX_TO_GSM.items(), key=lambda x: (len(x[0]), x[0])):
        meta = SAMPLE_METADATA[gsm]
        print(
            gsm,
            "patient_id=", meta["patient_id"],
            "suffix=", suffix,
            "cells=", suffix_counts.get(suffix, 0),
            "tumor_treatment=", meta["tumor_treatment"],
        )


def inspect_script_risks():
    print("\n=== Original script risk checks ===")
    print("- The script downloads GSE281120_counts.csv.gz but later searches for a .csv file.")
    print("- list_abspath_sample is used for gzip decompression but is not defined in this script.")
    print("- Decompression is unnecessary; the .csv.gz file can be read directly.")
    print("- This is one global matrix, not per-GSM files.")
    print("- Sample identity must be parsed from cell barcode suffix.")
    print("- Suffix 5 represents the normal/non-metastasis sample and should be excluded for the tumor cohort.")
    print("- sample_to_gsm has missing suffixes 4 and 5 by design; suffix 5 is explicitly excluded.")
    print("- Dense AnnData construction from a 734 MB CSV may be memory-heavy.")
    print("- The tool should use data size policy and optionally chunked reading for large global CSV files.")


def cleanup_raw_files():
    print("\n=== Cleaning raw downloaded files ===")

    removed = 0

    for filename in os.listdir(RAW_DIR):
        filepath = os.path.join(RAW_DIR, filename)

        if os.path.isfile(filepath):
            os.remove(filepath)
            removed += 1

    print(f"Removed {removed} raw files from: {RAW_DIR}")


def run_case_016_inspection(keep_raw=False):
    print("Case:", PAPER_ID)
    print("GEO accession:", GEO_ACCESSION)
    print("Raw directory:", RAW_DIR)
    print("Keep raw files:", keep_raw)
    print("Mode: scanpy-free global CSV count matrix inspection")

    ensure_directories()
    download_count_matrix()
    inspect_count_matrix()
    inspect_script_risks()

    print("\n=== Tool module implications ===")
    print("- src/global_matrix_metadata_standardizer.py should support global CSV count matrices.")
    print("- src/sample_id_parser.py should parse sample suffixes from cell barcodes.")
    print("- src/modality_filter.py or cohort_filter.py should exclude normal suffix 5.")
    print("- src/data_size_policy.py should handle large compressed CSV files.")
    print("- src/count_table_standardizer.py should support gene-by-cell CSV with first gene column.")
    print("- src/gene_mapping_and_deduplication.py should run after sample filtering.")

    if not keep_raw:
        cleanup_raw_files()

    print("\nInspection completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Inspect Case 016: 2025_Pounds / GSE281120"
    )
    parser.add_argument(
        "--keep-raw",
        action="store_true",
        help="Keep downloaded raw files after inspection."
    )

    args = parser.parse_args()
    run_case_016_inspection(keep_raw=args.keep_raw)
