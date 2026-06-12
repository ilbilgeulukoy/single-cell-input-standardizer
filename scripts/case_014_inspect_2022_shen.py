import os
import tarfile
import gzip
import argparse
from urllib.request import urlretrieve

import pandas as pd
from scipy.io import mmread


PAPER_ID = "2022_Shen"
GEO_ACCESSION = "GSE191301"

BASE_DIR = os.path.join("data", "case_014_2022_shen")
RAW_DIR = os.path.join(BASE_DIR, "raw")

RAW_TAR_FILENAME = "GSE191301_RAW.tar"
RAW_TAR_PATH = os.path.join(RAW_DIR, RAW_TAR_FILENAME)

GEO_DOWNLOAD_URL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE191301&format=file"


SAMPLE_METADATA = {
    "GSM5743307": {
        "patient_id": "Pre-NACT1A",
        "tumor_site": "Metastasis",
        "metastasis_site": "Peritoneum",
        "tumor_treatment": "No",
        "cancer_site_origin": "Fallopian tube",
        "tumour_grade": "HGSC",
        "tumour_stage": "IIIc",
        "histological_subtype": "Serous",
    },
    "GSM5743308": {
        "patient_id": "Pre-NACT1B",
        "tumor_site": "Metastasis",
        "metastasis_site": "Omentum",
        "tumor_treatment": "No",
        "cancer_site_origin": "Fallopian Tube",
        "tumour_grade": "HGSC",
        "tumour_stage": "IIIc",
        "histological_subtype": "Serous",
    },
    "GSM5743309": {
        "patient_id": "Pre-NACT1C",
        "tumor_site": "Metastasis",
        "metastasis_site": "Ascites",
        "tumor_treatment": "No",
        "cancer_site_origin": "Fallopian Tube",
        "tumour_grade": "HGSC",
        "tumour_stage": "IIIc",
        "histological_subtype": "Serous",
    },
    "GSM5743310": {
        "patient_id": "Post-NACT1D",
        "tumor_site": "Metastasis",
        "metastasis_site": "Peritoneum",
        "tumor_treatment": "Yes",
        "cancer_site_origin": "Fallopian Tube",
        "tumour_grade": "HGSC",
        "tumour_stage": "IIIc",
        "histological_subtype": "Serous",
    },
    "GSM5743311": {
        "patient_id": "Post-NACT1E",
        "tumor_site": "Metastasis",
        "metastasis_site": "Omentum",
        "tumor_treatment": "Yes",
        "cancer_site_origin": "Fallopian Tube",
        "tumour_grade": "HGSC",
        "tumour_stage": "IIIc",
        "histological_subtype": "Serous",
    },
    "GSM5743312": {
        "patient_id": "Post-NACT1F",
        "tumor_site": "Metastasis",
        "metastasis_site": "Ascites",
        "tumor_treatment": "Yes",
        "cancer_site_origin": "Fallopian Tube",
        "tumour_grade": "HGSC",
        "tumour_stage": "IIIc",
        "histological_subtype": "Serous",
    },
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
    matrix_files = [
        f for f in files
        if "_matrix" in f and f.endswith(".mtx.gz")
    ]
    feature_files = [
        f for f in files
        if ("_features" in f or "_genes" in f) and f.endswith(".tsv.gz")
    ]
    barcode_files = [
        f for f in files
        if "_barcodes" in f and f.endswith(".tsv.gz")
    ]
    other_files = [
        f for f in files
        if f not in matrix_files + feature_files + barcode_files
    ]

    print("\n=== Classified files ===")

    print("Matrix files:", len(matrix_files))
    for f in matrix_files:
        print("-", f)

    print("\nFeature/gene files:", len(feature_files))
    for f in feature_files:
        print("-", f)

    print("\nBarcode files:", len(barcode_files))
    for f in barcode_files:
        print("-", f)

    print("\nOther files:", len(other_files))
    for f in other_files:
        print("-", f)

    return {
        "matrix": matrix_files,
        "features": feature_files,
        "barcodes": barcode_files,
        "other": other_files,
    }


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


def build_triplets(classified):
    triplets = {}

    all_gsms = sorted(
        {
            f.split("_")[0].strip()
            for role in ["matrix", "features", "barcodes"]
            for f in classified[role]
        }
    )

    for gsm in all_gsms:
        matrix = [f for f in classified["matrix"] if f.startswith(gsm + "_")]
        features = [f for f in classified["features"] if f.startswith(gsm + "_")]
        barcodes = [f for f in classified["barcodes"] if f.startswith(gsm + "_")]

        triplets[gsm] = {
            "matrix": matrix,
            "features": features,
            "barcodes": barcodes,
        }

    print("\n=== Per-GSM triplet matching ===")

    for gsm, roles in triplets.items():
        print("\nGSM:", gsm)
        print("  matrix:", roles["matrix"])
        print("  features:", roles["features"])
        print("  barcodes:", roles["barcodes"])
        print(
            "  complete:",
            len(roles["matrix"]) == 1 and len(roles["features"]) == 1 and len(roles["barcodes"]) == 1,
        )

    incomplete = [
        gsm for gsm, roles in triplets.items()
        if not (len(roles["matrix"]) == 1 and len(roles["features"]) == 1 and len(roles["barcodes"]) == 1)
    ]

    print("\nIncomplete triplets:", incomplete)

    return triplets


def inspect_one_triplet(gsm, roles):
    if not (len(roles["matrix"]) == 1 and len(roles["features"]) == 1 and len(roles["barcodes"]) == 1):
        print(f"\nSkipping {gsm}, incomplete triplet.")
        return None

    matrix_file = roles["matrix"][0]
    features_file = roles["features"][0]
    barcodes_file = roles["barcodes"][0]

    matrix_path = os.path.join(RAW_DIR, matrix_file)
    features_path = os.path.join(RAW_DIR, features_file)
    barcodes_path = os.path.join(RAW_DIR, barcodes_file)

    print("\n=== Inspecting triplet ===")
    print("GSM:", gsm)
    print("Matrix:", matrix_file)
    print("Features:", features_file)
    print("Barcodes:", barcodes_file)

    with gzip.open(matrix_path, "rt") as f:
        matrix = mmread(f).tocsr()

    features = pd.read_csv(features_path, sep="\t", header=None, compression="gzip")
    barcodes = pd.read_csv(barcodes_path, sep="\t", header=None, compression="gzip")

    print("Matrix shape:", matrix.shape)
    print("Features shape:", features.shape)
    print("Barcodes shape:", barcodes.shape)

    if matrix.shape[0] == features.shape[0] and matrix.shape[1] == barcodes.shape[0]:
        orientation = "genes_x_cells_before_transpose"
        needs_transpose = True
        expected_anndata_shape = (matrix.shape[1], matrix.shape[0])
    elif matrix.shape[0] == barcodes.shape[0] and matrix.shape[1] == features.shape[0]:
        orientation = "cells_x_genes_already"
        needs_transpose = False
        expected_anndata_shape = matrix.shape
    else:
        orientation = "ambiguous"
        needs_transpose = None
        expected_anndata_shape = None

    gene_symbol_col = 1 if features.shape[1] > 1 else 0
    duplicated_gene_symbols = features.iloc[:, gene_symbol_col].astype(str).duplicated().sum()

    first_features = []
    for _, row in features.head(5).iterrows():
        first_features.append(" / ".join(str(x) for x in row.tolist()))

    first_barcodes = barcodes.iloc[:5, 0].astype(str).tolist()

    print("Matrix orientation:", orientation)
    print("Needs transpose:", needs_transpose)
    print("Expected AnnData shape:", expected_anndata_shape)
    print("Gene symbol column:", gene_symbol_col)
    print("Duplicated gene symbols:", duplicated_gene_symbols)
    print("First feature rows:")
    for x in first_features:
        print("-", x)
    print("First barcodes:", first_barcodes)

    return {
        "gsm": gsm,
        "matrix_file": matrix_file,
        "features_file": features_file,
        "barcodes_file": barcodes_file,
        "matrix_shape": matrix.shape,
        "features_shape": features.shape,
        "barcodes_shape": barcodes.shape,
        "orientation": orientation,
        "needs_transpose": needs_transpose,
        "expected_anndata_shape": expected_anndata_shape,
        "gene_symbol_col": gene_symbol_col,
        "duplicated_gene_symbols": duplicated_gene_symbols,
        "first_features": first_features,
        "first_barcodes": first_barcodes,
    }


def inspect_all_triplets(triplets):
    print("\n=== All triplet shape summaries ===")

    summaries = []

    for gsm in sorted(triplets):
        summary = inspect_one_triplet(gsm, triplets[gsm])
        if summary:
            summaries.append(summary)

    print("\n=== Compact summary table ===")
    total_cells = 0

    for s in summaries:
        n_cells = s["expected_anndata_shape"][0] if s["expected_anndata_shape"] else None
        n_genes = s["expected_anndata_shape"][1] if s["expected_anndata_shape"] else None
        total_cells += n_cells or 0

        sample_info = SAMPLE_METADATA.get(s["gsm"], {})
        print(
            s["gsm"],
            "patient_id=", sample_info.get("patient_id"),
            "treatment=", sample_info.get("tumor_treatment"),
            "site=", sample_info.get("metastasis_site"),
            "cells=", n_cells,
            "genes=", n_genes,
            "dup_gene_symbols=", s["duplicated_gene_symbols"],
        )

    print("\nTotal cells across readable triplets:", total_cells)
    print("Unique gene counts:", sorted(set(s["expected_anndata_shape"][1] for s in summaries if s["expected_anndata_shape"])))
    print("Number of complete triplets inspected:", len(summaries))

    return summaries


def inspect_script_risks():
    print("\n=== Original script risk checks ===")
    print("- Initial adata.obs_names line uses variable `gsm` before it is defined.")
    print("- Loop recomputes file lists but then uses mtx_files[0], feat_files[0], bc_files[0].")
    print("- This risks assigning the first triplet to every GSM.")
    print("- Correct rule: select files with filename starting with the current GSM prefix.")
    print("- File pattern may use '_matrix.' rather than '_matrix_', so detector should support both dot and underscore variants.")


def cleanup_raw_files():
    print("\n=== Cleaning raw downloaded files ===")

    removed = 0

    for filename in os.listdir(RAW_DIR):
        filepath = os.path.join(RAW_DIR, filename)

        if os.path.isfile(filepath):
            os.remove(filepath)
            removed += 1

    print(f"Removed {removed} raw files from: {RAW_DIR}")


def run_case_014_inspection(keep_raw=False):
    print("Case:", PAPER_ID)
    print("GEO accession:", GEO_ACCESSION)
    print("Raw directory:", RAW_DIR)
    print("Keep raw files:", keep_raw)
    print("Mode: scanpy-free per-sample 10x Matrix Market inspection")

    ensure_directories()
    download_raw_tar()
    extract_raw_tar()

    files = list_gsm_files()
    classified = classify_files(files)
    inspect_metadata_vs_files(files)
    triplets = build_triplets(classified)
    inspect_all_triplets(triplets)
    inspect_script_risks()

    print("\n=== Tool module implications ===")
    print("- src/tenx_mtx_standardizer.py should match triplets by GSM prefix.")
    print("- src/metadata_checker.py should validate pre/post NACT sample metadata.")
    print("- src/file_pattern_detector.py should support both '_matrix.' and '_matrix_' naming.")
    print("- src/gene_mapping_and_deduplication.py should run after concatenation.")

    if not keep_raw:
        cleanup_raw_files()

    print("\nInspection completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Inspect Case 014: 2022_Shen / GSE191301"
    )
    parser.add_argument(
        "--keep-raw",
        action="store_true",
        help="Keep downloaded raw files after inspection."
    )

    args = parser.parse_args()
    run_case_014_inspection(keep_raw=args.keep_raw)
