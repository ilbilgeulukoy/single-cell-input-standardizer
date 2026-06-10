import os
import tarfile
import gzip
import argparse
from urllib.request import urlretrieve

import pandas as pd
from scipy.io import mmread


PAPER_ID = "2022_Junfen"
GEO_ACCESSION = "GSE184880"

BASE_DIR = os.path.join("data", "case_006_2022_junfen")
RAW_DIR = os.path.join(BASE_DIR, "raw")
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")

RAW_TAR_FILENAME = "GSE184880_RAW.tar"
RAW_TAR_PATH = os.path.join(RAW_DIR, RAW_TAR_FILENAME)

GEO_DOWNLOAD_URL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE184880&format=file"


SAMPLE_METADATA = {
    "GSM5599225": {
        "patient_id": "Cancer1",
        "dataset_id": "GSE184880",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "HGSC",
        "tumour_stage": "IIIb",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    },
    "GSM5599226": {
        "patient_id": "Cancer2",
        "dataset_id": "GSE184880",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "HGSC",
        "tumour_stage": "IIb",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    },
    "GSM5599227": {
        "patient_id": "Cancer3",
        "dataset_id": "GSE184880",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "HGSC",
        "tumour_stage": "IC2",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    },
    "GSM5599228": {
        "patient_id": "Cancer4",
        "dataset_id": "GSE184880",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "HGSC",
        "tumour_stage": "IC2",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    },
    "GSM5599229": {
        "patient_id": "Cancer5",
        "dataset_id": "GSE184880",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "HGSC",
        "tumour_stage": "IIb",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    },
    "GSM5599230": {
        "patient_id": "Cancer6",
        "dataset_id": "GSE184880",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "HGSC",
        "tumour_stage": "IIIc",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    },
    "GSM5599231": {
        "patient_id": "Cancer7",
        "dataset_id": "GSE184880",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "HGSC",
        "tumour_stage": "IC2",
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

    cancer_files = [
        f for f in files
        if "_Cancer" in f and f.endswith((".mtx.gz", ".tsv.gz"))
    ]

    print("\n=== Detected all GSM files ===")
    print("Number of all GSM files:", len(files))

    print("\n=== Detected Cancer GSM files ===")
    print("Number of Cancer GSM files:", len(cancer_files))

    for f in cancer_files:
        print("-", f)

    return cancer_files


def group_files_by_gsm(files: list[str]) -> dict[str, dict[str, str]]:
    grouped = {}

    for f in files:
        gsm = f.split("_")[0]
        grouped.setdefault(gsm, {})

        if ".matrix." in f and f.endswith(".mtx.gz"):
            grouped[gsm]["matrix"] = f
        elif ".genes." in f and f.endswith(".tsv.gz"):
            grouped[gsm]["genes"] = f
        elif ".features." in f and f.endswith(".tsv.gz"):
            grouped[gsm]["features"] = f
        elif ".barcodes." in f and f.endswith(".tsv.gz"):
            grouped[gsm]["barcodes"] = f
        else:
            grouped[gsm].setdefault("other", [])
            grouped[gsm]["other"].append(f)

    print("\n=== Grouped Cancer files by GSM ===")

    for gsm, roles in sorted(grouped.items()):
        print(f"\n{gsm}")
        print("  matrix:", roles.get("matrix"))
        print("  genes:", roles.get("genes"))
        print("  features:", roles.get("features"))
        print("  barcodes:", roles.get("barcodes"))
        if roles.get("other"):
            print("  other:", roles.get("other"))

    return grouped


def inspect_first_triplet(grouped: dict[str, dict[str, str]]) -> None:
    complete_gsms = [
        gsm for gsm, roles in grouped.items()
        if roles.get("matrix") and roles.get("barcodes") and (roles.get("genes") or roles.get("features"))
    ]

    if not complete_gsms:
        raise FileNotFoundError("No complete Cancer matrix/genes-or-features/barcodes triplet found.")

    gsm = sorted(complete_gsms)[0]
    roles = grouped[gsm]

    matrix_path = os.path.join(RAW_DIR, roles["matrix"])
    gene_table_name = roles.get("genes") or roles.get("features")
    gene_table_role = "genes" if roles.get("genes") else "features"
    gene_table_path = os.path.join(RAW_DIR, gene_table_name)
    barcodes_path = os.path.join(RAW_DIR, roles["barcodes"])

    print("\n=== Inspecting first Cancer 10x-like triplet ===")
    print("GSM:", gsm)
    print("Matrix:", roles["matrix"])
    print("Gene table role:", gene_table_role)
    print("Gene table:", gene_table_name)
    print("Barcodes:", roles["barcodes"])

    with gzip.open(matrix_path, "rt") as f:
        matrix = mmread(f).tocsr()

    gene_table = pd.read_csv(gene_table_path, sep="\t", header=None, compression="gzip")
    barcodes = pd.read_csv(barcodes_path, header=None, compression="gzip")

    print("\nMatrix shape from Matrix Market:")
    print(matrix.shape)

    print("\nGene table shape:")
    print(gene_table.shape)

    print("\nBarcodes shape:")
    print(barcodes.shape)

    print("\nFirst gene table rows:")
    print(gene_table.head())

    print("\nFirst barcodes rows:")
    print(barcodes.head())

    if matrix.shape[0] == gene_table.shape[0] and matrix.shape[1] == barcodes.shape[0]:
        orientation = "genes_x_cells_before_transpose"
        needs_transpose = True
    elif matrix.shape[0] == barcodes.shape[0] and matrix.shape[1] == gene_table.shape[0]:
        orientation = "cells_x_genes_already"
        needs_transpose = False
    else:
        orientation = "ambiguous"
        needs_transpose = None

    gene_symbol_col = 1 if gene_table.shape[1] > 1 else 0
    duplicated_gene_symbols = gene_table.iloc[:, gene_symbol_col].astype(str).duplicated().sum()

    print("\n=== Interpretation ===")
    print("Detected format: Cancer-only 10x-like Matrix Market triplet")
    print("Gene table variant:", gene_table_role)
    print("Filename role separator pattern: .matrix. / .genes. / .barcodes.")
    print("Matrix orientation:", orientation)
    print("Needs transpose for AnnData:", needs_transpose)
    print("Gene symbol column:", gene_symbol_col)
    print("Duplicated gene symbols in first sample:", duplicated_gene_symbols)

    print("\nCorrect AnnData construction rule:")
    print("- Filter only Cancer samples when both cancer and non-malignant files exist.")
    print("- Group files by GSM accession.")
    print("- Accept genes.tsv.gz or features.tsv.gz as gene table.")
    print("- Support dot-separated role patterns such as .matrix. / .genes. / .barcodes.")
    print("- Read barcodes explicitly and prefix them with GSM.")
    print("- Transpose matrix if genes x cells.")
    print("- Apply gene mapping and duplicate-gene aggregation after concatenation if needed.")


def cleanup_raw_files() -> None:
    print("\n=== Cleaning raw downloaded files ===")

    removed = 0

    for filename in os.listdir(RAW_DIR):
        filepath = os.path.join(RAW_DIR, filename)

        if os.path.isfile(filepath):
            os.remove(filepath)
            removed += 1

    print(f"Removed {removed} raw files from: {RAW_DIR}")


def run_case_006_inspection(keep_raw: bool = False) -> None:
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

    print("\n=== Metadata Cancer samples ===")
    print("Curated sample metadata entries:", len(SAMPLE_METADATA))
    for gsm, meta in SAMPLE_METADATA.items():
        print(f"- {gsm}: {meta['patient_id']} / {meta['tumour_grade']} / {meta['tumour_stage']} / {meta['histological_subtype']}")

    print("\n=== Tool module implications ===")
    print("- src/tenx_mtx_standardizer.py should support dot-separated role patterns.")
    print("- src/tenx_mtx_standardizer.py should support Cancer-only sample filtering.")
    print("- src/metadata_checker.py should detect whether filtered GSM files have metadata entries.")
    print("- src/gene_mapping_and_deduplication.py may be applied after concatenation.")

    if not keep_raw:
        cleanup_raw_files()

    print("\nInspection completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Inspect Case 006: 2022_Junfen / GSE184880"
    )
    parser.add_argument(
        "--keep-raw",
        action="store_true",
        help="Keep downloaded raw files after inspection."
    )

    args = parser.parse_args()
    run_case_006_inspection(keep_raw=args.keep_raw)
