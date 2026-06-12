import os
import tarfile
import gzip
import argparse
from urllib.request import urlretrieve

import pandas as pd


PAPER_ID = "2024_Denisenko"
GEO_ACCESSION = "GSE211956"

BASE_DIR = os.path.join("data", "case_015_2024_denisenko")
RAW_DIR = os.path.join(BASE_DIR, "raw")

RAW_TAR_FILENAME = "GSE211956_RAW.tar"
RAW_TAR_PATH = os.path.join(RAW_DIR, RAW_TAR_FILENAME)

GEO_DOWNLOAD_URL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE211956&format=file"


SAMPLE_METADATA = {
    "GSM6506105": {
        "patient_id": "Y2",
        "dataset_id": "GSE211956",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "Yes",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "HGSC",
        "tumour_stage": "III-IV",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    },
    "GSM6506106": {
        "patient_id": "Y3",
        "dataset_id": "GSE211956",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "Yes",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "HGSC",
        "tumour_stage": "III-IV",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    },
    "GSM6506107": {
        "patient_id": "Y5",
        "dataset_id": "GSE211956",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "Yes",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "HGSC",
        "tumour_stage": "III-IV",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    },
    "GSM6506108": {
        "patient_id": "MJ10",
        "dataset_id": "GSE211956",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "Yes",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "HGSC",
        "tumour_stage": "III-IV",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    },
    "GSM6506109": {
        "patient_id": "MJ11",
        "dataset_id": "GSE211956",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "Yes",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "HGSC",
        "tumour_stage": "III-IV",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
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
    txt_gz_files = [
        f for f in files
        if f.endswith(".txt.gz")
    ]

    txt_files = [
        f for f in files
        if f.endswith(".txt") and not f.endswith(".txt.gz")
    ]

    other_files = [
        f for f in files
        if f not in txt_gz_files + txt_files
    ]

    print("\n=== Classified files ===")
    print("Compressed txt.gz count tables:", len(txt_gz_files))
    for f in txt_gz_files:
        print("-", f)

    print("\nUncompressed txt files:", len(txt_files))
    for f in txt_files:
        print("-", f)

    print("\nOther files:", len(other_files))
    for f in other_files:
        print("-", f)

    return {
        "txt_gz": txt_gz_files,
        "txt": txt_files,
        "other": other_files,
    }


def inspect_metadata_vs_files(files):
    detected_gsms = sorted({
        f.split("_")[0].strip()
        for f in files
        if f.startswith("GSM")
    })

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


def count_gzip_lines(filepath):
    line_count = 0
    with gzip.open(filepath, "rt") as handle:
        for _ in handle:
            line_count += 1
    return line_count


def inspect_txt_table(filename, nrows_preview=5):
    filepath = os.path.join(RAW_DIR, filename)
    gsm = filename.split("_")[0].strip()

    print("\n=== Inspecting TXT count table ===")
    print("File:", filename)
    print("GSM:", gsm)

    print("File size MB:", round(os.path.getsize(filepath) / 1024**2, 2))

    preview = pd.read_csv(filepath, sep="\t", compression="gzip", nrows=nrows_preview)

    print("\nPreview shape:")
    print(preview.shape)

    print("\nPreview columns first 10:")
    print(preview.columns[:10].tolist())

    print("\nPreview first rows:")
    print(preview.head())

    first_col = preview.columns[0]
    first_values = preview.iloc[:, 0].astype(str).head(10).tolist()

    print("\nFirst column name:", first_col)
    print("First column first values:", first_values)

    n_columns = len(preview.columns)
    n_cell_like_columns = n_columns - 1

    print("\nColumn count from header:", n_columns)
    print("Cell-like columns if first column is gene:", n_cell_like_columns)

    line_count = count_gzip_lines(filepath)
    n_data_rows = max(line_count - 1, 0)

    print("Line count including header:", line_count)
    print("Data rows:", n_data_rows)

    first_col_lower = str(first_col).lower()

    gene_like_values = any(
        value.startswith(("ENSG", "RP", "AL", "MT-", "MIR", "OR"))
        or value.isupper()
        for value in first_values
    )

    if first_col_lower in {"gene", "genes", "symbol", "feature", "features", "unnamed: 0"} or gene_like_values:
        orientation = "genes_x_cells_before_transpose"
        needs_transpose = True
        expected_anndata_shape = (n_cell_like_columns, n_data_rows)
    else:
        orientation = "ambiguous"
        needs_transpose = None
        expected_anndata_shape = None

    duplicate_gene_preview = preview.iloc[:, 0].astype(str).duplicated().sum()

    print("\n=== Interpretation for this file ===")
    print("Detected table type: compressed TXT count table")
    print("Likely orientation:", orientation)
    print("Needs transpose for AnnData:", needs_transpose)
    print("Expected AnnData shape from header/line count:", expected_anndata_shape)
    print("Duplicate first-column values in preview:", duplicate_gene_preview)

    patient_token = None
    parts = filename.split("_")
    if len(parts) >= 3:
        patient_token = parts[2].replace(".txt.gz", "").strip()

    print("\nFilename-derived patient token:", patient_token)

    if gsm in SAMPLE_METADATA:
        print("Curated patient_id:", SAMPLE_METADATA[gsm]["patient_id"])
        print("Patient token matches curated:", patient_token == SAMPLE_METADATA[gsm]["patient_id"])
    else:
        print("GSM not present in curated metadata.")

    return {
        "gsm": gsm,
        "filename": filename,
        "n_data_rows": n_data_rows,
        "n_columns": n_columns,
        "n_cell_like_columns": n_cell_like_columns,
        "orientation": orientation,
        "needs_transpose": needs_transpose,
        "expected_anndata_shape": expected_anndata_shape,
        "patient_token": patient_token,
    }


def inspect_all_txt_tables(txt_gz_files):
    print("\n=== All TXT count table summaries ===")

    summaries = []

    for filename in txt_gz_files:
        summary = inspect_txt_table(filename)
        summaries.append(summary)

    print("\n=== Compact summary table ===")

    total_cells = 0

    for s in summaries:
        n_cells = None
        n_genes = None

        if s["expected_anndata_shape"]:
            n_cells = s["expected_anndata_shape"][0]
            n_genes = s["expected_anndata_shape"][1]
            total_cells += n_cells

        sample_info = SAMPLE_METADATA.get(s["gsm"], {})

        print(
            s["gsm"],
            "filename=", s["filename"],
            "patient_id=", sample_info.get("patient_id"),
            "cells=", n_cells,
            "genes=", n_genes,
            "orientation=", s["orientation"],
            "patient_token=", s["patient_token"],
        )

    print("\nTotal cells estimated from headers:", total_cells)
    print("Unique gene counts estimated:", sorted(set(
        s["expected_anndata_shape"][1]
        for s in summaries
        if s["expected_anndata_shape"]
    )))
    print("Number of txt.gz files inspected:", len(summaries))


def inspect_script_risks():
    print("\n=== Original script risk checks ===")
    print("- gzip -dk creates uncompressed .txt files but the script still reads .txt.gz files.")
    print("- Decompression step is probably unnecessary and increases disk usage.")
    print("- sc.read_text(fp).T implies raw table is genes x cells.")
    print("- The NaN check after common-gene filtering only checks the last loop variable adata, not every sample.")
    print("- After concat, dense np.isnan(adata_combined.X) may be memory-heavy if X is large.")
    print("- Gene mapping files are written into directory_input instead of directory_annotated in the script.")
    print("- Final last-check section inspects `adata`, but after writing new_adata the script reloads into adata only at the end.")


def cleanup_raw_files():
    print("\n=== Cleaning raw downloaded files ===")

    removed = 0

    for filename in os.listdir(RAW_DIR):
        filepath = os.path.join(RAW_DIR, filename)

        if os.path.isfile(filepath):
            os.remove(filepath)
            removed += 1

    print(f"Removed {removed} raw files from: {RAW_DIR}")


def run_case_015_inspection(keep_raw=False):
    print("Case:", PAPER_ID)
    print("GEO accession:", GEO_ACCESSION)
    print("Raw directory:", RAW_DIR)
    print("Keep raw files:", keep_raw)
    print("Mode: scanpy-free compressed TXT count table inspection")

    ensure_directories()
    download_raw_tar()
    extract_raw_tar()

    files = list_gsm_files()
    classified = classify_files(files)
    inspect_metadata_vs_files(files)
    inspect_all_txt_tables(classified["txt_gz"])
    inspect_script_risks()

    print("\n=== Tool module implications ===")
    print("- src/count_table_standardizer.py should support compressed TXT count tables.")
    print("- src/data_size_policy.py should avoid unnecessary gzip decompression.")
    print("- src/metadata_checker.py should validate GSM and patient tokens from filenames.")
    print("- src/gene_mapping_and_deduplication.py should run after concatenation.")
    print("- src/validation.py should use sparse-safe NaN checks.")

    if not keep_raw:
        cleanup_raw_files()

    print("\nInspection completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Inspect Case 015: 2024_Denisenko / GSE211956"
    )
    parser.add_argument(
        "--keep-raw",
        action="store_true",
        help="Keep downloaded raw files after inspection."
    )

    args = parser.parse_args()
    run_case_015_inspection(keep_raw=args.keep_raw)
