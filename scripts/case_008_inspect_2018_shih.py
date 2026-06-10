import os
import tarfile
import argparse
from urllib.request import urlretrieve

import pandas as pd


PAPER_ID = "2018_Shih"
GEO_ACCESSION = "GSE118828"

BASE_DIR = os.path.join("data", "case_008_2018_shih")
RAW_DIR = os.path.join(BASE_DIR, "raw")
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")

RAW_TAR_FILENAME = "GSE118828_RAW.tar"
RAW_TAR_PATH = os.path.join(RAW_DIR, RAW_TAR_FILENAME)

GEO_DOWNLOAD_URL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE118828&format=file"


SAMPLE_METADATA = {
    "GSM3348303": {
        "patient_id": "553",
        "dataset_id": "GSE118828",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Peritoneal",
        "tumour_grade": "HGSC",
        "tumour_stage": "IV",
        "histological_subtype": "Serous",
        "patient_ethnicity": "White",
    },
    "GSM3348305": {
        "patient_id": "589",
        "dataset_id": "GSE118828",
        "cancer_type": "Ovarian",
        "tumor_site": "Metastasis",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "HGSC",
        "tumour_stage": "IIIC",
        "histological_subtype": "Serous",
        "patient_ethnicity": "White",
    },
    "GSM3348306": {
        "patient_id": "618",
        "dataset_id": "GSE118828",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "Yes",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "HGSC",
        "tumour_stage": "IIIC",
        "histological_subtype": "Serous",
        "patient_ethnicity": "African American",
    },
    "GSM3348308": {
        "patient_id": "568",
        "dataset_id": "GSE118828",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "LGSC",
        "tumour_stage": "IA",
        "histological_subtype": "Serous",
        "patient_ethnicity": "White",
    },
    "GSM3348309": {
        "patient_id": "580",
        "dataset_id": "GSE118828",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Fallopian Tube",
        "tumour_grade": "HGSC",
        "tumour_stage": "IIIB",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    },
    "GSM3348310": {
        "patient_id": "580",
        "dataset_id": "GSE118828",
        "cancer_type": "Ovarian",
        "tumor_site": "Metastasis",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "HGSC",
        "tumour_stage": "IIIB",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    },
    "GSM3348311": {
        "patient_id": "589",
        "dataset_id": "GSE118828",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "HGSC",
        "tumour_stage": "IIIC",
        "histological_subtype": "Serous",
        "patient_ethnicity": "White",
    },
    "GSM3348312": {
        "patient_id": "589",
        "dataset_id": "GSE118828",
        "cancer_type": "Ovarian",
        "tumor_site": "Metastasis",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "HGSC",
        "tumour_stage": "IIIC",
        "histological_subtype": "Serous",
        "patient_ethnicity": "White",
    },
    "GSM3348313": {
        "patient_id": "600",
        "dataset_id": "GSE118828",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "HGSC",
        "tumour_stage": "IIIA",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    },
    "GSM3348314": {
        "patient_id": "600",
        "dataset_id": "GSE118828",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "HGSC",
        "tumour_stage": "IIIA",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    },
    "GSM3348315": {
        "patient_id": "600",
        "dataset_id": "GSE118828",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "HGSC",
        "tumour_stage": "IIIA",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    },
    "GSM3348316": {
        "patient_id": "600",
        "dataset_id": "GSE118828",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "HGSC",
        "tumour_stage": "IIIA",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    },
    "GSM3348317": {
        "patient_id": "626",
        "dataset_id": "GSE118828",
        "cancer_type": "Ovarian",
        "tumor_site": "Metastasis",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "LGSC",
        "tumour_stage": "IIIC",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Asian",
    },
    "GSM3348318": {
        "patient_id": "626",
        "dataset_id": "GSE118828",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "LGSC",
        "tumour_stage": "IIIC",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Asian",
    },
    "GSM3348319": {
        "patient_id": "349",
        "dataset_id": "GSE118828",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "HGSC",
        "tumour_stage": "IIIA",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Unknown",
    },
    "GSM3348320": {
        "patient_id": "349",
        "dataset_id": "GSE118828",
        "cancer_type": "Ovarian",
        "tumor_site": "Metastasis",
        "metastasis_site": "Unknown",
        "tumor_treatment": "No",
        "cancer_site_origin": "Unknown",
        "tumour_grade": "HGSC",
        "tumour_stage": "IIIA",
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
    first_col_values = df_preview[first_col].astype(str)

    looks_like_cell_id_column = first_col.lower() in {"cellid", "cell_id", "cell", "barcode", "barcodes"}
    remaining_columns = list(df_preview.columns[1:])

    print("\n=== Interpretation ===")
    print("Detected format: compressed CSV count table")
    print("First column name:", first_col)
    print("First column looks like cell ID column:", looks_like_cell_id_column)
    print("Example first cell IDs:", list(first_col_values.head(5)))
    print("Remaining columns are treated as gene symbols.")
    print("Raw orientation: cells x genes")
    print("Needs transpose for AnnData: False")

    print("\nCorrect AnnData construction rule:")
    print("- Read csv.gz with pandas.")
    print("- Detect CellId column.")
    print("- Set obs_names from CellId.")
    print("- Drop CellId from expression matrix.")
    print("- Build AnnData with AnnData(X), no transpose.")
    print("- Assign var_names from remaining columns.")
    print("- Use gene symbols as dataset_gene_symbol for mapping.")


def cleanup_raw_files() -> None:
    print("\n=== Cleaning raw downloaded files ===")

    removed = 0

    for filename in os.listdir(RAW_DIR):
        filepath = os.path.join(RAW_DIR, filename)

        if os.path.isfile(filepath):
            os.remove(filepath)
            removed += 1

    print(f"Removed {removed} raw files from: {RAW_DIR}")


def run_case_008_inspection(keep_raw: bool = False) -> None:
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
        print(f"- {gsm}: patient {meta['patient_id']} / {meta['tumor_site']} / {meta['tumour_grade']} / {meta['tumour_stage']}")

    print("\n=== Tool module implications ===")
    print("- src/count_table_standardizer.py should support cell-by-gene CSV count tables.")
    print("- src/orientation_detector.py should detect CellId-column layouts and avoid transpose.")
    print("- src/metadata_checker.py should report skipped GSM files without curated metadata.")
    print("- src/gene_mapping_and_deduplication.py should use dataset_gene_symbol mapping for var_names.")

    if not keep_raw:
        cleanup_raw_files()

    print("\nInspection completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Inspect Case 008: 2018_Shih / GSE118828"
    )
    parser.add_argument(
        "--keep-raw",
        action="store_true",
        help="Keep downloaded raw files after inspection."
    )

    args = parser.parse_args()
    run_case_008_inspection(keep_raw=args.keep_raw)
