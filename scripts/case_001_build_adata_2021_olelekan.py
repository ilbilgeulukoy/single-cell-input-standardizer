import os
import argparse
import tarfile
from urllib.request import urlretrieve

import pandas as pd
import numpy as np
import anndata as ad
from anndata import AnnData
from scipy.sparse import issparse


# ============================================================
# Case 001: 2021_Olelekan / GSE147082
#
# Goal:
# Build a standardized AnnData/h5ad object from GEO-downloaded
# compressed CSV count tables.
#
# Raw input pattern:
#   GSE147082_RAW.tar
#       ├── GSM4416534_PT-3232.csv.gz
#       ├── GSM4416535_PT-5150.csv.gz
#       └── ...
#
# Key data engineering issue:
#   The raw CSV tables are genes x cells:
#       rows    = genes
#       columns = cells
#
#   AnnData expects cells x genes:
#       rows    = observations / cells
#       columns = variables / genes
#
#   Therefore, we create:
#       AnnData(df.T)
# ============================================================


# -----------------------------
# 1. Project variables
# -----------------------------
PAPER_ID = "2021_Olelekan"
GEO_ACCESSION = "GSE147082"

BASE_DIR = os.path.join("data", "case_001_2021_olelekan")
RAW_DIR = os.path.join(BASE_DIR, "raw")
ANNOTATED_DIR = os.path.join(BASE_DIR, "annotated")
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")

RAW_TAR_FILENAME = "GSE147082_RAW.tar"
RAW_TAR_PATH = os.path.join(RAW_DIR, RAW_TAR_FILENAME)

GEO_DOWNLOAD_URL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE147082&format=file"

OUTPUT_H5AD_PATH = os.path.join(
    OUTPUT_DIR,
    "2021_Olelekan_GSE147082_standardized.h5ad"
)


# -----------------------------
# 2. Metadata column names
# -----------------------------
COL_SAMPLE_ID = "sample_id"
COL_PATIENT_ID = "patient_id"
COL_DATASET_ID = "dataset_id"
COL_CANCER_TYPE = "cancer_type"
COL_TUMOR_SITE = "tumor_site"
COL_METASTASIS_SITE = "metastasis_site"
COL_TUMOR_TREATMENT = "tumor_treatment"
COL_CANCER_SITE_ORIGIN = "cancer_site_origin"
COL_TUMOUR_GRADE = "tumour_grade"
COL_TUMOUR_STAGE = "tumour_stage"
COL_HISTOLOGICAL_SUBTYPE = "histological_subtype"
COL_PATIENT_ETHNICITY = "patient_ethnicity"


# -----------------------------
# 3. Manually curated sample metadata
# -----------------------------
SAMPLE_METADATA = {
    "GSM4416534": {
        COL_PATIENT_ID: "PT-3232",
        COL_DATASET_ID: "GSE147082",
        COL_CANCER_TYPE: "Ovarian",
        COL_TUMOR_SITE: "Metastasis",
        COL_METASTASIS_SITE: "Omentum",
        COL_TUMOR_TREATMENT: "No",
        COL_CANCER_SITE_ORIGIN: "Left fallopian",
        COL_TUMOUR_GRADE: "HGSC",
        COL_TUMOUR_STAGE: "IIIc",
        COL_HISTOLOGICAL_SUBTYPE: "Serous",
        COL_PATIENT_ETHNICITY: "Asian",
    },
    "GSM4416535": {
        COL_PATIENT_ID: "PT-5150",
        COL_DATASET_ID: "GSE147082",
        COL_CANCER_TYPE: "Ovarian",
        COL_TUMOR_SITE: "Metastasis",
        COL_METASTASIS_SITE: "Omentum",
        COL_TUMOR_TREATMENT: "Yes",
        COL_CANCER_SITE_ORIGIN: "Unknown",
        COL_TUMOUR_GRADE: "HGSC",
        COL_TUMOUR_STAGE: "IVb",
        COL_HISTOLOGICAL_SUBTYPE: "Serous",
        COL_PATIENT_ETHNICITY: "White",
    },
    "GSM4416536": {
        COL_PATIENT_ID: "PT-6885",
        COL_DATASET_ID: "GSE147082",
        COL_CANCER_TYPE: "Ovarian",
        COL_TUMOR_SITE: "Metastasis",
        COL_METASTASIS_SITE: "Omentum",
        COL_TUMOR_TREATMENT: "Yes",
        COL_CANCER_SITE_ORIGIN: "Left fallopian",
        COL_TUMOUR_GRADE: "HGSC",
        COL_TUMOUR_STAGE: "IIIc",
        COL_HISTOLOGICAL_SUBTYPE: "Serous",
        COL_PATIENT_ETHNICITY: "Black",
    },
    "GSM4416537": {
        COL_PATIENT_ID: "PT-4806",
        COL_DATASET_ID: "GSE147082",
        COL_CANCER_TYPE: "Ovarian",
        COL_TUMOR_SITE: "Metastasis",
        COL_METASTASIS_SITE: "Omentum",
        COL_TUMOR_TREATMENT: "Yes",
        COL_CANCER_SITE_ORIGIN: "Left fallopian",
        COL_TUMOUR_GRADE: "HGSC",
        COL_TUMOUR_STAGE: "IIIc",
        COL_HISTOLOGICAL_SUBTYPE: "Serous",
        COL_PATIENT_ETHNICITY: "Black",
    },
    "GSM4416538": {
        COL_PATIENT_ID: "PT-3401",
        COL_DATASET_ID: "GSE147082",
        COL_CANCER_TYPE: "Ovarian",
        COL_TUMOR_SITE: "Metastasis",
        COL_METASTASIS_SITE: "Omentum",
        COL_TUMOR_TREATMENT: "No",
        COL_CANCER_SITE_ORIGIN: "Left fallopian",
        COL_TUMOUR_GRADE: "HGSC",
        COL_TUMOUR_STAGE: "IIIc",
        COL_HISTOLOGICAL_SUBTYPE: "Serous",
        COL_PATIENT_ETHNICITY: "White",
    },
    "GSM4416539": {
        COL_PATIENT_ID: "PT-2834",
        COL_DATASET_ID: "GSE147082",
        COL_CANCER_TYPE: "Ovarian",
        COL_TUMOR_SITE: "Metastasis",
        COL_METASTASIS_SITE: "Omentum",
        COL_TUMOR_TREATMENT: "No",
        COL_CANCER_SITE_ORIGIN: "Fallopian",
        COL_TUMOUR_GRADE: "HGSC",
        COL_TUMOUR_STAGE: "IIIc",
        COL_HISTOLOGICAL_SUBTYPE: "Serous",
        COL_PATIENT_ETHNICITY: "Asian",
    },
}


# -----------------------------
# 4. Directory setup
# -----------------------------
def ensure_directories() -> None:
    """Create all required directories for this case."""
    os.makedirs(RAW_DIR, exist_ok=True)
    os.makedirs(ANNOTATED_DIR, exist_ok=True)
    os.makedirs(OUTPUT_DIR, exist_ok=True)


# -----------------------------
# 5. Download GEO archive
# -----------------------------
def download_raw_tar() -> None:
    """
    Download the GEO supplementary RAW tar archive if it is not already present.
    """
    if os.path.exists(RAW_TAR_PATH):
        print(f"RAW tar already exists: {RAW_TAR_PATH}")
        return

    print("Downloading GEO supplementary archive...")
    print("URL:", GEO_DOWNLOAD_URL)
    print("Output:", RAW_TAR_PATH)

    urlretrieve(GEO_DOWNLOAD_URL, RAW_TAR_PATH)

    print("Download completed.")


# -----------------------------
# 6. Extract GEO archive
# -----------------------------
def extract_raw_tar() -> None:
    """
    Extract the GEO RAW tar archive into RAW_DIR.
    """
    extracted_csv_gz = [
        f for f in os.listdir(RAW_DIR)
        if f.startswith("GSM") and f.endswith(".csv.gz")
    ]

    if extracted_csv_gz:
        print(f"RAW archive already appears extracted. Found {len(extracted_csv_gz)} GSM csv.gz files.")
        return

    print("Extracting RAW tar archive...")

    with tarfile.open(RAW_TAR_PATH, "r") as tar:
        tar.extractall(RAW_DIR)

    print("Extraction completed.")


# -----------------------------
# 7. List sample files
# -----------------------------
def list_gsm_csv_files() -> list[str]:
    """
    List GSM-level compressed CSV count tables.
    """
    sample_files = sorted([
        os.path.join(RAW_DIR, f)
        for f in os.listdir(RAW_DIR)
        if f.startswith("GSM") and f.endswith(".csv.gz")
    ])

    print("\n=== Detected GSM CSV files ===")
    print("Number of files:", len(sample_files))

    for fp in sample_files:
        print("-", os.path.basename(fp))

    return sample_files


# -----------------------------
# 8. Parse sample ID and patient ID from filename
# -----------------------------
def parse_ids_from_filename(filepath: str) -> tuple[str, str]:
    """
    Extract sample_id and patient_id from filenames such as:
        GSM4416534_PT-3232.csv.gz
    """
    basename = os.path.basename(filepath)
    basename = basename.replace(".csv.gz", "")

    parts = basename.split("_")

    if len(parts) < 2:
        raise ValueError(f"Could not parse sample/patient IDs from filename: {basename}")

    sample_id = parts[0].strip()
    patient_id_from_filename = parts[1].strip()

    return sample_id, patient_id_from_filename


# -----------------------------
# 9. Inspect first CSV
# -----------------------------
def inspect_first_csv(sample_files: list[str]) -> None:
    """
    Inspect the first compressed CSV file to infer raw table structure.
    """
    if not sample_files:
        raise FileNotFoundError("No GSM csv.gz files found. Cannot inspect input structure.")

    first_file = sample_files[0]

    print("\n=== Inspecting first CSV file ===")
    print("File:", first_file)

    df_preview = pd.read_csv(first_file, compression="gzip", nrows=5)

    print("\nPreview shape:")
    print(df_preview.shape)

    print("\nFirst columns:")
    print(list(df_preview.columns[:10]))

    print("\nPreview rows:")
    print(df_preview.head())

    print("\nInterpretation:")
    print("- The first column contains gene symbols.")
    print("- The remaining columns correspond to cells.")
    print("- This indicates raw matrix orientation: genes x cells.")
    print("- AnnData expects cells x genes.")
    print("- Therefore transpose is required: AnnData(df.T).")


# -----------------------------
# 10. Validate raw count dataframe
# -----------------------------
def validate_count_dataframe(df: pd.DataFrame, filepath: str) -> None:
    """
    Basic sanity checks before AnnData creation.
    """
    if df.empty:
        raise ValueError(f"Input dataframe is empty: {filepath}")

    if df.shape[1] < 2:
        raise ValueError(
            f"Input dataframe should contain one gene column and at least one cell column: {filepath}"
        )

    first_col = df.columns[0]
    if first_col is None:
        raise ValueError(f"First column name is missing in: {filepath}")

    duplicated_genes = df.iloc[:, 0].duplicated().sum()
    if duplicated_genes > 0:
        print(f"Warning: {duplicated_genes} duplicated gene names detected in {os.path.basename(filepath)}")

    numeric_part = df.iloc[:, 1:]
    non_numeric_columns = numeric_part.columns[
        ~numeric_part.apply(lambda col: pd.api.types.is_numeric_dtype(col))
    ]

    if len(non_numeric_columns) > 0:
        raise ValueError(
            f"Non-numeric cell count columns detected in {filepath}: {list(non_numeric_columns[:10])}"
        )


# -----------------------------
# 11. Build one AnnData object from one CSV.gz
# -----------------------------
def build_adata_from_count_csv(filepath: str) -> AnnData:
    """
    Build one AnnData object from one compressed CSV count table.

    Raw file:
        rows    = genes
        columns = cells

    AnnData:
        rows    = cells
        columns = genes

    Therefore:
        AnnData(df.T)
    """
    sample_id, patient_id_from_filename = parse_ids_from_filename(filepath)

    if sample_id not in SAMPLE_METADATA:
        raise KeyError(f"{sample_id} is missing from SAMPLE_METADATA dictionary.")

    sample_metadata = SAMPLE_METADATA[sample_id]
    patient_id_from_metadata = sample_metadata[COL_PATIENT_ID]

    if patient_id_from_filename != patient_id_from_metadata:
        print(
            "Warning: patient ID mismatch:",
            sample_id,
            "filename =",
            patient_id_from_filename,
            "metadata =",
            patient_id_from_metadata,
        )

    print(f"\n--- Building AnnData for {sample_id} / {patient_id_from_metadata} ---")

    # 1. Read raw CSV.gz
    df = pd.read_csv(filepath, compression="gzip")
    print("Raw dataframe shape:", df.shape)

    # 2. Validate before transformation
    validate_count_dataframe(df, filepath)

    # 3. Rename first column as gene
    original_first_column = df.columns[0]
    df = df.rename(columns={original_first_column: "gene"})

    # 4. Set gene names as index
    df = df.set_index("gene")

    # 5. Make index and columns string-safe
    df.index = df.index.astype(str)
    df.columns = df.columns.astype(str)

    # 6. Report interpreted structure
    print("Interpreted raw matrix orientation: genes x cells")
    print("Number of genes:", df.shape[0])
    print("Number of cells:", df.shape[1])
    print("First 5 genes:", list(df.index[:5]))
    print("First 5 cells:", list(df.columns[:5]))

    # 7. Create AnnData
    # AnnData expects observations/cells as rows and variables/genes as columns.
    # Since df is genes x cells, df.T is cells x genes.
    adata = AnnData(df.T)

    # 8. Assign names explicitly
    adata.obs_names = df.columns.astype(str)
    adata.var_names = df.index.astype(str)

    # 9. Add gene-level annotation
    adata.var["gene_symbol"] = adata.var_names

    # 10. Add sample-level metadata to every cell
    adata.obs[COL_SAMPLE_ID] = sample_id

    for column_name, value in sample_metadata.items():
        adata.obs[column_name] = value

    # 11. Add provenance columns
    adata.obs["source_file"] = os.path.basename(filepath)
    adata.obs["paper_id"] = PAPER_ID
    adata.obs["geo_accession"] = GEO_ACCESSION

    # 12. Final per-sample checks
    print("AnnData shape:", adata.shape)
    print("obs columns:", list(adata.obs.columns))
    print("var columns:", list(adata.var.columns))

    return adata


# -----------------------------
# 12. Build all sample AnnData objects
# -----------------------------
def build_all_sample_adatas(sample_files: list[str]) -> list[AnnData]:
    """
    Build one AnnData object per GSM CSV file.
    """
    list_adata = []

    for filepath in sample_files:
        adata_sample = build_adata_from_count_csv(filepath)
        list_adata.append(adata_sample)

    print("\n=== Built sample AnnData objects ===")
    print("Number of AnnData objects:", len(list_adata))

    for adata_sample in list_adata:
        sample_id = adata_sample.obs[COL_SAMPLE_ID].iloc[0]
        patient_id = adata_sample.obs[COL_PATIENT_ID].iloc[0]
        print(f"- {sample_id} / {patient_id}: {adata_sample.n_obs} cells x {adata_sample.n_vars} genes")

    return list_adata


# -----------------------------
# 13. Check NaN values
# -----------------------------
def check_nan_in_adata(adata_obj: AnnData, label: str) -> None:
    """
    Check NaN values in AnnData X, obs, and var.
    """
    print(f"\n=== NaN check: {label} ===")

    if issparse(adata_obj.X):
        nan_count_x = np.isnan(adata_obj.X.data).sum()
    else:
        nan_count_x = np.isnan(adata_obj.X).sum()

    print("NaN count in X:", nan_count_x)
    print("NaN exists in obs:", adata_obj.obs.isna().values.any())
    print("NaN exists in var:", adata_obj.var.isna().values.any())



# -----------------------------
# 14. Align genes across samples
# -----------------------------
def align_common_genes(list_adata: list[AnnData]) -> list[AnnData]:
    """
    Keep only genes shared across all samples.

    Why:
        Each sample CSV may not contain exactly the same gene set.
        If we concatenate with different gene sets, missing values can appear.
        To build a clean combined object, we first restrict all samples to
        the intersection of genes.
    """
    if not list_adata:
        raise ValueError("No AnnData objects provided.")

    common_genes = list_adata[0].var_names

    for adata_sample in list_adata[1:]:
        common_genes = common_genes.intersection(adata_sample.var_names)

    print("\n=== Common gene alignment ===")
    print("Genes in first sample:", list_adata[0].n_vars)
    print("Common genes across all samples:", len(common_genes))

    aligned = []

    for adata_sample in list_adata:
        sample_id = adata_sample.obs[COL_SAMPLE_ID].iloc[0]
        aligned_sample = adata_sample[:, common_genes].copy()
        aligned.append(aligned_sample)
        print(f"- {sample_id}: {aligned_sample.n_obs} cells x {aligned_sample.n_vars} common genes")

    return aligned


# -----------------------------
# 15. Concatenate samples
# -----------------------------
def concatenate_samples(list_adata: list[AnnData]) -> AnnData:
    """
    Concatenate all sample-level AnnData objects into one combined AnnData.

    axis=0 means:
        stack cells from all samples on top of each other.

    join="inner" means:
        keep only variables/genes shared across all samples.
    """
    print("\n=== Concatenating samples ===")

    adata_combined = ad.concat(
        list_adata,
        axis=0,
        join="inner",
        index_unique="-"
    )

    adata_combined.obs_names_make_unique()

    # Re-add gene-level annotation after concatenation.
    # anndata.concat may drop var columns depending on merge behavior.
    adata_combined.var["gene_symbol"] = adata_combined.var_names.astype(str)

    print("Combined AnnData shape:", adata_combined.shape)
    print("Number of cells:", adata_combined.n_obs)
    print("Number of genes:", adata_combined.n_vars)
    print("obs columns:", list(adata_combined.obs.columns))
    print("var columns:", list(adata_combined.var.columns))

    return adata_combined


# -----------------------------
# 16. Write final h5ad
# -----------------------------
def write_final_h5ad(adata_combined: AnnData) -> None:
    """
    Save the standardized combined AnnData object as h5ad.
    """
    print("\n=== Writing final h5ad ===")
    adata_combined.write_h5ad(OUTPUT_H5AD_PATH)
    print("Saved:", OUTPUT_H5AD_PATH)


# -----------------------------
# 17. Clean raw files
# -----------------------------
def cleanup_raw_files() -> None:
    """
    Remove downloaded raw files from the local computer.

    This keeps the project lightweight:
        - raw TAR is deleted
        - extracted GSM csv.gz files are deleted

    The workflow remains reproducible because the script can download
    the GEO archive again when needed.
    """
    print("\n=== Cleaning raw downloaded files ===")

    if not os.path.exists(RAW_DIR):
        print("Raw directory does not exist. Nothing to clean.")
        return

    removed = 0

    for filename in os.listdir(RAW_DIR):
        filepath = os.path.join(RAW_DIR, filename)

        if os.path.isfile(filepath):
            os.remove(filepath)
            removed += 1

    print(f"Removed {removed} raw files from: {RAW_DIR}")


# -----------------------------
# 18. Final report
# -----------------------------
def print_final_report(adata_combined: AnnData, keep_raw: bool) -> None:
    """
    Print a human-readable summary of the standardization case.
    """
    print("\n=== Final case report ===")
    print("Case:", PAPER_ID)
    print("GEO accession:", GEO_ACCESSION)
    print("Detected input format: compressed CSV count tables")
    print("Raw orientation: genes x cells")
    print("AnnData conversion: AnnData(df.T)")
    print("Metadata strategy: manually curated sample-level GEO metadata")
    print("Final object:", adata_combined)
    print("Final h5ad:", OUTPUT_H5AD_PATH)

    if keep_raw:
        print("Raw files were kept because --keep-raw was used.")
    else:
        print("Raw files were deleted after processing.")


# -----------------------------
# 19. Main runner
# -----------------------------
def run_case_001(keep_raw: bool = False) -> None:
    """
    Run the full Case 001 workflow.

    Default behavior:
        raw downloaded files are deleted after processing.

    Debug behavior:
        use --keep-raw to keep raw files locally.
    """
    print("Case:", PAPER_ID)
    print("GEO accession:", GEO_ACCESSION)
    print("Raw directory:", RAW_DIR)
    print("Output h5ad path:", OUTPUT_H5AD_PATH)
    print("Number of manually curated samples:", len(SAMPLE_METADATA))
    print("Keep raw files:", keep_raw)

    ensure_directories()
    download_raw_tar()
    extract_raw_tar()

    sample_files = list_gsm_csv_files()
    inspect_first_csv(sample_files)

    list_adata = build_all_sample_adatas(sample_files)

    for adata_sample in list_adata:
        sample_id = adata_sample.obs[COL_SAMPLE_ID].iloc[0]
        check_nan_in_adata(adata_sample, label=sample_id)

    list_adata_common = align_common_genes(list_adata)
    adata_combined = concatenate_samples(list_adata_common)

    check_nan_in_adata(adata_combined, label="combined")
    write_final_h5ad(adata_combined)

    if not keep_raw:
        cleanup_raw_files()

    print_final_report(adata_combined, keep_raw=keep_raw)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Build standardized AnnData/h5ad for Case 001: 2021_Olelekan / GSE147082"
    )
    parser.add_argument(
        "--keep-raw",
        action="store_true",
        help="Keep downloaded raw TAR and extracted CSV.gz files after processing."
    )

    args = parser.parse_args()
    run_case_001(keep_raw=args.keep_raw)
