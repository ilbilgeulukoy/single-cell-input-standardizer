import os
import argparse
from collections import defaultdict


PAPER_ID = "2021_Regner"
GEO_ACCESSION = "GSE173682"

BASE_DIR = os.path.join("data", "case_002_2021_regner")
RAW_DIR = os.path.join(BASE_DIR, "raw")
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")

# GEO reports the full supplementary TAR archive as approximately 15.8 GB.
# For this project, the local personal-computer workflow should not download it by default.
GEO_ARCHIVE_NAME = "GSE173682_RAW.tar"
GEO_ARCHIVE_SIZE = "15.8 GB"
GEO_DOWNLOAD_URL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE173682&format=file"


# This is a known expected file pattern from the user's previous working script.
# The full archive contains per-GSM 10x-like RNA Matrix Market triplets:
#   - matrix.mtx.gz
#   - features.tsv.gz
#   - barcodes.tsv.gz
# and ATAC fragment files for matched ATAC samples.
EXPECTED_RNA_GSM_IDS = [
    "GSM5276933",
    "GSM5276934",
    "GSM5276935",
    "GSM5276936",
    "GSM5276937",
    "GSM5276938",
    "GSM5276939",
    "GSM5276940",
    "GSM5276941",
    "GSM5276943",
]

SAMPLE_METADATA = {
    "GSM5276933": {
        "patient_id": "3533EL",
        "dataset_id": "GSE173682",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "Unknown",
        "cancer_site_origin": "Endometrium",
        "tumour_grade": "2",
        "tumour_stage": "IA",
        "histological_subtype": "Endometrioid",
        "patient_ethnicity": "African American",
    },
    "GSM5276934": {
        "patient_id": "3571DL",
        "dataset_id": "GSE173682",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "Unknown",
        "cancer_site_origin": "Endometrium",
        "tumour_grade": "3",
        "tumour_stage": "IA",
        "histological_subtype": "Endometrioid",
        "patient_ethnicity": "Caucasian",
    },
    "GSM5276935": {
        "patient_id": "36186L",
        "dataset_id": "GSE173682",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "Unknown",
        "cancer_site_origin": "Endometrium",
        "tumour_grade": "2",
        "tumour_stage": "IA",
        "histological_subtype": "Endometrioid",
        "patient_ethnicity": "Caucasian",
    },
    "GSM5276936": {
        "patient_id": "36639L",
        "dataset_id": "GSE173682",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "Unknown",
        "cancer_site_origin": "Endometrium",
        "tumour_grade": "1",
        "tumour_stage": "IA",
        "histological_subtype": "Endometrioid",
        "patient_ethnicity": "Caucasian",
    },
    "GSM5276937": {
        "patient_id": "366C5L",
        "dataset_id": "GSE173682",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "Unknown",
        "cancer_site_origin": "Endometrium",
        "tumour_grade": "1",
        "tumour_stage": "IA",
        "histological_subtype": "Endometrioid",
        "patient_ethnicity": "Caucasian",
    },
    "GSM5276938": {
        "patient_id": "37EACL",
        "dataset_id": "GSE173682",
        "cancer_type": "Ovarian",
        "tumor_site": "Metastasis",
        "metastasis_site": "Unknown",
        "tumor_treatment": "Unknown",
        "cancer_site_origin": "Ovary",
        "tumour_grade": "HGSC",
        "tumour_stage": "IIIA",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Caucasian",
    },
    "GSM5276939": {
        "patient_id": "38FE7L",
        "dataset_id": "GSE173682",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "Unknown",
        "cancer_site_origin": "Endometrium",
        "tumour_grade": "1",
        "tumour_stage": "IA",
        "histological_subtype": "Endometrioid",
        "patient_ethnicity": "Caucasian",
    },
    "GSM5276940": {
        "patient_id": "3BAE2L",
        "dataset_id": "GSE173682",
        "cancer_type": "Ovarian",
        "tumor_site": "Metastasis",
        "metastasis_site": "Unknown",
        "tumor_treatment": "Unknown",
        "cancer_site_origin": "Ovary",
        "tumour_grade": "HGSC",
        "tumour_stage": "IIB",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Caucasian",
    },
    "GSM5276941": {
        "patient_id": "3E5CFL",
        "dataset_id": "GSE173682",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "Unknown",
        "cancer_site_origin": "Ovary",
        "tumour_grade": "HGSC",
        "tumour_stage": "IIIC",
        "histological_subtype": "Serous",
        "patient_ethnicity": "Asian",
    },
    "GSM5276943": {
        "patient_id": "3E5CFL",
        "dataset_id": "GSE173682",
        "cancer_type": "Ovarian",
        "tumor_site": "Primary",
        "metastasis_site": "Unknown",
        "tumor_treatment": "Unknown",
        "cancer_site_origin": "Ovary",
        "tumour_grade": "Unknown",
        "tumour_stage": "IVB",
        "histological_subtype": "Carcinosarcoma",
        "patient_ethnicity": "Asian",
    },
}


def ensure_directories() -> None:
    os.makedirs(RAW_DIR, exist_ok=True)
    os.makedirs(OUTPUT_DIR, exist_ok=True)


def build_expected_file_manifest() -> list[dict[str, str]]:
    """
    Build a lightweight expected manifest without downloading the 15.8 GB GEO archive.

    The exact filenames inside GEO follow a GSM-based naming pattern in the user's
    previous working script:
        GSM*_matrix-*.mtx.gz
        GSM*_features-*.tsv.gz
        GSM*_barcodes-*.tsv.gz
        GSM*_ATAC_fragments*.tsv.gz
    """
    manifest = []

    for gsm in EXPECTED_RNA_GSM_IDS:
        metadata = SAMPLE_METADATA[gsm]

        manifest.append({
            "gsm_id": gsm,
            "patient_id": metadata["patient_id"],
            "expected_role": "matrix",
            "expected_pattern": f"{gsm}_*_matrix-*.mtx.gz",
            "modality": "RNA",
        })
        manifest.append({
            "gsm_id": gsm,
            "patient_id": metadata["patient_id"],
            "expected_role": "features",
            "expected_pattern": f"{gsm}_*_features-*.tsv.gz",
            "modality": "RNA",
        })
        manifest.append({
            "gsm_id": gsm,
            "patient_id": metadata["patient_id"],
            "expected_role": "barcodes",
            "expected_pattern": f"{gsm}_*_barcodes-*.tsv.gz",
            "modality": "RNA",
        })

    return manifest


def print_case_002_inspection_report() -> None:
    manifest = build_expected_file_manifest()

    print("Case:", PAPER_ID)
    print("GEO accession:", GEO_ACCESSION)
    print("GEO archive:", GEO_ARCHIVE_NAME)
    print("GEO archive size:", GEO_ARCHIVE_SIZE)
    print("Download URL:", GEO_DOWNLOAD_URL)

    print("\n=== Data-size-aware decision ===")
    print("The full GEO supplementary archive is too large for the default local workflow.")
    print("Recommended mode: server/HPC processing or manifest-based inspection.")
    print("Local raw download: skipped by default.")

    print("\n=== Detected expected input pattern ===")
    print("Format: per-sample 10x-like Matrix Market triplets")
    print("Required RNA files per sample:")
    print("- matrix.mtx.gz")
    print("- features.tsv.gz")
    print("- barcodes.tsv.gz")
    print("Additional files:")
    print("- ATAC fragments may be present but are not processed in the initial RNA AnnData workflow.")

    print("\n=== Expected RNA samples ===")
    for gsm in EXPECTED_RNA_GSM_IDS:
        metadata = SAMPLE_METADATA[gsm]
        print(
            f"- {gsm} / patient {metadata['patient_id']} / "
            f"{metadata['histological_subtype']} / stage {metadata['tumour_stage']}"
        )

    print("\n=== Lightweight expected manifest ===")
    for row in manifest:
        print(
            f"{row['gsm_id']},{row['patient_id']},"
            f"{row['expected_role']},{row['expected_pattern']},{row['modality']}"
        )

    print("\n=== Recommended AnnData strategy ===")
    print("For each GSM sample:")
    print("1. Locate matrix, features, and barcodes files.")
    print("2. Read matrix.mtx.gz as a sparse Matrix Market file.")
    print("3. Transpose matrix to cells x genes.")
    print("4. Read barcodes.tsv.gz and assign obs_names.")
    print("5. Prefix cell barcodes with sample_id to avoid duplicated cell names.")
    print("6. Read features.tsv.gz and assign gene symbols as var_names.")
    print("7. Merge duplicated gene symbols by summing counts.")
    print("8. Add sample-level metadata to adata.obs.")
    print("9. Align common genes across samples.")
    print("10. Concatenate all sample AnnData objects.")
    print("11. Validate X, obs, and var.")
    print("12. Write standardized h5ad.")

    print("\n=== Tool design implication ===")
    print("This case should create a reusable 10x-like Matrix Market standardizer module.")
    print("Suggested module: src/tenx_mtx_standardizer.py")


def run_case_002_inspection() -> None:
    ensure_directories()
    print_case_002_inspection_report()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Inspect Case 002 without downloading the full 15.8 GB GEO archive."
    )
    parser.parse_args()
    run_case_002_inspection()
