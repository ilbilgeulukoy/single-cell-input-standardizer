import os
import tarfile
import argparse
from urllib.request import urlretrieve

import h5py
import pandas as pd


PAPER_ID = "2022_Loret"
GEO_ACCESSION = "GSE201047"

BASE_DIR = os.path.join("data", "case_017_2022_loret")
RAW_DIR = os.path.join(BASE_DIR, "raw")

RAW_TAR_FILENAME = "GSE201047_RAW.tar"
RAW_TAR_PATH = os.path.join(RAW_DIR, RAW_TAR_FILENAME)

GEO_DOWNLOAD_URL = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE201047&format=file"


SAMPLE_METADATA = {
    "GSM6049610": {"patient_id": "N_OT_PT1", "tumor_site": "Metastasis", "metastasis_site": "Omentum", "tumor_treatment": "No", "cancer_site_origin": "Omentum", "tumour_stage": "IIIc"},
    "GSM6049611": {"patient_id": "N_A_PT1", "tumor_site": "Metastasis", "metastasis_site": "Ascites", "tumor_treatment": "No", "cancer_site_origin": "Omentum", "tumour_stage": "IIIc"},
    "GSM6049612": {"patient_id": "N_PER_PT1", "tumor_site": "Metastasis", "metastasis_site": "Peritoneum", "tumor_treatment": "No", "cancer_site_origin": "Omentum", "tumour_stage": "IIIc"},
    "GSM6049613": {"patient_id": "N_OM_PT1", "tumor_site": "Metastasis", "metastasis_site": "Omentum", "tumor_treatment": "No", "cancer_site_origin": "Omentum", "tumour_stage": "IIIc"},
    "GSM6049614": {"patient_id": "N_BL_PT1", "tumor_site": "Primary", "metastasis_site": "Bladder", "tumor_treatment": "No", "cancer_site_origin": "Ovary", "tumour_stage": "IIIc"},
    "GSM6049615": {"patient_id": "T_OT_PT1", "tumor_site": "Primary", "metastasis_site": "Omentum", "tumor_treatment": "Yes", "cancer_site_origin": "Ovary", "tumour_stage": "IIIc"},
    "GSM6049616": {"patient_id": "T_PER_PT1", "tumor_site": "Primary", "metastasis_site": "Peritoneum", "tumor_treatment": "Yes", "cancer_site_origin": "Ovary", "tumour_stage": "IIIc"},
    "GSM6049617": {"patient_id": "T_OM_PT1", "tumor_site": "Primary", "metastasis_site": "Omentum", "tumor_treatment": "Yes", "cancer_site_origin": "Ovary", "tumour_stage": "IIIc"},
    "GSM6049618": {"patient_id": "T_A_PT1", "tumor_site": "Primary", "metastasis_site": "Ascites", "tumor_treatment": "Yes", "cancer_site_origin": "Ovary", "tumour_stage": "IIIc"},
    "GSM6049619": {"patient_id": "N_OT_PT2", "tumor_site": "Primary", "metastasis_site": "Omentum", "tumor_treatment": "No", "cancer_site_origin": "Ovary", "tumour_stage": "IIIc"},
    "GSM6049620": {"patient_id": "N_A_PT2", "tumor_site": "Primary", "metastasis_site": "Ascites", "tumor_treatment": "No", "cancer_site_origin": "Ovary", "tumour_stage": "IIIc"},
    "GSM6049621": {"patient_id": "N_OM_PT2", "tumor_site": "Primary", "metastasis_site": "Omentum", "tumor_treatment": "No", "cancer_site_origin": "Ovary", "tumour_stage": "IIIc"},
    "GSM6049622": {"patient_id": "T_OT_PT2", "tumor_site": "Metastasis", "metastasis_site": "Omentum", "tumor_treatment": "Yes", "cancer_site_origin": "Omentum", "tumour_stage": "IIIc"},
    "GSM6049623": {"patient_id": "T_OM_PT2", "tumor_site": "Metastasis", "metastasis_site": "Omentum", "tumor_treatment": "Yes", "cancer_site_origin": "Omentum", "tumour_stage": "IIIc"},
    "GSM6049624": {"patient_id": "T_PER_PT2", "tumor_site": "Primary", "metastasis_site": "Peritoneum", "tumor_treatment": "Yes", "cancer_site_origin": "Ovary", "tumour_stage": "IIIc"},
    "GSM6049625": {"patient_id": "N_OT_PT3", "tumor_site": "Primary", "metastasis_site": "Omentum", "tumor_treatment": "No", "cancer_site_origin": "Ovary", "tumour_stage": "IIIb"},
    "GSM6049626": {"patient_id": "N_A_PT3", "tumor_site": "Primary", "metastasis_site": "Ascites", "tumor_treatment": "No", "cancer_site_origin": "Ovary", "tumour_stage": "IIIb"},
    "GSM6049627": {"patient_id": "N_PER_PT3", "tumor_site": "Primary", "metastasis_site": "Peritoneum", "tumor_treatment": "No", "cancer_site_origin": "Ovary", "tumour_stage": "IIIb"},
    "GSM6049628": {"patient_id": "T_OT_PT3", "tumor_site": "Primary", "metastasis_site": "Omentum", "tumor_treatment": "Yes", "cancer_site_origin": "Ovary", "tumour_stage": "IIIb"},
    "GSM6049629": {"patient_id": "T_A_PT3", "tumor_site": "Primary", "metastasis_site": "Ascites", "tumor_treatment": "Yes", "cancer_site_origin": "Ovary", "tumour_stage": "IIIb"},
    "GSM6049630": {"patient_id": "T_OM_PT3", "tumor_site": "Primary", "metastasis_site": "Omentum", "tumor_treatment": "Yes", "cancer_site_origin": "Ovary", "tumour_stage": "IIIb"},
    "GSM6049631": {"patient_id": "T_PER_PT3", "tumor_site": "Primary", "metastasis_site": "Peritoneum", "tumor_treatment": "Yes", "cancer_site_origin": "Ovary", "tumour_stage": "IIIb"},
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


def parse_filename_tokens(filename):
    clean = filename.replace(".h5", "")
    parts = clean.split("_")
    return parts


def inspect_10x_h5(filename):
    filepath = os.path.join(RAW_DIR, filename)
    gsm = filename.split("_")[0].strip()

    print("\n=== Inspecting 10x h5 file ===")
    print("File:", filename)
    print("GSM:", gsm)

    tokens = parse_filename_tokens(filename)
    print("Filename tokens:", tokens)

    meta = SAMPLE_METADATA.get(gsm)

    if meta:
        print("Curated patient_id:", meta["patient_id"])
        print("Curated tumor_site:", meta["tumor_site"])
        print("Curated metastasis_site:", meta["metastasis_site"])
        print("Curated tumor_treatment:", meta["tumor_treatment"])
        print("Curated tumour_stage:", meta["tumour_stage"])
    else:
        print("No curated metadata for this GSM.")

    with h5py.File(filepath, "r") as f:
        print("Top-level keys:", list(f.keys()))

        duplicated_names = None
        first_features = []

        if "matrix" in f:
            matrix = f["matrix"]
            format_version = "10x_h5_v3_matrix_group"

            print("Detected h5 layout:", format_version)
            print("matrix subkeys:", list(matrix.keys()))

            shape = tuple(int(x) for x in matrix["shape"][:])
            n_genes = shape[0]
            n_cells = shape[1]

            print("10x matrix shape from /matrix/shape:", shape)
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

                ids = decode_array(features["id"], 5) if "id" in features else []
                names = decode_array(features["name"], 5) if "name" in features else []
                feature_types = decode_array(features["feature_type"], 5) if "feature_type" in features else []

                print("First feature IDs:", ids)
                print("First feature names:", names)
                print("First feature types:", feature_types)

                if "name" in features:
                    all_names = [
                        x.decode("utf-8") if isinstance(x, bytes) else str(x)
                        for x in features["name"][:]
                    ]
                    duplicated_names = pd.Index(all_names).duplicated().sum()
                    print("Duplicated feature names:", duplicated_names)

                for i in range(min(len(ids), len(names), len(feature_types))):
                    first_features.append(f"{ids[i]} / {names[i]} / {feature_types[i]}")

        else:
            genome_keys = [k for k in f.keys() if isinstance(f[k], h5py.Group)]
            if not genome_keys:
                print("WARNING: no /matrix group and no genome group found.")
                return None

            genome_key = genome_keys[0]
            matrix = f[genome_key]
            format_version = "10x_h5_v2_genome_group"

            print("Detected h5 layout:", format_version)
            print("Genome group:", genome_key)
            print("genome group subkeys:", list(matrix.keys()))

            required = ["data", "indices", "indptr", "shape", "barcodes", "genes", "gene_names"]
            missing = [x for x in required if x not in matrix]
            print("Missing old-format required keys:", missing)

            if missing:
                return None

            shape = tuple(int(x) for x in matrix["shape"][:])
            n_genes = shape[0]
            n_cells = shape[1]

            print("10x matrix shape from genome/shape:", shape)
            print("Interpreted raw orientation: genes x cells")
            print("AnnData expected shape after read_10x_h5:", (n_cells, n_genes))

            print("data shape:", matrix["data"].shape)
            print("indices shape:", matrix["indices"].shape)
            print("indptr shape:", matrix["indptr"].shape)
            print("barcodes shape:", matrix["barcodes"].shape)
            print("genes shape:", matrix["genes"].shape)
            print("gene_names shape:", matrix["gene_names"].shape)

            print("First barcodes:", decode_array(matrix["barcodes"], 5))
            print("First gene IDs:", decode_array(matrix["genes"], 5))
            print("First gene names:", decode_array(matrix["gene_names"], 5))

            all_names = [
                x.decode("utf-8") if isinstance(x, bytes) else str(x)
                for x in matrix["gene_names"][:]
            ]
            duplicated_names = pd.Index(all_names).duplicated().sum()
            print("Duplicated gene names:", duplicated_names)

            ids = decode_array(matrix["genes"], 5)
            names = decode_array(matrix["gene_names"], 5)
            for i in range(min(len(ids), len(names))):
                first_features.append(f"{ids[i]} / {names[i]}")

    return {
        "gsm": gsm,
        "filename": filename,
        "n_cells": n_cells,
        "n_genes": n_genes,
        "duplicated_feature_names": duplicated_names,
        "tokens": tokens,
        "metadata": meta,
        "first_features": first_features,
    }


def inspect_all_h5_files(h5_files):
    print("\n=== All h5 sample summaries ===")

    summaries = []

    for filename in h5_files:
        summary = inspect_10x_h5(filename)
        if summary:
            summaries.append(summary)

    print("\n=== Compact summary table ===")
    total_cells = 0

    for s in summaries:
        total_cells += s["n_cells"]

        meta = s["metadata"] or {}

        print(
            s["gsm"],
            "patient_id=", meta.get("patient_id"),
            "cells=", s["n_cells"],
            "genes=", s["n_genes"],
            "tumor_site=", meta.get("tumor_site"),
            "metastasis_site=", meta.get("metastasis_site"),
            "treatment=", meta.get("tumor_treatment"),
            "stage=", meta.get("tumour_stage"),
            "dup_gene_names=", s["duplicated_feature_names"],
        )

    print("\nTotal cells across h5 files:", total_cells)
    print("Unique gene counts:", sorted(set(s["n_genes"] for s in summaries)))
    print("Number of h5 files inspected:", len(summaries))

    print("\n=== Metadata distribution summary ===")

    for field in ["tumor_site", "metastasis_site", "tumor_treatment", "cancer_site_origin", "tumour_stage"]:
        counts = {}
        for s in summaries:
            meta = s["metadata"] or {}
            value = meta.get(field, "MISSING")
            counts[value] = counts.get(value, 0) + 1

        print(field + ":", counts)

    return summaries


def inspect_script_risks():
    print("\n=== Original script risk checks ===")
    print("- df.isna() is called but df is not defined.")
    print("- patient_id is parsed from filename token split('_')[4] but this parsed variable is not used.")
    print("- sc.read_10x_h5 requires scanpy; h5py-only inspection is useful for lightweight validation.")
    print("- obs_names are not explicitly prefixed by GSM before concatenation; barcode collisions are possible.")
    print("- Dense np.isnan(adata_combined.X) may be memory-heavy or invalid for sparse matrices.")
    print("- Multiple metadata dimensions exist: N/T label, anatomical site, patient PT1-PT3, treatment status and stage.")
    print("- Metadata sanity checker should validate consistency between sample label prefixes and tumor_treatment.")


def cleanup_raw_files():
    print("\n=== Cleaning raw downloaded files ===")

    removed = 0

    for filename in os.listdir(RAW_DIR):
        filepath = os.path.join(RAW_DIR, filename)

        if os.path.isfile(filepath):
            os.remove(filepath)
            removed += 1

    print(f"Removed {removed} raw files from: {RAW_DIR}")


def run_case_017_inspection(keep_raw=False):
    print("Case:", PAPER_ID)
    print("GEO accession:", GEO_ACCESSION)
    print("Raw directory:", RAW_DIR)
    print("Keep raw files:", keep_raw)
    print("Mode: h5py-only per-sample 10x h5 inspection")

    ensure_directories()
    download_raw_tar()
    extract_raw_tar()

    files = list_gsm_files()
    h5_files, other_files = classify_files(files)
    inspect_metadata_vs_files(files)
    inspect_all_h5_files(h5_files)
    inspect_script_risks()

    print("\n=== Tool module implications ===")
    print("- src/tenx_h5_standardizer.py should support many per-sample 10x h5 files.")
    print("- src/metadata_checker.py should validate curated multi-site patient metadata.")
    print("- src/sample_id_parser.py should parse N/T, site and PT tokens from patient_id labels.")
    print("- src/gene_mapping_and_deduplication.py should run after concatenation.")
    print("- src/validation.py should use sparse-safe NaN checks.")

    if not keep_raw:
        cleanup_raw_files()

    print("\nInspection completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Inspect Case 017: 2022_Loret / GSE201047"
    )
    parser.add_argument(
        "--keep-raw",
        action="store_true",
        help="Keep downloaded raw files after inspection."
    )

    args = parser.parse_args()
    run_case_017_inspection(keep_raw=args.keep_raw)
