import os
import argparse
from urllib.request import urlretrieve

import pandas as pd


PAPER_ID = "2022_Zhang"
GEO_ACCESSION = "GSE165897"

BASE_DIR = os.path.join("data", "case_009_2022_zhang")
RAW_DIR = os.path.join(BASE_DIR, "raw")
OUTPUT_DIR = os.path.join(BASE_DIR, "outputs")

COUNT_FILENAME = "GSE165897_UMIcounts_HGSOC.tsv.gz"
META_FILENAME = "GSE165897_cellInfo_HGSOC.tsv.gz"

COUNT_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165897/suppl/GSE165897%5FUMIcounts%5FHGSOC%2Etsv%2Egz"
META_URL = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE165nnn/GSE165897/suppl/GSE165897%5FcellInfo%5FHGSOC%2Etsv%2Egz"

COUNT_PATH = os.path.join(RAW_DIR, COUNT_FILENAME)
META_PATH = os.path.join(RAW_DIR, META_FILENAME)


SAMPLE_METADATA = {
    "GSM5057576": {"sample_label": "EOC1005_primary", "treatment": "No", "origin": "Peritoneum", "stage": "IVA"},
    "GSM5057577": {"sample_label": "EOC1005_interval", "treatment": "Yes", "origin": "Omentum", "stage": "IVA"},
    "GSM5057578": {"sample_label": "EOC136_primary", "treatment": "No", "origin": "Mesentery", "stage": "IVA"},
    "GSM5057579": {"sample_label": "EOC136_interval", "treatment": "Yes", "origin": "Omentum", "stage": "IVA"},
    "GSM5057580": {"sample_label": "EOC153_primary", "treatment": "No", "origin": "Omentum", "stage": "IVA"},
    "GSM5057581": {"sample_label": "EOC153_interval", "treatment": "Yes", "origin": "Omentum", "stage": "IVA"},
    "GSM5057582": {"sample_label": "EOC227_primary", "treatment": "No", "origin": "Omentum", "stage": "IVA"},
    "GSM5057583": {"sample_label": "EOC227_interval", "treatment": "Yes", "origin": "Omentum", "stage": "IVA"},
    "GSM5057584": {"sample_label": "EOC349_primary", "treatment": "No", "origin": "Omentum", "stage": "IVB"},
    "GSM5057585": {"sample_label": "EOC349_interval", "treatment": "Yes", "origin": "Omentum", "stage": "IVB"},
    "GSM5057586": {"sample_label": "EOC372_primary", "treatment": "No", "origin": "Peritoneum", "stage": "IIIC"},
    "GSM5057587": {"sample_label": "EOC372_interval", "treatment": "Yes", "origin": "Peritoneum", "stage": "IIIC"},
    "GSM5057588": {"sample_label": "EOC3_primary", "treatment": "No", "origin": "Peritoneum", "stage": "IVA"},
    "GSM5057589": {"sample_label": "EOC3_interval", "treatment": "Yes", "origin": "Omentum", "stage": "IVA"},
    "GSM5057590": {"sample_label": "EOC443_primary", "treatment": "No", "origin": "Omentum", "stage": "IVA"},
    "GSM5057591": {"sample_label": "EOC443_interval", "treatment": "Yes", "origin": "Omentum", "stage": "IVA"},
    "GSM5057592": {"sample_label": "EOC540_primary", "treatment": "No", "origin": "Omentum", "stage": "IIIC"},
    "GSM5057593": {"sample_label": "EOC540_interval", "treatment": "Yes", "origin": "Omentum", "stage": "IIIC"},
    "GSM5057594": {"sample_label": "EOC733_primary", "treatment": "No", "origin": "Peritoneum", "stage": "IVA"},
    "GSM5057595": {"sample_label": "EOC733_interval", "treatment": "Yes", "origin": "Omentum", "stage": "IVA"},
    "GSM5057596": {"sample_label": "EOC87_primary", "treatment": "No", "origin": "Peritoneum", "stage": "IIIC"},
    "GSM5057597": {"sample_label": "EOC87_interval", "treatment": "Yes", "origin": "Omentum", "stage": "IIIC"},
}


def ensure_directories() -> None:
    os.makedirs(RAW_DIR, exist_ok=True)
    os.makedirs(OUTPUT_DIR, exist_ok=True)


def download_file(url: str, path: str) -> None:
    if os.path.exists(path):
        print(f"File already exists: {path}")
        return

    print("Downloading:", url)
    print("Output:", path)

    urlretrieve(url, path)

    print("Download completed.")


def download_inputs() -> None:
    download_file(COUNT_URL, COUNT_PATH)
    download_file(META_URL, META_PATH)


def inspect_count_matrix() -> None:
    print("\n=== Inspecting UMI count matrix ===")
    print("File:", COUNT_PATH)

    df_preview = pd.read_csv(COUNT_PATH, sep="\t", compression="gzip", nrows=5)

    print("\nPreview shape:")
    print(df_preview.shape)

    print("\nFirst columns:")
    print(list(df_preview.columns[:10]))

    print("\nPreview rows:")
    print(df_preview.head())

    first_col = df_preview.columns[0]

    print("\nCount matrix interpretation:")
    print("First column:", first_col)
    print("First column is treated as gene symbols.")
    print("Remaining columns are cell IDs.")
    print("Raw orientation: genes x cells")
    print("Required transpose for AnnData: True")


def inspect_metadata() -> pd.DataFrame:
    print("\n=== Inspecting cell metadata ===")
    print("File:", META_PATH)

    df_meta_preview = pd.read_csv(META_PATH, sep="\t", compression="gzip", nrows=10)

    print("\nMetadata preview shape:")
    print(df_meta_preview.shape)

    print("\nMetadata columns:")
    print(list(df_meta_preview.columns))

    print("\nMetadata preview rows:")
    print(df_meta_preview.head())

    if "cell" in df_meta_preview.columns:
        print("\nMetadata contains cell column: True")
    else:
        print("\nMetadata contains cell column: False")

    if "sample" in df_meta_preview.columns:
        sample_ids = (
            df_meta_preview["sample"]
            .astype(str)
            .str.strip()
            .str.split("_")
            .str[:2]
            .str.join("_")
        )
        print("\nParsed sample IDs from metadata sample column:")
        print(sample_ids.value_counts())

    return df_meta_preview


def inspect_header_and_join_logic() -> None:
    print("\n=== Inspecting full count matrix header and metadata join logic ===")

    count_header = pd.read_csv(COUNT_PATH, sep="\t", compression="gzip", nrows=0).columns.astype(str)
    cell_columns = [c for c in count_header if c != count_header[0]]

    print("Number of cell columns in count matrix:", len(cell_columns))
    print("First count matrix cell IDs:")
    for cell in cell_columns[:10]:
        print("-", cell)

    df_meta = pd.read_csv(META_PATH, sep="\t", compression="gzip", usecols=["cell", "sample"])
    metadata_cells = set(df_meta["cell"].astype(str))

    matched = sum(1 for c in cell_columns if c in metadata_cells)
    print("\nCell IDs matched in metadata:", matched)
    print("Cell IDs not matched in metadata:", len(cell_columns) - matched)

    parsed_sample_ids = (
        df_meta["sample"]
        .astype(str)
        .str.strip()
        .str.split("_")
        .str[:2]
        .str.join("_")
    )

    print("\nParsed sample IDs across full metadata:")
    print(parsed_sample_ids.value_counts().sort_index())

    sample_to_gsm = {v["sample_label"]: k for k, v in SAMPLE_METADATA.items()}
    observed_sample_ids = set(parsed_sample_ids)
    missing_from_dict = sorted(observed_sample_ids - set(sample_to_gsm))

    print("\nParsed sample IDs missing from sample_to_gsm dictionary:")
    print(missing_from_dict)

    print("\nJoin rule:")
    print("- count matrix columns are cell IDs")
    print("- metadata has a cell column")
    print("- AnnData obs should be joined by exact cell ID")
    print("- sample_id should be parsed from metadata sample column")
    print("- sample_id should be mapped to GSM via sample_to_gsm")


def cleanup_raw_files() -> None:
    print("\n=== Cleaning raw downloaded files ===")

    removed = 0

    for filename in os.listdir(RAW_DIR):
        filepath = os.path.join(RAW_DIR, filename)

        if os.path.isfile(filepath):
            os.remove(filepath)
            removed += 1

    print(f"Removed {removed} raw files from: {RAW_DIR}")


def run_case_009_inspection(keep_raw: bool = False) -> None:
    print("Case:", PAPER_ID)
    print("GEO accession:", GEO_ACCESSION)
    print("Raw directory:", RAW_DIR)
    print("Keep raw files:", keep_raw)

    ensure_directories()
    download_inputs()

    inspect_count_matrix()
    inspect_metadata()
    inspect_header_and_join_logic()

    print("\n=== Curated metadata entries ===")
    print("Number of GSM/sample mappings:", len(SAMPLE_METADATA))
    for gsm, meta in SAMPLE_METADATA.items():
        print(f"- {gsm}: {meta['sample_label']} / treatment={meta['treatment']} / origin={meta['origin']} / stage={meta['stage']}")

    print("\n=== Tool module implications ===")
    print("- src/global_matrix_metadata_standardizer.py should support one global UMI matrix plus one cell metadata table.")
    print("- src/cell_metadata_joiner.py should join obs by exact cell IDs.")
    print("- src/sample_id_parser.py should parse sample_id from metadata sample strings.")
    print("- src/gene_mapping_and_deduplication.py should use dataset_gene_symbol mapping for var_names.")

    if not keep_raw:
        cleanup_raw_files()

    print("\nInspection completed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Inspect Case 009: 2022_Zhang / GSE165897"
    )
    parser.add_argument(
        "--keep-raw",
        action="store_true",
        help="Keep downloaded raw files after inspection."
    )

    args = parser.parse_args()
    run_case_009_inspection(keep_raw=args.keep_raw)
