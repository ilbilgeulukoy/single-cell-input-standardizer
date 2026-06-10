# Single-Cell Input Standardizer

A Python-based bioinformatics tool for inspecting messy single-cell input files and generating reproducible AnnData / h5ad standardization workflows.

This project focuses on the data ingestion step before preprocessing, quality control, clustering, annotation, or downstream biological analysis.

## What problem does this tool solve?

Single-cell datasets downloaded from GEO, supplementary files, publications, or lab servers are often not directly ready for Scanpy or Seurat workflows.

Common issues include:

- unknown file formats
- custom CSV or TSV count matrices
- missing or incomplete metadata
- unclear matrix orientation
- genes x cells matrices that need transposition
- cell barcodes stored as column names
- gene symbols stored in the first column
- sample IDs hidden in filenames
- patient metadata requiring manual curation
- inconsistent gene sets across samples

This tool aims to help users move from messy downloaded files to a standardized AnnData / h5ad object.

## Core idea

The tool inspects input files, detects their structure, identifies risks, and recommends or executes an AnnData conversion strategy.

Typical workflow:

1. Inspect downloaded single-cell files.
2. Detect the input format.
3. Identify whether the data is 10x Matrix Market, 10x HDF5, existing h5ad, custom count table, or metadata table.
4. Detect matrix orientation when possible.
5. Identify gene IDs, cell barcodes, sample IDs, and metadata fields.
6. Build or recommend a standardized AnnData object.
7. Validate X, obs, and var.
8. Save a clean h5ad output.
9. Avoid storing large raw files permanently.

## Current supported workflow

The first implemented workflow supports custom compressed CSV count tables where:

- the first column contains gene symbols
- the remaining columns contain cell barcodes
- the raw matrix orientation is genes x cells
- AnnData must be created using AnnData(df.T)

This workflow can:

- download temporary example data
- inspect compressed CSV count tables
- create one AnnData object per sample
- add sample-level metadata to adata.obs
- store gene symbols in adata.var
- align samples using common genes
- concatenate samples into a final AnnData object
- validate missing values
- write a standardized h5ad file
- delete raw downloaded files after processing

## Example validation case

The current implementation was validated using a real GEO-derived single-cell RNA-seq example containing multiple compressed CSV count tables.

The example demonstrates a common real-world ingestion problem:

- non-standard single-cell file structure
- gene-by-cell count matrices
- metadata reconstructed from sample information
- sample-level AnnData creation
- common gene alignment
- final h5ad generation

The validation case is documented in:

- validation case YAML file
- markdown report
- reproducible build script

The raw data is not stored in the repository.

## Repository structure

single-cell-input-standardizer/
  cases/
    case_001_2021_olelekan_gse147082.yaml

  reports/
    case_001_2021_olelekan_report.md

  scripts/
    case_001_build_adata_2021_olelekan.py
    inspect_h5ad.py

  src/
    __init__.py

  README.md
  requirements.txt
  .gitignore

## Raw data policy

Raw downloaded files are treated as temporary files.

The tool can download example data, process it, generate a standardized output, and then remove the raw archive and extracted count matrices from local storage.

This keeps the project lightweight and avoids storing large omics files in GitHub.

Ignored file types include:

- raw archives
- compressed count matrices
- h5ad outputs
- h5 files
- Matrix Market files
- generated local outputs

## Installation

Create and activate a Python environment.

Example with conda:

conda create -n sc-standardizer python=3.11 -y
conda activate sc-standardizer

Install dependencies:

pip install -r requirements.txt

## Running the current example workflow

Run the current standardization workflow:

python scripts/case_001_build_adata_2021_olelekan.py

By default, raw downloaded files are deleted after processing.

For debugging only, raw files can be kept temporarily:

python scripts/case_001_build_adata_2021_olelekan.py --keep-raw

Inspect a generated h5ad object:

python scripts/inspect_h5ad.py path/to/file.h5ad

## Planned general tool interface

The project will be extended toward a more general interface where users can provide:

- a folder path
- a file manifest
- a count matrix
- a metadata table
- a 10x Matrix Market folder
- a 10x HDF5 file
- an existing h5ad file

The tool will return:

- detected format
- detected matrix orientation
- required transformation
- metadata warnings
- AnnData construction strategy
- validation summary
- recommended h5ad output path

