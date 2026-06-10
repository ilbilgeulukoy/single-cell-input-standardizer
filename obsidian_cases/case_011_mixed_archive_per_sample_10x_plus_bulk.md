# Case 011 - Mixed archive with per-sample 10x triplets, pooled single-cell triplets and bulk RNA-seq

## Case identity

- Case ID: case_011_2023_hippen_gse217517
- Source example: GSE217517
- Dataset context: high-grade serous ovarian cancer
- Input category: mixed GEO archive
- Selected modality: per-sample single-cell RNA-seq

## Why this case matters

This archive contains several different data types together.

It includes:

- 8 per-sample single-cell 10x Matrix Market triplets
- 2 pooled single-cell 10x Matrix Market triplets
- 24 bulk RNA-seq STAR tables

The tool must not treat every file in the GEO archive as the same data type.

## Detected files

Total GSM files detected:

54

Selected per-sample single-cell triplets:

- GSM6720925 single-cell triplet for patient 2251
- GSM6720926 single-cell triplet for patient 2267
- GSM6720927 single-cell triplet for patient 2283
- GSM6720928 single-cell triplet for patient 2293
- GSM6720929 single-cell triplet for patient 2380
- GSM6720930 single-cell triplet for patient 2428
- GSM6720931 single-cell triplet for patient 2467
- GSM6720932 single-cell triplet for patient 2497

Pooled single-cell triplets detected but excluded:

- GSM6720933 pooled single-cell 12162021
- GSM6720934 pooled single-cell 01132022

Bulk RNA-seq STAR tables detected but excluded:

- bulk_chunk_ribo
- bulk_dissociated_ribo
- bulk_dissociated_polyA

## First selected triplet

Sample:

GSM6720925

Files:

- GSM6720925_single_cell_matrix_2251.mtx.gz
- GSM6720925_single_cell_features_2251.tsv.gz
- GSM6720925_single_cell_barcodes_2251.tsv.gz

Matrix shape:

36,601 genes x 13,240 cells

Features shape:

36,601 rows x 3 columns

Barcodes shape:

13,240 rows x 1 column

The matrix is genes x cells before transpose.

AnnData construction requires transpose.

## Feature table structure

The feature table has 3 columns:

- gene_ensembl_id
- gene_symbol
- feature_type

The gene symbol column should be used as var_names.

Example rows:

- ENSG00000243485 / MIR1302-2HG / Gene Expression
- ENSG00000237613 / FAM138A / Gene Expression
- ENSG00000186092 / OR4F5 / Gene Expression

## Main bug detected in original script

The original script risks selecting the first matrix, features and barcode files for every GSM.

Risky logic:

mtx_files = next((f for f in mtx_paths if os.path.basename(f)), None)

This does not filter by GSM accession.

Correct logic:

Match matrix, features and barcodes by GSM prefix.

For example:

- GSM6720925 matrix must be paired with GSM6720925 features and GSM6720925 barcodes
- GSM6720926 matrix must be paired with GSM6720926 features and GSM6720926 barcodes

## Correct file matching rule

For each selected GSM:

1. Find matrix file that starts with the GSM accession and contains `_single_cell_matrix_`.
2. Find features file that starts with the GSM accession and contains `_single_cell_features_`.
3. Find barcode file that starts with the GSM accession and contains `_single_cell_barcodes_`.
4. Exclude pooled files.
5. Exclude bulk files.

## Metadata strategy

The script uses curated sample-level metadata.

Each GSM maps to a patient ID:

- GSM6720925 -> 2251
- GSM6720926 -> 2267
- GSM6720927 -> 2283
- GSM6720928 -> 2293
- GSM6720929 -> 2380
- GSM6720930 -> 2428
- GSM6720931 -> 2467
- GSM6720932 -> 2497

Constant metadata fields:

- dataset_id: GSE217517
- cancer_type: Ovarian
- tumor_site: Unknown
- metastasis_site: Unknown
- tumor_treatment: No
- cancer_site_origin: Unknown
- tumour_grade: HGSC
- tumour_stage: Unknown
- histological_subtype: Serous
- patient_ethnicity: Unknown

## Correct standardization recipe

1. Download and extract GSE217517_RAW.tar.
2. Detect all file types in the archive.
3. Keep only per-sample single-cell triplets for GSM6720925 to GSM6720932.
4. Exclude pooled single-cell triplets.
5. Exclude bulk RNA-seq STAR tables.
6. Match matrix, features and barcodes by GSM prefix.
7. Read each matrix.mtx.gz.
8. Transpose genes x cells to cells x genes.
9. Prefix cell barcodes with GSM accession.
10. Use gene symbols from the feature table as var_names.
11. Add curated sample-level metadata to obs.
12. Concatenate selected samples.
13. Apply gene mapping.
14. Remove invalid genes.
15. Deduplicate gene symbols using sparse aggregation.
16. Validate output.
17. Write standardized h5ad.
18. Delete raw files after processing.

## Reusable rule

A GEO archive can contain several modalities.

The tool must classify files before processing.

For this case:

selected modality = per-sample single-cell Matrix Market triplets

excluded modalities:

- pooled single-cell
- bulk RNA-seq STAR tables

## Product lesson

This case teaches the standardizer modality filtering.

It is not enough to detect 10x files.

The tool must also decide which 10x files belong to the intended analysis scope.
