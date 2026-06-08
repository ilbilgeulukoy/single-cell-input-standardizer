# Case 001: 2021_Olelekan / GSE147082

## Goal

This case demonstrates how to standardize GEO-downloaded single-cell RNA-seq count tables into a clean AnnData / h5ad object.

The dataset was provided as compressed CSV count tables rather than a standard 10x Matrix Market folder or an existing h5ad file.

## Input source

- Dataset: GSE147082
- Case name: 2021_Olelekan
- Source: GEO supplementary files
- Archive: GSE147082_RAW.tar
- Raw file policy: temporary download only, deleted after processing

## Detected raw files

The GEO archive contained six sample-level compressed CSV count tables:

- GSM4416534_PT-3232.csv.gz
- GSM4416535_PT-5150.csv.gz
- GSM4416536_PT-6885.csv.gz
- GSM4416537_PT-4806.csv.gz
- GSM4416538_PT-3401.csv.gz
- GSM4416539_PT-2834.csv.gz

## Input format interpretation

The first inspected file had the following structure:

- First column: gene symbols
- Remaining columns: cell barcodes
- Raw matrix orientation: genes x cells

Example:

| gene | cell_1 | cell_2 | cell_3 |
|---|---:|---:|---:|
| TSPAN6 | 0 | 0 | 1 |
| DPM1 | 9 | 2 | 0 |
| SCYL3 | 0 | 0 | 0 |

AnnData expects:

- rows = cells / observations
- columns = genes / variables

Therefore, the matrix must be transposed before AnnData creation.

Key conversion:

AnnData(df.T)

## Per-sample AnnData objects

| Sample ID | Patient ID | Cells | Genes before alignment |
|---|---|---:|---:|
| GSM4416534 | PT-3232 | 1102 | 19297 |
| GSM4416535 | PT-5150 | 1244 | 24776 |
| GSM4416536 | PT-6885 | 1071 | 24063 |
| GSM4416537 | PT-4806 | 1108 | 26499 |
| GSM4416538 | PT-3401 | 3451 | 32024 |
| GSM4416539 | PT-2834 | 1909 | 26221 |

## Metadata strategy

Sample-level metadata was manually curated from GEO and expanded to each cell in `adata.obs`.

Final obs columns:

- sample_id
- patient_id
- dataset_id
- cancer_type
- tumor_site
- metastasis_site
- tumor_treatment
- cancer_site_origin
- tumour_grade
- tumour_stage
- histological_subtype
- patient_ethnicity
- source_file
- paper_id
- geo_accession

## Gene alignment

The six samples had different gene sets.

To avoid missing values during concatenation, all samples were restricted to the intersection of shared genes.

- Genes in first sample: 19297
- Common genes across all samples: 16041
- Concatenation strategy: inner join on common genes

## Final AnnData object

Final standardized object:

- Cells: 9885
- Genes: 16041
- Samples: 6
- Format: h5ad

Final shape:

9885 cells x 16041 genes

## Validation

Validation checks:

| Check | Result |
|---|---|
| NaN in X | 0 |
| NaN in obs | False |
| NaN in var | False |
| gene_symbol in var | Present |
| raw files deleted | True |

## Standardization recipe

1. Download `GSE147082_RAW.tar` from GEO.
2. Extract the temporary archive.
3. Detect GSM-level `.csv.gz` count tables.
4. Inspect the first CSV file.
5. Identify the first column as gene symbols.
6. Identify remaining columns as cell barcodes.
7. Interpret raw matrix orientation as genes x cells.
8. Rename the first column as `gene`.
9. Set `gene` as the dataframe index.
10. Create AnnData using `AnnData(df.T)`.
11. Assign cell barcodes to `adata.obs_names`.
12. Assign gene symbols to `adata.var_names`.
13. Add `gene_symbol` to `adata.var`.
14. Add manually curated sample-level metadata to `adata.obs`.
15. Build one AnnData object per sample.
16. Validate X, obs, and var.
17. Align samples to common genes.
18. Concatenate samples into one AnnData object.
19. Save the standardized h5ad file.
20. Delete temporary raw downloaded files.

## Key lesson

This case shows why single-cell input standardization is necessary before preprocessing.

The dataset was not provided in a ready-to-use 10x or h5ad format. It required file inspection, matrix orientation detection, metadata reconstruction, gene alignment, validation, and controlled raw file cleanup.
