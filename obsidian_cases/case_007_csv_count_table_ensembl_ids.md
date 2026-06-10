# Case 007 - CSV count tables with Ensembl feature IDs

## Case identity

- Case ID: case_007_2023_guo_gse181955
- Source example: GSE181955
- Dataset context: high-grade serous ovarian cancer single-cell RNA-seq
- Input category: compressed CSV count tables
- Tool lesson: detect whether the first column contains gene symbols or Ensembl gene IDs

## Input pattern

The GEO archive contained 8 `.csv.gz` count tables.

Detected files:

- GSM5514787_N1_normal_ovary.matrix.csv.gz
- GSM5514788_OMT-1_CD45_POS.matrix.csv.gz
- GSM5514789_OMT-1_CD45_NEG.matrix.csv.gz
- GSM5514790_OMT-3_CD45_POS.matrix.csv.gz
- GSM5514791_OMT-3_CD45_NEG.matrix.csv.gz
- GSM5514792_T1.matrix.csv.gz
- GSM5514793_T6.matrix.csv.gz
- GSM5514794_T6_CD45_NEG.matrix.csv.gz

Only 6 files had curated metadata in the working script and were selected.

Skipped files:

- GSM5514787_N1_normal_ovary.matrix.csv.gz
- GSM5514793_T6.matrix.csv.gz

## First selected file inspection

For GSM5514788_OMT-1_CD45_POS.matrix.csv.gz:

- preview shape: 5 rows x 8951 columns
- first column name: `Unnamed: 0`
- first column values: Ensembl-like IDs starting with `ENSG`
- remaining columns: cell barcodes

Example first feature IDs:

- ENSG00000243485
- ENSG00000237613
- ENSG00000186092
- ENSG00000238009
- ENSG00000239945

## Main problem

The file is a gene/feature-by-cell count table and must be transposed for AnnData.

However, unlike earlier CSV count table cases, the first column does not contain gene symbols.

It contains Ensembl gene IDs.

This changes the gene mapping strategy.

## Detection logic

The first column should be classified as Ensembl-like when most values start with:

ENSG

In that case, the mapping key should be:

dataset_gene_ensembl_id

not:

dataset_gene_symbol

## Correct AnnData recipe

1. Read each selected `.csv.gz` count table with pandas.
2. Rename the first column to `feature_id`.
3. Set `feature_id` as the dataframe index.
4. Detect feature identifier type.
5. Transpose to cells x features.
6. Create AnnData with `AnnData(df.T)`.
7. Assign obs_names from cell barcode columns.
8. Assign var_names from Ensembl IDs.
9. Add sample-level metadata to obs.
10. Align common features.
11. Concatenate selected samples.
12. Generate gene mapping table using `dataset_gene_ensembl_id`.
13. Add standardized gene metadata to var.
14. Remove genes flagged as invalid.
15. Aggregate duplicated gene symbols after mapping.
16. Validate X, obs and var.
17. Write standardized h5ad.
18. Delete temporary raw files.

## Metadata lesson

This case also has a metadata selection issue.

The archive contains 8 CSV files, but only 6 were used because the working metadata dictionary contained 6 GSM entries.

The tool should report skipped GSM files clearly.

## Reusable rule

A CSV count table first column may contain either:

- gene symbols
- Ensembl gene IDs
- ambiguous feature identifiers

The standardizer should infer the feature identifier type before choosing the gene mapping strategy.

## Product lesson

This case supports a new module:

`feature_identifier_type_detector`

It also extends:

- `count_table_standardizer`
- `gene_mapping_and_deduplication`
- `metadata_checker`
