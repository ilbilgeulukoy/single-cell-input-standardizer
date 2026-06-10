# Case 004 - Local 10x-like Matrix Market triplets with duplicated gene symbols

## Case identity

- Case ID: case_004_2021_hippen_gse158937
- Source example: GSE158937
- Dataset context: high-grade serous ovarian cancer single-cell RNA-seq
- Input category: per-sample 10x-like Matrix Market triplets
- Tool lesson: correct barcode assignment and sparse duplicated gene aggregation

## Input pattern

The GEO archive contains 3 samples and 9 files total.

Each sample has one complete triplet:

- matrix.mtx.gz
- features.tsv.gz
- barcodes.tsv.gz

Detected files:

- GSM4816045_matrix_16030X2_HJVMLDMXX.mtx.gz
- GSM4816045_features_16030X2_HJVMLDMXX.tsv.gz
- GSM4816045_barcodes_16030X2_HJVMLDMXX.tsv.gz
- GSM4816046_matrix_16030X3_HJTWLDMXX.mtx.gz
- GSM4816046_features_16030X3_HJTWLDMXX.tsv.gz
- GSM4816046_barcodes_16030X3_HJTWLDMXX.tsv.gz
- GSM4816047_matrix_16030X4_HJTWLDMXX.mtx.gz
- GSM4816047_features_16030X4_HJTWLDMXX.tsv.gz
- GSM4816047_barcodes_16030X4_HJTWLDMXX.tsv.gz

## First sample inspection

For GSM4816045:

- matrix shape: 36601 x 7123
- features shape: 36601 x 3
- barcodes shape: 7123 x 1

This means:

- rows of matrix match features
- columns of matrix match barcodes
- raw orientation is genes x cells
- AnnData requires cells x genes

Therefore, the matrix must be transposed before AnnData construction.

## Features structure

The features file has three columns:

1. Ensembl gene ID
2. gene symbol
3. feature type

Example rows:

- ENSG00000243485 / MIR1302-2HG / Gene Expression
- ENSG00000237613 / FAM138A / Gene Expression
- ENSG00000186092 / OR4F5 / Gene Expression

The previous working workflow used the second column as gene symbols.

## Main problem

This case has three important standardization issues:

1. Matrix orientation must be detected.
2. Cell barcodes must be read explicitly from barcodes.tsv.gz.
3. Duplicated gene symbols must be merged.

The first sample contained 10 duplicated gene symbols.

## Correct AnnData recipe

1. Read matrix.mtx.gz as sparse Matrix Market.
2. Read features.tsv.gz.
3. Read barcodes.tsv.gz.
4. Check matrix shape against features and barcodes.
5. Transpose matrix if rows match features and columns match barcodes.
6. Assign GSM-prefixed real barcodes to obs_names.
7. Assign gene symbols from features column 1 to var_names.
8. Add Ensembl IDs and feature types to var.
9. Merge duplicated gene symbols.
10. Add sample-level metadata to obs.
11. Concatenate all samples.
12. Validate missing values.
13. Write standardized h5ad.
14. Delete temporary raw files.

## Important correction from old workflow

The old script prefixed default obs_names after reading the matrix.

A safer rule is:

read barcodes.tsv.gz explicitly, then set obs_names from those barcodes.

This preserves the real cell barcode identity.

## Sparse duplicate gene aggregation

The old workflow converted X to dense using toarray before groupby.

This can work for small datasets but is unsafe for large single-cell matrices.

The reusable tool should merge duplicated genes using a sparse aggregation matrix instead.

## Reusable rule

If a Matrix Market triplet has:

- matrix rows equal to number of features
- matrix columns equal to number of barcodes

then the matrix is genes x cells and must be transposed for AnnData.

## Product lesson

This case supports:

- `tenx_mtx_standardizer`
- `gene_mapping_and_deduplication`
- `metadata_checker`

It is the locally processable version of the large Matrix Market archive pattern.
