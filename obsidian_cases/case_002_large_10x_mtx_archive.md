# Case 002 - Large GEO archive with 10x-like Matrix Market triplets

## Case identity

- Case ID: case_002_large_10x_mtx_archive
- Source example: GSE173682
- Dataset context: gynecologic malignancy single-cell RNA-seq / ATAC-related files
- Input category: large GEO archive with per-sample 10x-like triplets
- Tool lesson: data-size-aware inspection and 10x Matrix Market standardization

## Input pattern

The GEO archive is very large, approximately 15.8 GB.

Expected per-GSM RNA files:

- matrix.mtx.gz
- features.tsv.gz
- barcodes.tsv.gz

Additional ATAC fragment files may also be present:

- ATAC_fragments.tsv.gz

## Main problem

The archive is too large for default local personal-computer workflows.

A local download attempt may fail or consume unnecessary storage.

Therefore, the tool should support:

- manifest-based inspection
- filename-pattern-based detection
- server/HPC recommendation
- raw data cleanup policy

## Detection logic

A dataset can be classified as 10x-like Matrix Market when each sample has:

- one matrix file
- one features or genes file
- one barcodes file

Expected pattern:

- `*_matrix-*.mtx.gz`
- `*_features-*.tsv.gz`
- `*_barcodes-*.tsv.gz`

## Standardization recipe

For each sample:

1. Locate matrix, features, and barcodes files.
2. Read matrix.mtx.gz as a sparse Matrix Market object.
3. Transpose matrix to cells x genes if needed.
4. Read barcodes.tsv.gz and assign cell IDs.
5. Prefix cell barcodes with sample ID.
6. Read features.tsv.gz and assign gene symbols.
7. Detect duplicated gene symbols.
8. Merge duplicated genes by summing counts.
9. Add sample-level metadata to obs.
10. Align common genes across samples.
11. Concatenate samples.
12. Validate X, obs, and var.
13. Write standardized h5ad.

## Reusable rule

If files include matrix.mtx.gz, features.tsv.gz, and barcodes.tsv.gz, use a 10x-like Matrix Market standardizer.

If archive size is large, do not download by default on local machines.

## Product lesson

This case supports two modules:

- `tenx_mtx_standardizer`
- `data_size_policy`

The tool should warn users when a public archive is too large and recommend server/HPC mode.
