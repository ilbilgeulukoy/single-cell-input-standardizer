# Case 008 - CSV cell-by-gene count tables with CellId column

## Case identity

- Case ID: case_008_2018_shih_gse118828
- Source example: GSE118828
- Dataset context: serous epithelial ovarian cancer single-cell RNA-seq
- Input category: compressed CSV count tables
- Tool lesson: detect cell-by-gene CSV layouts using a CellId column

## Input pattern

The GEO archive contained 18 `.csv.gz` count tables.

The selected working subset contained 16 files with curated metadata.

Two files were skipped because they did not have curated metadata in the working script.

Skipped files:

- GSM3348304_565_Cystadenoma_S1.counts.umiCounts.aboveBackground.table.csv.gz
- GSM3348307_TB10040568_NORMAL_S1.counts.umiCounts.table.csv.gz

## First selected file inspection

For GSM3348303_553_Perit_S1.counts.umiCounts.aboveBackground.table.csv.gz:

- preview shape: 5 rows x 26365 columns
- first column name: `CellId`
- first column values: cell identifiers
- remaining columns: gene symbols

Example CellId values:

- aacagctgaattagcacg
- aagccactaggtaacagc
- aagccatcgcctgagctt
- aagtatcagacttccaag
- acaaggtgagacagatgt

## Main problem

This case is different from earlier CSV count table cases.

Earlier CSV cases had:

genes x cells

This case has:

cells x genes

Therefore, the matrix must not be transposed.

## Detection logic

A CSV count table should be classified as cell-by-gene when:

- the first column is named `CellId`, `cell_id`, `cell`, `barcode` or similar
- rows correspond to cells
- remaining columns correspond to gene symbols
- numeric expression values occupy the remaining columns

## Correct AnnData recipe

1. Read the `.csv.gz` table with pandas.
2. Detect the `CellId` column.
3. Set `obs_names` from `CellId`.
4. Drop `CellId` from the expression matrix.
5. Create AnnData directly with `AnnData(X)`.
6. Do not transpose.
7. Assign `var_names` from the remaining gene columns.
8. Add sample-level metadata to obs.
9. Align common genes across samples.
10. Concatenate samples.
11. Generate gene mapping table using `dataset_gene_symbol`.
12. Add standardized gene metadata to var.
13. Remove genes flagged as invalid.
14. Aggregate duplicated gene symbols after mapping.
15. Validate X, obs and var.
16. Write standardized h5ad.
17. Delete temporary raw files.

## Metadata lesson

The archive contains 18 CSV files, but only 16 were selected because the metadata dictionary contained 16 GSM entries.

The tool should report skipped GSM files clearly.

## Reusable rule

If a CSV count table contains a `CellId` column, the standardizer should treat the table as cells x genes and avoid transposition.

## Product lesson

This case extends:

- `count_table_standardizer`
- `orientation_detector`
- `metadata_checker`
- `gene_mapping_and_deduplication`

It teaches the assistant not to blindly transpose all CSV count tables.
