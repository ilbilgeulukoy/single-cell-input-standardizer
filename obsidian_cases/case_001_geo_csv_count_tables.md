# Case 001 - GEO compressed CSV count tables

## Case identity

- Case ID: case_001_geo_csv_count_tables
- Source example: GSE147082
- Dataset context: ovarian cancer single-cell RNA-seq
- Input category: custom compressed CSV count tables
- Tool lesson: gene-by-cell table standardization

## Input pattern

The downloaded GEO archive contained multiple sample-level `.csv.gz` count tables.

Each file corresponded to one GSM sample.

Observed structure:

- first column: gene symbols
- remaining columns: cell barcodes
- raw orientation: genes x cells

## Main problem

AnnData expects:

- rows = cells / observations
- columns = genes / variables

The raw file was:

- rows = genes
- columns = cells

Therefore, the matrix had to be transposed before AnnData creation.

## Detection logic

A file can be classified as a gene-by-cell count table when:

- the file is `.csv`, `.csv.gz`, `.tsv`, or `.tsv.gz`
- the first column contains gene-like strings
- most remaining columns are numeric
- remaining column names look like cell barcodes or cell IDs

## Standardization recipe

1. Read the count table with pandas.
2. Rename the first column as `gene`.
3. Set `gene` as the dataframe index.
4. Convert index and columns to string.
5. Create AnnData using `AnnData(df.T)`.
6. Set `adata.obs_names` from cell barcode columns.
7. Set `adata.var_names` from gene index.
8. Add `gene_symbol` to `adata.var`.
9. Add sample-level metadata to `adata.obs`.
10. Repeat for each sample.
11. Align all samples to common genes.
12. Concatenate AnnData objects.
13. Validate X, obs, and var.
14. Write standardized h5ad.
15. Delete temporary raw files.

## Reusable rule

If a count table is detected as `genes x cells`, the recommended AnnData constructor is:

AnnData(df.T)

## Validation result

Final standardized object:

- cells: 9885
- genes: 16041
- samples: 6
- NaN in X: 0
- NaN in obs: False
- NaN in var: False

## Product lesson

This case supports the `count_table_standardizer` module.

The tool should be able to detect custom count tables, infer orientation, recommend transposition, and produce an AnnData recipe.
