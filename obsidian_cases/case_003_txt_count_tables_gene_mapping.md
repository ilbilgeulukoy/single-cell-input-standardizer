# Case 003 - TXT count tables with gene mapping and duplicate gene aggregation

## Case identity

- Case ID: case_003_2020_tongton_gse130000
- Source example: GSE130000
- Dataset context: ovarian cancer single-cell RNA-seq
- Technology: Drop-seq
- Input category: compressed TXT count tables
- Tool lesson: gene mapping and duplicate gene aggregation after AnnData construction

## Input pattern

The GEO supplementary archive contains 8 GSM-level `.txt.gz` count tables:

- GSM3729170_P1_dge.txt.gz
- GSM3729171_P2_dge.txt.gz
- GSM3729172_P3_dge.txt.gz
- GSM3729173_P4_dge.txt.gz
- GSM3729174_M1_dge.txt.gz
- GSM3729175_M2_dge.txt.gz
- GSM3729176_R1_dge.txt.gz
- GSM3729177_R2_dge.txt.gz

The first inspected file had:

- first column: `GENE`
- remaining columns: cell barcodes
- preview shape: 5 rows x 9926 columns

Example first genes:

- A1BG
- A1BG-AS1
- A2M
- A2M-AS1
- A2MP1

## Main problem

This case has two standardization layers.

First, the raw count table is oriented as:

genes x cells

AnnData expects:

cells x genes

Therefore, the object must be transposed.

Second, the gene identifiers need additional cleaning and harmonization.

The previous working script performed:

- gene mapping
- removal of genes flagged as invalid or unwanted
- duplicate gene symbol aggregation by summing counts

## Detection logic

A file belongs to this case type when:

- file extension is `.txt.gz`
- file starts with a GSM accession
- file is tab-delimited
- first column is named `GENE`
- remaining columns are numeric cell barcode count columns
- raw orientation is genes x cells

## AnnData recipe

1. Read each TXT count table using tab delimiter.
2. Treat the first column as gene symbols.
3. Treat remaining columns as cell barcodes.
4. Transpose to cells x genes.
5. Set obs_names and var_names as strings.
6. Add sample-level metadata.
7. Build one AnnData object per sample.
8. Align common genes.
9. Concatenate samples.
10. Validate missing values.

## Gene mapping recipe

After concatenation:

1. Create a gene feature table from `adata.var_names`.
2. Map dataset gene symbols to standardized gene identifiers.
3. Add gene metadata to `adata.var`.
4. Identify genes flagged for removal.
5. Remove invalid or unwanted genes.
6. Use the selected standardized gene symbol field for deduplication.
7. Aggregate duplicated gene symbols by summing counts.
8. Create a new AnnData object with unique gene symbols.

## Reusable rule

If gene symbols are inconsistent, duplicated, or contain non-standard entries, the tool should run a gene mapping and deduplication step before final h5ad export.

## Product lesson

This case supports a new module:

`gene_mapping_and_deduplication`

It also expands `count_table_standardizer` because TXT/TSV gene-by-cell tables should be treated as part of the same general count-table family as CSV count tables.
