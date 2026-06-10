# Case 009 - Global UMI matrix with cellInfo metadata table

## Case identity

- Case ID: case_009_2022_zhang_gse165897
- Source example: GSE165897
- Dataset context: high-grade serous ovarian cancer single-cell RNA-seq
- Input category: global TSV UMI count matrix plus cell metadata table
- Tool lesson: join global count matrix to cell metadata using exact cell IDs

## Input pattern

This case has two global files:

- GSE165897_UMIcounts_HGSOC.tsv.gz
- GSE165897_cellInfo_HGSOC.tsv.gz

Unlike per-sample cases, there is not one count file per sample.

Instead, all cells are stored in one global matrix and metadata is stored in one global cellInfo table.

## Count matrix structure

The count matrix preview had:

- preview shape: 5 rows x 51787 columns
- first column name: `Unnamed: 0`
- first column values: gene symbols
- remaining columns: cell IDs
- number of cell columns: 51786

Example cell IDs:

- AAACCTGCAGGTTTCA-EOC372_pPer
- AAACCTGGTCCGAATT-EOC372_pPer
- AAAGATGCATCTGGTA-EOC372_pPer
- AAAGTAGTCGCTTAGA-EOC372_pPer
- AAATGCCAGGTGCACA-EOC372_pPer

The raw orientation is:

genes x cells

Therefore, the matrix must be transposed before AnnData creation.

## Cell metadata structure

The metadata table has a `cell` column.

Metadata columns detected:

- cell
- sample
- patient_id
- treatment_phase
- anatomical_location
- cell_type
- cell_subtype
- nCount_RNA
- nFeature_RNA
- percent.mt

## Join result

All matrix cell IDs matched metadata cell IDs.

- matrix cell columns: 51786
- matched metadata cells: 51786
- unmatched cells: 0

This is an ideal global matrix plus metadata join case.

## Sample parsing

The metadata `sample` column contains strings such as:

EOC372_primary_Peritoneum

The script parses sample_id by:

1. converting sample to string
2. stripping whitespace
3. splitting by underscore
4. keeping the first two tokens
5. joining them with underscore

Example:

EOC372_primary_Peritoneum -> EOC372_primary

The parsed sample_id is then mapped to GSM ID using the curated sample_to_gsm dictionary.

## Parsed sample distribution

- EOC1005_interval: 1288
- EOC1005_primary: 3071
- EOC136_interval: 3602
- EOC136_primary: 2501
- EOC153_interval: 2401
- EOC153_primary: 1380
- EOC227_interval: 2015
- EOC227_primary: 2178
- EOC349_interval: 1471
- EOC349_primary: 1483
- EOC372_interval: 4671
- EOC372_primary: 711
- EOC3_interval: 2503
- EOC3_primary: 3142
- EOC443_interval: 4463
- EOC443_primary: 2122
- EOC540_interval: 2043
- EOC540_primary: 1921
- EOC733_interval: 4077
- EOC733_primary: 1742
- EOC87_interval: 1491
- EOC87_primary: 1510

## Correct AnnData recipe

1. Read the global UMI count matrix.
2. Use the first column as gene symbols.
3. Use remaining columns as cell IDs.
4. Transpose the matrix to cells x genes.
5. Create AnnData.
6. Set obs_names from matrix cell IDs.
7. Set var_names from gene symbols.
8. Read the cellInfo metadata table.
9. Join metadata to obs by exact cell ID.
10. Parse sample_id from the metadata sample column.
11. Map sample_id to GSM ID using curated sample dictionary.
12. Merge curated sample-level metadata into obs.
13. Generate gene mapping table using dataset_gene_symbol.
14. Add standardized gene metadata to var.
15. Remove genes flagged as invalid.
16. Aggregate duplicated gene symbols using sparse aggregation.
17. Convert object obs columns to string before writing h5ad.
18. Validate X, obs and var.
19. Write standardized h5ad.
20. Delete temporary raw files.

## Main product lesson

Not all datasets arrive as one file per sample.

Some datasets arrive as:

one global matrix
+
one global metadata table

The standardizer needs a dedicated global-matrix workflow.

## Reusable modules needed

- `global_matrix_metadata_standardizer`
- `cell_metadata_joiner`
- `sample_id_parser`
- `gene_mapping_and_deduplication`

## Reusable rule

When a metadata table contains a `cell` column and the count matrix columns are cell IDs, AnnData obs must be joined by exact cell ID.

The tool should report the number of matched and unmatched cells before producing h5ad.
