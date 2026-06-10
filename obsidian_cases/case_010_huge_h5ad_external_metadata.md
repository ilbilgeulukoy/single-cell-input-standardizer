# Case 010 - Huge h5ad-like matrix with external rich metadata

## Case identity

- Case ID: case_010_2022_ignacio_gse180661
- Source example: GSE180661
- Dataset context: high-grade serous ovarian cancer single-cell RNA-seq
- Input category: huge existing h5ad-like matrix plus external rich cell metadata
- Tool lesson: external metadata can define the authoritative cell subset

## Why this case matters

This is a large metadata-rich cohort dataset.

It is not a simple per-sample count table.

The input consists of:

- a very large matrix file: GSE180661_matrix.h5
- an external rich metadata table: GSE180661_GEO_cells.tsv.gz

The matrix file is too large for laptop processing and should be processed on server or HPC.

## Server location

The matrix file is already present on Derrick:

/data/ulukoy/ovarian_cancer/data/2022_Ignacio/a_data_raw/downloaded/GSE180661_matrix.h5

The metadata file is:

/data/ulukoy/ovarian_cancer/data/2022_Ignacio/a_data_raw/downloaded/GSE180661_GEO_cells.tsv.gz

## Matrix inspection

The HDF5 file has AnnData-like top-level keys:

- X
- obs
- uns
- var

It can be read with AnnData backed mode.

Observed matrix shape:

- 1,376,121 cells
- 32,223 genes

The matrix has no existing obs columns.

The var table includes:

- gene_ids
- feature_types

First var names include:

- MIR1302-2HG
- OR4F5
- AL627309.1
- AL627309.3
- AL627309.4

## Metadata inspection

The metadata table has 929,690 rows.

Metadata columns include:

- cell_id
- sample
- cell_type
- percent.mt
- nCount_RNA
- nFeature_RNA
- umap50_1
- umap50_2
- cluster_label
- cluster_label_sub
- cell_type_super
- patient_id
- tumor_subsite
- tumor_site
- tumor_supersite
- sort_parameters
- therapy
- surgery

## Cell matching result

Metadata cells found in matrix:

929,690 / 929,690

Metadata cells missing from matrix:

0

Matrix cells without metadata:

446,431

This means the matrix contains more cells than the curated external metadata table.

## Main decision

The external metadata should define the authoritative cell subset.

The standardized output should keep only cells present in GSE180661_GEO_cells.tsv.gz.

Matrix-only cells without metadata should be dropped unless the user explicitly requests a raw matrix universe output.

## Critical bug prevention

The metadata table must be indexed by `cell_id` before subsetting AnnData.

Wrong pattern:

Use numeric metadata index to subset AnnData.

Correct pattern:

1. Set metadata index to cell_id.
2. Keep only AnnData cells that exist in metadata.
3. Reorder metadata with df_meta.loc[adata.obs_names].
4. Attach metadata to adata.obs.

Correct logic:

df_meta = df_meta.set_index("cell_id")
common_cells = adata.obs_names.intersection(df_meta.index)
adata = adata[common_cells, :].copy()
adata.obs = df_meta.loc[adata.obs_names].copy()

## Clinical metadata harmonization

The source column is:

tumor_supersite

Derived tumor_site rule:

- Adnexa -> Primary
- Ascites -> Ascites
- everything else -> Metastasis

Derived metastasis_site rule:

- if tumor_site is Metastasis, keep tumor_supersite
- otherwise set empty string

Derived cancer_site_origin rule:

- if tumor_site is Primary, keep tumor_supersite
- otherwise set empty string

Constant standardized fields:

- dataset_id: GSE180661
- cancer_type: Ovarian
- tumor_treatment: No
- tumour_grade: HGSC
- tumour_stage: Unknown
- histological_subtype: Serous
- patient_ethnicity: Unknown

## Correct standardization recipe

1. Use the existing server-side matrix.h5 file.
2. Read matrix.h5 as AnnData-like h5ad object.
3. Read GSE180661_GEO_cells.tsv.gz metadata.
4. Set metadata index to cell_id.
5. Validate cell matching.
6. Keep only cells present in metadata.
7. Drop matrix-only cells without metadata.
8. Reorder metadata to match AnnData obs_names.
9. Attach metadata to obs.
10. Harmonize clinical metadata.
11. Select useful metadata columns or preserve rich metadata depending on output mode.
12. Use gene symbols as mapping key.
13. Generate gene mapping table.
14. Remove invalid genes.
15. Aggregate duplicated gene symbols.
16. Validate output.
17. Write standardized h5ad on server or HPC.

## Reusable rule

When an h5ad-like matrix contains more cells than an external metadata table, the tool should not assume this is an error.

It should report:

- matrix cells
- metadata cells
- matched cells
- metadata cells missing from matrix
- matrix cells missing from metadata

Then it should require a cell subset policy.

For this case:

metadata_authoritative_subset

## Product lesson

This case expands the standardizer from format detection into cohort curation.

The tool must understand that rich metadata can define the usable biological dataset.
