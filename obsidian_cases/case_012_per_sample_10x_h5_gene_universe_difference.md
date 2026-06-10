# Case 012 - Per-sample 10x h5 files with gene universe differences

## Case identity

- Case ID: case_012_2024_brand_gse233615
- Source example: GSE233615
- Paper label: 2024_Brand
- Input category: per-sample 10x Genomics filtered_feature_bc_matrix.h5
- Selected modality: single-cell RNA-seq

## Why this case matters

This case introduces per-sample 10x h5 files.

Unlike Matrix Market triplets, each sample is stored as a single `.h5` file.

The tool must support 10x h5 detection and h5py-only inspection when scanpy is not available.

## Detected files

The archive contains 9 GSM files.

All 9 files are 10x h5 files:

- GSM7431434_E1_filtered_feature_bc_matrix.h5
- GSM7431435_E2_filtered_feature_bc_matrix.h5
- GSM7431436_E3_filtered_feature_bc_matrix.h5
- GSM7431437_E4_filtered_feature_bc_matrix.h5
- GSM7431438_E7_filtered_feature_bc_matrix.h5
- GSM7431439_E9_filtered_feature_bc_matrix.h5
- GSM7431440_E10_filtered_feature_bc_matrix.h5
- GSM7431441_E11_filtered_feature_bc_matrix.h5
- GSM7431442_E12_filtered_feature_bc_matrix.h5

No other files were detected.

## H5 structure

Each file has a top-level `matrix` group.

The matrix group contains:

- barcodes
- data
- features
- indices
- indptr
- shape

The features group contains:

- _all_tag_keys
- feature_type
- genome
- id
- name

This is standard 10x Genomics h5 structure.

## Orientation

The h5 matrix shape is genes x cells.

After `scanpy.read_10x_h5`, the AnnData shape should be cells x genes.

## Sample summary

- GSM7431434 / E1: 3,808 cells x 36,601 genes
- GSM7431435 / E2: 4,726 cells x 36,601 genes
- GSM7431436 / E3: 4,408 cells x 36,601 genes
- GSM7431437 / E4: 4,470 cells x 36,601 genes
- GSM7431438 / E7: 13,062 cells x 36,601 genes
- GSM7431439 / E9: 5,420 cells x 36,601 genes
- GSM7431440 / E10: 7,406 cells x 36,601 genes
- GSM7431441 / E11: 8,936 cells x 36,601 genes
- GSM7431442 / E12: 13,179 cells x 33,538 genes

Total cells:

65,415

Unique gene counts:

- 33,538
- 36,601

## Main standardization issue

One sample has a smaller gene universe than the others.

Most samples have 36,601 genes.

GSM7431442 / E12 has 33,538 genes.

The tool must not assume that all 10x h5 files in the same GEO archive have identical feature sets.

## Metadata validation

The filename contains both GSM accession and patient token.

Example:

GSM7431434_E1_filtered_feature_bc_matrix.h5

- GSM accession: GSM7431434
- patient token: E1

All filename patient tokens matched the curated metadata dictionary.

## Metadata strategy

The script uses curated sample-level metadata.

Sample-specific fields include:

- patient_id
- tumor_site
- metastasis_site
- tumour_grade
- patient_ethnicity

Constant fields include:

- dataset_id: GSE233615
- cancer_type: Ovarian
- tumor_treatment: No
- cancer_site_origin: Fallopian Tube
- tumour_stage: Unknown
- histological_subtype: Serous

This case includes both LGSC and HGSC samples, so tumour_grade must remain sample-specific.

## Duplicate gene names

Duplicated feature names were detected in every sample.

Most samples had 10 duplicated feature names.

Sample E12 had 24 duplicated feature names.

Gene deduplication is required after gene mapping.

## Correct standardization recipe

1. Download and extract GSE233615_RAW.tar.
2. Detect per-sample 10x h5 files.
3. Validate all GSM files against curated metadata.
4. Validate filename patient tokens against curated patient IDs.
5. Read each file with scanpy.read_10x_h5 during final processing.
6. Use h5py-only inspection when scanpy is unavailable.
7. Add GSM prefix to cell barcodes if needed.
8. Add curated sample-level metadata to obs.
9. Validate gene universe sizes across samples.
10. Apply explicit gene harmonization policy.
11. Concatenate samples.
12. Generate gene mapping table using dataset_gene_symbol.
13. Remove genes flagged as to_remove.
14. Deduplicate gene symbols using sparse aggregation.
15. Validate X, obs, var and sample counts.
16. Write standardized h5ad.
17. Delete raw files after processing.

## Product lesson

This case teaches the standardizer to support 10x h5 input files.

It also teaches that samples from the same dataset can have different gene universes.

A robust tool must detect and report this before concatenation.
