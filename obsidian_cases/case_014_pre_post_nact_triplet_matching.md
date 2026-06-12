# Case 014 - Pre/Post NACT per-sample Matrix Market triplets

## Case identity

- Case ID: case_014_2022_shen_gse191301
- Paper label: 2022_Shen
- GEO accession: GSE191301
- Input category: per-sample 10x Matrix Market triplets
- Special lesson: GSM-prefix triplet matching and treatment metadata validation

## Why this case matters

This case contains six per-sample 10x Matrix Market triplets.

The biological design is treatment-aware:

- 3 Pre-NACT samples
- 3 Post-NACT samples

The tool must correctly match matrix, features and barcodes by GSM prefix.

The original script risks reusing the first triplet for all samples.

## Detected files

Total GSM files:

18

Matrix files:

6

Feature files:

6

Barcode files:

6

Other files:

0

Complete triplets:

6

## Selected samples

- GSM5743307 / Pre-NACT1A / No treatment / Peritoneum / 6,456 cells
- GSM5743308 / Pre-NACT1B / No treatment / Omentum / 8,372 cells
- GSM5743309 / Pre-NACT1C / No treatment / Ascites / 7,850 cells
- GSM5743310 / Post-NACT1D / Yes treatment / Peritoneum / 5,431 cells
- GSM5743311 / Post-NACT1E / Yes treatment / Omentum / 5,784 cells
- GSM5743312 / Post-NACT1F / Yes treatment / Ascites / 7,887 cells

Total cells:

41,780

Unique gene count:

27,984

## First triplet

Sample:

GSM5743307 / Pre-NACT1A

Files:

- GSM5743307_Pre-NACT1A_matrix.mtx.gz
- GSM5743307_Pre-NACT1A_features.tsv.gz
- GSM5743307_Pre-NACT1A_barcodes.tsv.gz

Matrix shape:

27,984 genes x 6,456 cells

Expected AnnData shape:

6,456 cells x 27,984 genes

Feature table shape:

27,984 rows x 3 columns

Barcode table shape:

6,456 rows x 1 column

## Orientation

The matrix is genes x cells before transpose.

AnnData construction requires transpose.

## Feature table

Feature columns:

- gene_ensembl_id
- gene_symbol
- feature_type

First feature examples:

- ENSG00000243485 / MIR1302-2HG / Gene Expression
- ENSG00000237613 / FAM138A / Gene Expression
- ENSG00000186092 / OR4F5 / Gene Expression
- ENSG00000238009 / AL627309.1 / Gene Expression
- ENSG00000239945 / AL627309.3 / Gene Expression

Gene symbol column index:

1

Duplicated gene symbols per sample:

20

## Original script risk

The original script loops over GSM IDs but selects:

mtx_files[0]
feat_files[0]
bc_files[0]

This can silently assign the first matrix/features/barcodes triplet to every sample.

Correct behavior:

For each GSM, select the matrix, features and barcodes file whose filename starts with that GSM.

Example:

GSM5743307 must use:

- GSM5743307_Pre-NACT1A_matrix.mtx.gz
- GSM5743307_Pre-NACT1A_features.tsv.gz
- GSM5743307_Pre-NACT1A_barcodes.tsv.gz

GSM5743312 must use:

- GSM5743312_Post-NACT1F_matrix.mtx.gz
- GSM5743312_Post-NACT1F_features.tsv.gz
- GSM5743312_Post-NACT1F_barcodes.tsv.gz

## Metadata validation

Pre-NACT samples should have:

tumor_treatment = No

Post-NACT samples should have:

tumor_treatment = Yes

This rule is useful for detecting metadata inconsistencies.

## Correct standardization recipe

1. Download and extract GSE191301_RAW.tar.
2. Detect all GSM files.
3. Classify matrix, features and barcode files.
4. Build triplets by GSM prefix.
5. Validate all triplets are complete.
6. Validate all GSMs exist in curated metadata.
7. Validate Pre/Post NACT treatment metadata.
8. Read each Matrix Market file.
9. Transpose genes x cells to cells x genes.
10. Read barcodes explicitly.
11. Prefix obs_names with GSM accession.
12. Use gene symbol column from features.tsv.gz.
13. Attach curated sample-level metadata.
14. Concatenate selected samples.
15. Generate gene mapping table.
16. Remove invalid genes.
17. Deduplicate gene symbols using sparse aggregation.
18. Validate output.
19. Write standardized h5ad.
20. Delete raw files after processing.

## Reusable rule

Never assign files by list position inside a sample loop.

Always match files by sample identifier.

## Product lesson

This case strengthens the standardizer's triplet matching logic.

It also adds treatment-aware metadata validation for pre/post therapy datasets.
