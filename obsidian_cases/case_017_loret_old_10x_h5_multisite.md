# Case 017 - Old 10x h5 genome-group layout with multi-site metadata

## Case identity

- Case ID: case_017_2022_loret_gse201047
- Paper label: 2022_Loret
- GEO accession: GSE201047
- Input category: per-sample old 10x h5 files
- Special lesson: old 10x h5 genome-group layout and rich multi-site metadata

## Why this case matters

This case is the final and most complete h5 format case.

It contains 22 per-sample h5 files.

The files are not modern 10x h5 files with a top-level `/matrix` group.

Instead, they use an older 10x h5 layout with a top-level genome group:

hg19

The standardizer must support both layouts.

## Detected files

Total GSM files:

22

Selected h5 files:

22

Other files:

0

All detected GSMs matched the curated metadata dictionary.

## H5 layout

Top-level key:

hg19

Genome group subkeys:

- barcodes
- data
- gene_names
- genes
- indices
- indptr
- shape

This is an old 10x h5 genome-group layout.

The matrix shape is stored in:

hg19/shape

Gene IDs are stored in:

hg19/genes

Gene names are stored in:

hg19/gene_names

Barcodes are stored in:

hg19/barcodes

## Orientation

Raw h5 shape:

32,738 genes x cells

AnnData orientation after reading:

cells x 32,738 genes

## Feature examples

First gene IDs and names:

- ENSG00000243485 / MIR1302-10
- ENSG00000237613 / FAM138A
- ENSG00000186092 / OR4F5
- ENSG00000238009 / RP11-34P13.7
- ENSG00000239945 / RP11-34P13.8

Duplicated gene names per sample:

95

## Sample summary

Total cells:

106,169

Unique gene count:

32,738

Number of samples:

22

## Metadata dimensions

This case has rich sample-level metadata.

Important fields:

- GSM accession
- patient_id label
- N/T label
- anatomical site code
- PT patient code
- tumor_site
- metastasis_site
- tumor_treatment
- cancer_site_origin
- tumour_stage

## Metadata distribution

Tumor site:

- Metastasis: 6 samples
- Primary: 16 samples

Metastasis site:

- Omentum: 11 samples
- Ascites: 5 samples
- Peritoneum: 5 samples
- Bladder: 1 sample

Tumor treatment:

- No: 11 samples
- Yes: 11 samples

Cancer site origin:

- Omentum: 6 samples
- Ovary: 16 samples

Tumour stage:

- IIIc: 15 samples
- IIIb: 7 samples

## Filename token structure

Example filename:

GSM6049610_1_N_OT_PT1_filtered_gene_bc_matrices_h5.h5

Token interpretation:

- token 0: GSM accession
- token 1: numeric sample order
- token 2: N/T label
- token 3: anatomical site code
- token 4: patient code

Curated patient_id contains the full N/T/site/PT label.

Example:

N_OT_PT1

## Correct h5 handling rule

The standardizer should first check whether `/matrix` exists.

If yes:

process modern 10x h5 v3 layout.

If no:

look for a genome group such as `/hg19`.

Then process:

- genome_group/barcodes
- genome_group/data
- genome_group/genes
- genome_group/gene_names
- genome_group/indices
- genome_group/indptr
- genome_group/shape

## Original script risks

The original script calls:

df.isna()

But `df` is not defined.

It parses patient_id from filename token 4, but this parsed variable is not used.

It does not explicitly prefix obs_names with GSM before concatenation.

This may allow barcode collisions across samples.

It also uses dense NaN checks on adata_combined.X, which can be memory-heavy or invalid for sparse matrices.

## Correct standardization recipe

1. Download and extract GSE201047_RAW.tar.
2. Detect 22 per-sample h5 files.
3. Validate all GSMs against curated metadata.
4. Detect h5 layout for each file.
5. Support both `/matrix` and `/hg19` style h5 files.
6. Read old h5 keys from `/hg19`.
7. Interpret matrix as genes x cells.
8. Read with scanpy.read_10x_h5 during final processing.
9. Prefix barcodes with GSM before concatenation.
10. Add curated multi-site sample metadata.
11. Validate N/T, site and PT filename tokens against metadata.
12. Harmonize common genes if required.
13. Concatenate all samples.
14. Generate gene mapping table.
15. Remove invalid genes.
16. Deduplicate gene symbols using sparse aggregation.
17. Convert metastasis_site to string if needed.
18. Validate output using sparse-safe checks.
19. Write standardized h5ad.
20. Delete raw files.

## Product lesson

This case teaches that 10x h5 is not one single layout.

A robust standardizer must support both modern `/matrix` h5 files and older genome-group h5 files such as `/hg19`.

It also teaches multi-dimensional sample metadata validation.
