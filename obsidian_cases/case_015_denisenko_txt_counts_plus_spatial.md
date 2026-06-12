# Case 015 - Denisenko mixed scRNA-seq TXT counts plus spatial transcriptomics archive

## Case identity

- Case ID: case_015_2024_denisenko_gse211956
- Paper label: 2024_Denisenko
- GEO accession: GSE211956
- Input category: mixed archive with scRNA-seq TXT count tables and spatial Visium files
- Selected modality: single-cell RNA-seq

## Why this case matters

This case contains both scRNA-seq and spatial transcriptomics files.

The selected scRNA-seq files are compressed TXT count tables.

The spatial files are Matrix Market triplets with spatial.zip folders and should be excluded from this scRNA-seq standardization case.

## Detected files

Total GSM files:

37

Selected scRNA-seq TXT count tables:

5

Spatial Visium file groups:

8

Spatial files:

32

## Selected scRNA-seq files

- GSM6506105_counts_Y2.txt.gz
- GSM6506106_counts_Y3.txt.gz
- GSM6506107_counts_Y5.txt.gz
- GSM6506108_counts_MJ10.txt.gz
- GSM6506109_counts_MJ11.txt.gz

## Excluded spatial files

The archive also contains spatial files for GSM6506110 to GSM6506117.

Each spatial sample has:

- barcodes.tsv.gz
- features.tsv.gz
- matrix.mtx.gz
- spatial.zip

These are excluded from the scRNA-seq case.

## TXT count table structure

The compressed TXT files are gene-by-cell count tables.

Rows are genes.

Columns are cell barcodes.

In the pandas preview, the gene names appeared as row index, not as an explicit first column.

Example gene rows:

- AL627309.1
- AL669831.5
- LINC00115
- FAM41C
- NOC2L

Example cell barcode columns:

- AAACGGGTCACTCTTA_3
- AAAGCAAGTGCAACTT_3
- AAAGTAGCAGGTTTCA_3

## Orientation

Raw orientation:

genes x cells

AnnData orientation requires transpose:

cells x genes

## Sample summary

- GSM6506105 / Y2: 328 cells x 19,915 genes
- GSM6506106 / Y3: 6,384 cells x 19,915 genes
- GSM6506107 / Y5: 7,588 cells x 19,915 genes
- GSM6506108 / MJ10: 426 cells x 19,915 genes
- GSM6506109 / MJ11: 1,956 cells x 19,915 genes

Total cells:

16,682

Unique gene count:

19,915

## Metadata strategy

The script uses curated sample-level metadata.

Sample IDs and patient IDs:

- GSM6506105 -> Y2
- GSM6506106 -> Y3
- GSM6506107 -> Y5
- GSM6506108 -> MJ10
- GSM6506109 -> MJ11

Metadata constants:

- dataset_id: GSE211956
- cancer_type: Ovarian
- tumor_site: Primary
- metastasis_site: Unknown
- tumor_treatment: Yes
- cancer_site_origin: Unknown
- tumour_grade: HGSC
- tumour_stage: III-IV
- histological_subtype: Serous
- patient_ethnicity: Unknown

## Original script risks

The original script decompresses all `.txt.gz` files using gzip -dk.

But it later reads the `.txt.gz` files directly.

This decompression step is unnecessary and increases disk usage.

The original script also uses dense NaN checks after concatenation.

This may be memory-heavy for larger datasets.

## Correct standardization recipe

1. Download and extract GSE211956_RAW.tar.
2. Classify all files by modality.
3. Select only GSM6506105-GSM6506109 `.txt.gz` scRNA-seq count tables.
4. Exclude spatial Visium files.
5. Read `.txt.gz` directly.
6. Treat rows/index as genes and columns as cells.
7. Transpose to cells x genes.
8. Prefix cell barcodes with GSM if needed.
9. Add curated sample-level metadata.
10. Validate patient token from filename.
11. Harmonize common genes across samples.
12. Concatenate selected samples.
13. Generate gene mapping table.
14. Remove invalid genes.
15. Deduplicate gene symbols using sparse aggregation.
16. Use sparse-safe NaN checks.
17. Write standardized h5ad.
18. Delete raw files.

## Product lesson

This case teaches the standardizer to handle compressed TXT gene-by-cell count tables in a mixed scRNA-seq plus spatial transcriptomics archive.

It also teaches that gene identifiers may be stored in the row index, not in an explicit first column.

The orientation detector must therefore inspect both columns and row index.
