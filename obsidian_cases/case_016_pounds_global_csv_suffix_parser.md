# Case 016 - Global CSV matrix with cell barcode suffix sample IDs

## Case identity

- Case ID: case_016_2025_pounds_gse281120
- Paper label: 2025_Pounds
- GEO accession: GSE281120
- Input category: single global compressed CSV count matrix
- Special lesson: sample assignment from cell barcode suffix and normal sample exclusion

## Why this case matters

This case does not provide one file per sample.

Instead, all cells are stored in one large compressed CSV count matrix.

Sample identity is encoded in the cell barcode suffix.

The standardizer must parse this suffix and map it to GSM accessions.

## Raw file

File:

GSE281120_counts.csv.gz

Size:

733.99 MB

Total cell columns:

57,490

Gene rows:

26,199

Raw matrix orientation:

genes x cells

Expected AnnData orientation:

cells x genes

## CSV structure

The first column contains gene symbols.

Pandas reads the first column as:

Unnamed: 0

The remaining columns are cell barcodes.

Example genes:

- AL627309.1
- AL669831.5
- FAM87B
- LINC00115
- FAM41C

Example cell barcodes:

- AAACCCAAGACGAGCT.1_1
- AAACCCATCTGGTGGC.1_1
- AAACGAAAGAGAGAAC.1_1

## Quoted header issue

Manual gzip header parsing preserved quotes around cell names.

This produced suffixes such as:

- 1"
- 2"
- 10"

Pandas cleaned the column names correctly.

The reusable parser should strip quotes before parsing suffixes.

Correct suffix parsing rule:

1. Strip quotes.
2. Split cell barcode on underscore.
3. Take the final token.

## Cell suffix distribution

- suffix 1: 5,220 cells
- suffix 2: 2,271 cells
- suffix 3: 3,170 cells
- suffix 5: 20,363 cells
- suffix 6: 5,985 cells
- suffix 7: 4,344 cells
- suffix 8: 5,747 cells
- suffix 9: 4,055 cells
- suffix 10: 6,335 cells

Total cells before exclusion:

57,490

## Normal sample exclusion

Suffix 5 is excluded.

Reason:

normal / non-metastasis sample

Excluded cells:

20,363

Final selected tumor cells:

37,127

Final expected shape:

37,127 cells x 26,199 genes

## Suffix to GSM mapping

- suffix 1 -> GSM8611246 / Pt1
- suffix 2 -> GSM8611247 / Pt2
- suffix 3 -> GSM8611248 / Pt3
- suffix 6 -> GSM8611251 / Pt6
- suffix 7 -> GSM8611252 / Pt7
- suffix 8 -> GSM8611254 / Pt8
- suffix 9 -> GSM8611255 / Pt9
- suffix 10 -> GSM8611257 / Pt10

Suffix 4 is absent from the selected mapping.

Suffix 5 is intentionally excluded.

## Per-sample selected cell counts

- GSM8611246 / Pt1: 5,220 cells
- GSM8611247 / Pt2: 2,271 cells
- GSM8611248 / Pt3: 3,170 cells
- GSM8611251 / Pt6: 5,985 cells
- GSM8611252 / Pt7: 4,344 cells
- GSM8611254 / Pt8: 5,747 cells
- GSM8611255 / Pt9: 4,055 cells
- GSM8611257 / Pt10: 6,335 cells

## Original script risks

The original script downloads:

GSE281120_counts.csv.gz

But then searches for a `.csv` file.

It also uses `list_abspath_sample` for decompression, but that variable is not defined.

Decompression is not needed.

The compressed CSV can be read directly.

This case is also large enough that dense reading should trigger a data size policy warning.

## Correct standardization recipe

1. Download GSE281120_counts.csv.gz.
2. Inspect the header first.
3. Read the compressed CSV directly.
4. Treat the first column as gene symbols.
5. Treat all remaining columns as cell barcodes.
6. Strip quotes from cell barcode names.
7. Parse sample suffix after the final underscore.
8. Remove suffix 5.
9. Map remaining suffixes to GSM accessions.
10. Assign sample-level metadata.
11. Transpose genes x cells to cells x genes.
12. Validate selected cell counts.
13. Generate gene mapping table.
14. Remove invalid genes.
15. Deduplicate gene symbols using sparse aggregation.
16. Validate output using memory-safe checks.
17. Write standardized h5ad.
18. Remove raw file.

## Product lesson

This case teaches that sample identity may exist only inside cell barcode names.

The standardizer must support global matrices where the sample split is not file-based.

It also teaches that quoted CSV headers can break naive parsers.
