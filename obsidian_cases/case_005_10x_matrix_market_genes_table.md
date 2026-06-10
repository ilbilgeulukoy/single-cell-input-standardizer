# Case 005 - 10x-like Matrix Market triplets with genes.tsv.gz gene table

## Case identity

- Case ID: case_005_2020_geistlinger_gse154600
- Source example: GSE154600
- Dataset context: high-grade serous ovarian cancer single-cell RNA-seq
- Technology: 10x Genomics
- Input category: per-sample 10x-like Matrix Market triplets
- Tool lesson: support genes.tsv.gz as an alternative to features.tsv.gz

## Input pattern

The GEO archive contains 5 samples and 15 GSM files.

Each sample has one triplet:

- matrix.mtx.gz
- genes.tsv.gz
- barcodes.tsv.gz

Detected files:

- GSM4675273_T59_matrix.mtx.gz
- GSM4675273_T59_genes.tsv.gz
- GSM4675273_T59_barcodes.tsv.gz
- GSM4675274_T76_matrix.mtx.gz
- GSM4675274_T76_genes.tsv.gz
- GSM4675274_T76_barcodes.tsv.gz
- GSM4675275_T77_matrix.mtx.gz
- GSM4675275_T77_genes.tsv.gz
- GSM4675275_T77_barcodes.tsv.gz
- GSM4675276_T89_matrix.mtx.gz
- GSM4675276_T89_genes.tsv.gz
- GSM4675276_T89_barcodes.tsv.gz
- GSM4675277_T90_matrix.mtx.gz
- GSM4675277_T90_genes.tsv.gz
- GSM4675277_T90_barcodes.tsv.gz

## First sample inspection

For GSM4675273:

- matrix shape: 33538 x 16280
- gene table shape: 33538 x 3
- barcodes shape: 16280 x 1

This means:

- rows of matrix match genes
- columns of matrix match barcodes
- raw orientation is genes x cells
- AnnData requires cells x genes

Therefore, the matrix must be transposed.

## Gene table variant

This case uses:

genes.tsv.gz

instead of:

features.tsv.gz

The gene table still has three columns:

1. Ensembl gene ID
2. gene symbol
3. feature type

Example rows:

- ENSG00000243485 / MIR1302-2HG / Gene Expression
- ENSG00000237613 / FAM138A / Gene Expression
- ENSG00000186092 / OR4F5 / Gene Expression

## Main problem

This case is a Matrix Market triplet dataset, but the gene annotation file is named `genes.tsv.gz`.

A robust detector should not only search for `features.tsv.gz`.

It should accept both:

- features.tsv.gz
- genes.tsv.gz

## Correct AnnData recipe

1. Read matrix.mtx.gz as sparse Matrix Market.
2. Read barcodes.tsv.gz explicitly.
3. Read genes.tsv.gz as the gene table.
4. Check matrix shape against gene table and barcodes.
5. Transpose matrix if rows match genes and columns match barcodes.
6. Assign GSM-prefixed real barcodes to obs_names.
7. Assign gene symbols from column 1 to var_names.
8. Add Ensembl IDs and feature types to var.
9. Detect duplicated gene symbols.
10. Add sample-level metadata.
11. Concatenate all samples.
12. Apply gene mapping and cleanup.
13. Merge duplicated gene symbols using sparse aggregation.
14. Validate missing values.
15. Write standardized h5ad.
16. Delete temporary raw files.

## Reusable rule

A 10x-like Matrix Market detector should define the gene table role as:

features.tsv.gz OR genes.tsv.gz

If column 1 exists, use it as gene symbols. Otherwise, fall back to column 0.

## Product lesson

This case extends:

`tenx_mtx_standardizer`

It teaches the standardizer to support older or exported 10x-like datasets where gene annotations are stored in genes.tsv.gz rather than features.tsv.gz.
