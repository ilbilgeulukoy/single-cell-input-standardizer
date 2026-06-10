# Case 006 - Cancer-only dot-separated 10x-like Matrix Market triplets

## Case identity

- Case ID: case_006_2022_junfen_gse184880
- Source example: GSE184880
- Dataset context: high-grade serous ovarian cancer single-cell RNA-seq
- Input category: Cancer-only 10x-like Matrix Market triplets
- Tool lesson: sample filtering and dot-separated Matrix Market role patterns

## Input pattern

The archive contained 36 GSM files in total.

Only files containing `_Cancer` were selected for this case.

The selected cancer subset contained 21 files:

- 7 cancer samples
- 3 files per sample
- matrix file
- genes file
- barcodes file

## Selected files

- GSM5599225_Cancer1.matrix.mtx.gz
- GSM5599225_Cancer1.genes.tsv.gz
- GSM5599225_Cancer1.barcodes.tsv.gz
- GSM5599226_Cancer2.matrix.mtx.gz
- GSM5599226_Cancer2.genes.tsv.gz
- GSM5599226_Cancer2.barcodes.tsv.gz
- GSM5599227_Cancer3.matrix.mtx.gz
- GSM5599227_Cancer3.genes.tsv.gz
- GSM5599227_Cancer3.barcodes.tsv.gz
- GSM5599228_Cancer4.matrix.mtx.gz
- GSM5599228_Cancer4.genes.tsv.gz
- GSM5599228_Cancer4.barcodes.tsv.gz
- GSM5599229_Cancer5.matrix.mtx.gz
- GSM5599229_Cancer5.genes.tsv.gz
- GSM5599229_Cancer5.barcodes.tsv.gz
- GSM5599230_Cancer6.matrix.mtx.gz
- GSM5599230_Cancer6.genes.tsv.gz
- GSM5599230_Cancer6.barcodes.tsv.gz
- GSM5599231_Cancer7.matrix.mtx.gz
- GSM5599231_Cancer7.genes.tsv.gz
- GSM5599231_Cancer7.barcodes.tsv.gz

## Filename grammar

This case uses dot-separated role names:

- `.matrix.mtx.gz`
- `.genes.tsv.gz`
- `.barcodes.tsv.gz`

This differs from earlier cases where the role was encoded using underscores.

## First sample inspection

For GSM5599225:

- matrix shape: 27984 x 8823
- gene table shape: 27984 x 3
- barcodes shape: 8823 x 1

This means:

- matrix rows match genes
- matrix columns match barcodes
- raw orientation is genes x cells
- AnnData requires cells x genes

Therefore, the matrix must be transposed.

## Gene table

The gene table has three columns:

1. Ensembl gene ID
2. gene symbol
3. feature type

The standardizer should use column 1 as the gene symbol when available.

## Main problem

This case adds two new detection requirements:

1. Filter files by sample type keyword, here `_Cancer`.
2. Detect dot-separated Matrix Market role patterns.

A hardcoded detector looking only for `_matrix_`, `_genes_` or `_barcodes_` would miss this case.

## Correct AnnData recipe

1. Extract the GEO archive temporarily.
2. List all GSM files.
3. Filter files containing `_Cancer`.
4. Group files by GSM accession.
5. Detect `.matrix.`, `.genes.`, and `.barcodes.` role patterns.
6. Read matrix.mtx.gz as sparse Matrix Market.
7. Read genes.tsv.gz as the gene table.
8. Read barcodes.tsv.gz explicitly.
9. Check matrix shape against gene table and barcodes.
10. Transpose matrix if genes x cells.
11. Assign GSM-prefixed real barcodes to obs_names.
12. Assign gene symbols from gene table column 1 to var_names.
13. Add sample-level metadata to obs.
14. Concatenate samples.
15. Apply gene mapping and cleanup.
16. Merge duplicated gene symbols using sparse aggregation.
17. Validate X, obs and var.
18. Write standardized h5ad.
19. Delete temporary raw files.

## Reusable rule

The 10x-like detector must be pattern-flexible.

It should detect role names across variants such as:

- `_matrix_`
- `_matrix.`
- `.matrix.`
- `.matrix.mtx.gz`

The same applies to genes/features and barcodes.

## Product lesson

This case extends:

`tenx_mtx_standardizer`

It teaches the tool that filename grammar is not stable across GEO supplementary files.

The product should support user-defined or automatically inferred sample filters, such as Cancer-only selection.
