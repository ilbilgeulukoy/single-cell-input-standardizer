# Case 003 - TXT count tables with gene mapping and duplicate gene aggregation

## Case identity

- Case ID: case_003_2020_tongton_gse130000
- Source example: GSE130000
- Dataset context: ovarian cancer single-cell RNA-seq
- Technology: Drop-seq
- Input category: compressed TXT count tables
- Tool lesson: gene mapping and duplicate gene aggregation after AnnData construction

## Input pattern

The GEO supplementary archive contains GSM-level `.txt.gz` count tables.

The previous working script loaded each file using:

scanpy.read_text(fp, delimiter="\t")

Then transposed the AnnData object:

adata = adata.T

This indicates that the raw table was treated as genes x cells and had to be converted to cells x genes.

## Main problem

This case has two layers of standardization:

1. Count table orientation and AnnData construction.
2. Gene identifier cleanup and duplicated gene symbol aggregation.

The second layer is new compared with earlier cases.

## Detection logic

A file belongs to this case type when:

- file extension is `.txt.gz`
- file starts with GSM accession
- tab-delimited count table structure is detected
- previous or inferred orientation is genes x cells
- gene names require downstream standardization

## AnnData recipe

1. Read each TXT count table.
2. Transpose to cells x genes.
3. Set obs_names and var_names as strings.
4. Add sample-level metadata.
5. Build one AnnData object per sample.
6. Align common genes.
7. Concatenate samples.
8. Validate missing values.

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

This module should generalize the previous project-specific `SampleGeneMapGenerator` logic.
