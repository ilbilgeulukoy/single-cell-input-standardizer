from src.count_table_standardizer import (
    read_count_table,
    infer_count_table_structure,
    build_anndata_from_gene_by_cell_table,
    check_nan_in_adata,
)

filepath = "examples/example_gene_by_cell_counts.csv"

df = read_count_table(filepath)
print("DataFrame:")
print(df)

structure = infer_count_table_structure(df)
print("\nDetected structure:")
print(structure)

adata = build_anndata_from_gene_by_cell_table(
    df=df,
    sample_id="sample_example",
    sample_metadata={
        "dataset_id": "demo_dataset",
        "condition": "example_condition",
    },
    source_file=filepath,
)

print("\nAnnData:")
print(adata)
print("\nobs:")
print(adata.obs)
print("\nvar:")
print(adata.var)

print("\nNaN check:")
print(check_nan_in_adata(adata))
