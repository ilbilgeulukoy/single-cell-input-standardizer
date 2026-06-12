import numpy as np
from scipy import sparse

from src.gene_mapping_and_deduplication import (
    build_gene_mapping_table,
    aggregate_duplicate_gene_symbols,
    detect_gene_id_type,
)


def test_detect_gene_id_type():
    assert detect_gene_id_type(["ENSG1", "ENSG2", "GeneA"]) == "ensembl_id"
    assert detect_gene_id_type(["GeneA", "GeneB"]) == "gene_symbol_or_custom_feature"


def test_build_gene_mapping_table():
    result = build_gene_mapping_table(
        gene_ids=["ENSG1", "ENSG2", "ENSG3"],
        gene_symbols=["GeneA", "GeneB", "GeneA"],
    )

    print("\n=== Gene mapping table ===")
    print(result.mapping_table)

    assert result.n_input_genes == 3
    assert result.n_unique_gene_symbols == 2
    assert result.n_duplicate_gene_symbols == 1
    assert "GeneA" in result.duplicate_examples
    assert result.warnings


def test_aggregate_duplicate_gene_symbols_sparse():
    matrix = sparse.csr_matrix(
        np.array(
            [
                [1, 2, 3],
                [4, 5, 6],
            ]
        )
    )

    result = aggregate_duplicate_gene_symbols(
        matrix,
        ["GeneA", "GeneB", "GeneA"],
    )

    print("\n=== Aggregated matrix ===")
    print(result.matrix.toarray())
    print(result.aggregation_report)

    assert result.matrix.shape == (2, 2)
    assert result.gene_symbols == ["GeneA", "GeneB"]
    assert result.matrix.toarray().tolist() == [[4, 2], [10, 5]]
    assert result.evidence["n_aggregated_symbols"] == 1


def main():
    test_detect_gene_id_type()
    test_build_gene_mapping_table()
    test_aggregate_duplicate_gene_symbols_sparse()
    print("\nAll gene mapping and deduplication tests passed.")


if __name__ == "__main__":
    main()
