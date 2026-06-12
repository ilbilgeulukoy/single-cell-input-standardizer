from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any
from collections import Counter

import pandas as pd
import numpy as np
from scipy import sparse


@dataclass
class GeneMappingResult:
    n_input_genes: int
    n_unique_gene_symbols: int
    n_duplicate_gene_symbols: int
    duplicate_examples: list[str]
    mapping_table: pd.DataFrame
    warnings: list[str] = field(default_factory=list)
    evidence: dict[str, Any] = field(default_factory=dict)


@dataclass
class GeneDeduplicationResult:
    matrix: Any
    gene_symbols: list[str]
    aggregation_report: pd.DataFrame
    warnings: list[str] = field(default_factory=list)
    evidence: dict[str, Any] = field(default_factory=dict)


def detect_gene_id_type(values: list[str]) -> str:
    if not values:
        return "unknown"

    ensembl = sum(str(v).startswith(("ENSG", "ENSMUSG")) for v in values)
    if ensembl / len(values) >= 0.5:
        return "ensembl_id"

    return "gene_symbol_or_custom_feature"


def build_gene_mapping_table(
    gene_ids: list[str] | None,
    gene_symbols: list[str],
) -> GeneMappingResult:
    if gene_ids is None:
        gene_ids = [None] * len(gene_symbols)

    if len(gene_ids) != len(gene_symbols):
        raise ValueError("gene_ids and gene_symbols must have the same length.")

    counts = Counter(gene_symbols)
    duplicates = sorted([symbol for symbol, count in counts.items() if count > 1])

    table = pd.DataFrame(
        {
            "gene_id": gene_ids,
            "gene_symbol": gene_symbols,
            "is_duplicate_symbol": [counts[symbol] > 1 for symbol in gene_symbols],
            "feature_id_type": detect_gene_id_type([g for g in gene_ids if g is not None]),
        }
    )

    warnings = []
    if duplicates:
        warnings.append("Duplicate gene symbols detected. Sparse aggregation may be required.")

    return GeneMappingResult(
        n_input_genes=len(gene_symbols),
        n_unique_gene_symbols=len(counts),
        n_duplicate_gene_symbols=len(duplicates),
        duplicate_examples=duplicates[:10],
        mapping_table=table,
        warnings=warnings,
        evidence={
            "feature_id_type": table["feature_id_type"].iloc[0] if len(table) else "unknown",
        },
    )


def aggregate_duplicate_gene_symbols(matrix: Any, gene_symbols: list[str]) -> GeneDeduplicationResult:
    if len(gene_symbols) != matrix.shape[1]:
        raise ValueError("gene_symbols length must match matrix columns.")

    unique_symbols = []
    symbol_to_index = {}

    for symbol in gene_symbols:
        if symbol not in symbol_to_index:
            symbol_to_index[symbol] = len(unique_symbols)
            unique_symbols.append(symbol)

    row_indices = []
    col_indices = []
    data = []

    for old_idx, symbol in enumerate(gene_symbols):
        row_indices.append(old_idx)
        col_indices.append(symbol_to_index[symbol])
        data.append(1)

    aggregation_matrix = sparse.csr_matrix(
        (data, (row_indices, col_indices)),
        shape=(len(gene_symbols), len(unique_symbols)),
    )

    if sparse.issparse(matrix):
        aggregated = matrix @ aggregation_matrix
    else:
        aggregated = sparse.csr_matrix(matrix) @ aggregation_matrix

    counts = Counter(gene_symbols)
    report = pd.DataFrame(
        {
            "gene_symbol": unique_symbols,
            "n_original_features": [counts[symbol] for symbol in unique_symbols],
        }
    )

    duplicated = report[report["n_original_features"] > 1]

    warnings = []
    if len(duplicated):
        warnings.append("Duplicate gene symbols were aggregated by summing counts.")

    return GeneDeduplicationResult(
        matrix=aggregated,
        gene_symbols=unique_symbols,
        aggregation_report=report,
        warnings=warnings,
        evidence={
            "n_input_genes": len(gene_symbols),
            "n_output_genes": len(unique_symbols),
            "n_aggregated_symbols": int(len(duplicated)),
        },
    )


if __name__ == "__main__":
    print("Gene mapping and deduplication utilities.")
