from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml


COMMON_REQUIRED_FIELDS = [
    "dataset_id",
    "script_archetype",
    "curated_sample_metadata",
]


ARCHETYPE_REQUIRED_FIELDS = {
    "global_matrix_with_cell_metadata": [
        "input_format",
        "matrix_path",
        "output_h5ad",
    ],
    "per_sample_count_tables_gene_by_cell": [
        "downloaded_dir",
        "sample_file_pattern",
        "sample_id_regex",
    ],
    "tenx_mtx_triplet_single_or_multi_sample": [
        "downloaded_dir",
        "sample_id_regex",
    ],
}


def load_recipe(recipe_path: str | Path) -> dict[str, Any]:
    path = Path(recipe_path)

    if not path.exists():
        raise FileNotFoundError(f"Recipe file not found: {path}")

    recipe = yaml.safe_load(path.read_text(encoding="utf-8"))

    if recipe is None:
        raise ValueError(f"Recipe is empty: {path}")

    missing_common = [
        field for field in COMMON_REQUIRED_FIELDS
        if field not in recipe
    ]

    if missing_common:
        raise ValueError(f"Recipe is missing required fields: {missing_common}")

    archetype = recipe["script_archetype"]

    if archetype not in ARCHETYPE_REQUIRED_FIELDS:
        supported = ", ".join(sorted(ARCHETYPE_REQUIRED_FIELDS))
        raise ValueError(
            f"Unsupported script_archetype: {archetype}. "
            f"Supported archetypes: {supported}"
        )

    missing_archetype = [
        field for field in ARCHETYPE_REQUIRED_FIELDS[archetype]
        if field not in recipe
    ]

    if missing_archetype:
        raise ValueError(
            f"Recipe for archetype '{archetype}' is missing required fields: "
            f"{missing_archetype}"
        )

    return recipe


def get_obs_mappings(recipe: dict[str, Any]) -> dict[str, str]:
    return recipe.get("obs_mappings", {})


def validate_expected_obs_columns(
    observed_columns: list[str],
    expected_columns: list[str] | None,
) -> list[str]:
    if expected_columns is None:
        return []

    observed = set(observed_columns)
    return [
        column for column in expected_columns
        if column not in observed
    ]
