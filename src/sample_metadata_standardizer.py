from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import pandas as pd

from src.standard_schema import STANDARD_OBS_COLUMNS, COL_SAMPLE_ID


@dataclass
class SampleMetadataStandardizationResult:
    passed: bool
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    evidence: dict[str, Any] = field(default_factory=dict)


def normalize_curated_sample_metadata(
    curated_sample_metadata: dict[str, dict[str, Any]],
) -> pd.DataFrame:
    if not isinstance(curated_sample_metadata, dict):
        raise ValueError("curated_sample_metadata must be a dictionary.")

    records = []

    for sample_id, metadata in curated_sample_metadata.items():
        if not isinstance(metadata, dict):
            raise ValueError(f"Metadata for sample {sample_id} must be a dictionary.")

        record = {COL_SAMPLE_ID: str(sample_id)}
        for column in STANDARD_OBS_COLUMNS:
            if column == COL_SAMPLE_ID:
                continue
            record[column] = metadata.get(column, "Unknown")

        records.append(record)

    df = pd.DataFrame(records)

    if df.empty:
        return pd.DataFrame(columns=STANDARD_OBS_COLUMNS).set_index(COL_SAMPLE_ID, drop=False)

    return df[STANDARD_OBS_COLUMNS].set_index(COL_SAMPLE_ID, drop=False)


def apply_curated_sample_metadata_to_obs(
    obs: pd.DataFrame,
    curated_sample_metadata: dict[str, dict[str, Any]],
    sample_id_column: str = COL_SAMPLE_ID,
    fill_missing: str = "Unknown",
    keep_extra_obs_columns: bool = True,
) -> tuple[pd.DataFrame, SampleMetadataStandardizationResult]:
    errors = []
    warnings = []

    if sample_id_column not in obs.columns:
        errors.append(f"sample_id_column '{sample_id_column}' is missing from obs.")
        return obs, SampleMetadataStandardizationResult(
            passed=False,
            errors=errors,
            warnings=warnings,
            evidence={"obs_columns": list(obs.columns)},
        )

    metadata_df = normalize_curated_sample_metadata(curated_sample_metadata)

    obs = obs.copy()
    obs[sample_id_column] = obs[sample_id_column].astype(str)

    observed_samples = sorted(obs[sample_id_column].dropna().astype(str).unique().tolist())
    curated_samples = sorted(metadata_df.index.astype(str).tolist())

    missing_from_dictionary = sorted(set(observed_samples) - set(curated_samples))
    dictionary_without_cells = sorted(set(curated_samples) - set(observed_samples))

    if missing_from_dictionary:
        errors.append(
            "Observed samples missing from curated_sample_metadata: "
            + ", ".join(missing_from_dictionary)
        )

    if dictionary_without_cells:
        warnings.append(
            "Curated samples without observed cells: "
            + ", ".join(dictionary_without_cells)
        )

    for column in STANDARD_OBS_COLUMNS:
        if column not in obs.columns:
            obs[column] = fill_missing

    for sample_id in observed_samples:
        if sample_id not in metadata_df.index:
            continue

        mask = obs[sample_id_column] == sample_id
        for column in STANDARD_OBS_COLUMNS:
            obs.loc[mask, column] = metadata_df.loc[sample_id, column]

    if keep_extra_obs_columns:
        obs = obs[STANDARD_OBS_COLUMNS + [c for c in obs.columns if c not in STANDARD_OBS_COLUMNS]]
    else:
        obs = obs[STANDARD_OBS_COLUMNS]

    missing_standard_columns = [c for c in STANDARD_OBS_COLUMNS if c not in obs.columns]

    if missing_standard_columns:
        errors.append(
            "Standard obs columns missing after metadata standardization: "
            + ", ".join(missing_standard_columns)
        )

    return obs, SampleMetadataStandardizationResult(
        passed=len(errors) == 0,
        errors=errors,
        warnings=warnings,
        evidence={
            "sample_id_column": sample_id_column,
            "n_observed_samples": len(observed_samples),
            "n_curated_samples": len(curated_samples),
            "observed_samples": observed_samples,
            "curated_samples": curated_samples,
            "missing_from_dictionary": missing_from_dictionary,
            "dictionary_without_cells": dictionary_without_cells,
            "standard_obs_columns": STANDARD_OBS_COLUMNS,
        },
    )


def validate_standard_obs_schema(obs: pd.DataFrame) -> SampleMetadataStandardizationResult:
    missing = [column for column in STANDARD_OBS_COLUMNS if column not in obs.columns]

    return SampleMetadataStandardizationResult(
        passed=len(missing) == 0,
        errors=[] if not missing else ["Missing standard obs columns: " + ", ".join(missing)],
        warnings=[],
        evidence={
            "obs_columns": list(obs.columns),
            "standard_obs_columns": STANDARD_OBS_COLUMNS,
            "missing": missing,
        },
    )
