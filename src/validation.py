from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Mapping
import math

import numpy as np


@dataclass
class ValidationResult:
    check_label: str
    passed: bool
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    evidence: dict[str, Any] = field(default_factory=dict)


def _is_sparse_matrix(x: Any) -> bool:
    return hasattr(x, "tocoo") and hasattr(x, "data")


def _safe_shape(x: Any) -> tuple[int, int] | None:
    shape = getattr(x, "shape", None)
    if shape is None or len(shape) != 2:
        return None
    return int(shape[0]), int(shape[1])


def check_matrix_shape(
    x: Any,
    expected_n_obs: int | None = None,
    expected_n_vars: int | None = None,
) -> ValidationResult:
    shape = _safe_shape(x)
    errors = []
    warnings = []

    if shape is None:
        errors.append("Matrix does not expose a valid 2D shape.")
        return ValidationResult(
            check_label="matrix_shape",
            passed=False,
            errors=errors,
            warnings=warnings,
            evidence={"shape": None},
        )

    n_obs, n_vars = shape

    if n_obs == 0 or n_vars == 0:
        errors.append("Matrix has zero cells or zero genes.")

    if expected_n_obs is not None and n_obs != expected_n_obs:
        errors.append(f"Expected {expected_n_obs} observations, found {n_obs}.")

    if expected_n_vars is not None and n_vars != expected_n_vars:
        errors.append(f"Expected {expected_n_vars} variables, found {n_vars}.")

    return ValidationResult(
        check_label="matrix_shape",
        passed=len(errors) == 0,
        errors=errors,
        warnings=warnings,
        evidence={
            "n_obs": n_obs,
            "n_vars": n_vars,
            "expected_n_obs": expected_n_obs,
            "expected_n_vars": expected_n_vars,
        },
    )


def check_obs_var_alignment(
    x: Any,
    obs_names: list[str],
    var_names: list[str],
) -> ValidationResult:
    shape = _safe_shape(x)
    errors = []

    if shape is None:
        errors.append("Matrix does not expose a valid 2D shape.")
        return ValidationResult(
            check_label="obs_var_alignment",
            passed=False,
            errors=errors,
            evidence={"shape": None},
        )

    n_obs, n_vars = shape

    if len(obs_names) != n_obs:
        errors.append(f"obs_names length {len(obs_names)} does not match matrix rows {n_obs}.")

    if len(var_names) != n_vars:
        errors.append(f"var_names length {len(var_names)} does not match matrix columns {n_vars}.")

    return ValidationResult(
        check_label="obs_var_alignment",
        passed=len(errors) == 0,
        errors=errors,
        evidence={
            "matrix_shape": shape,
            "n_obs_names": len(obs_names),
            "n_var_names": len(var_names),
        },
    )


def check_duplicate_names(
    names: list[str],
    axis_label: str,
    allow_duplicates: bool = False,
) -> ValidationResult:
    seen = set()
    duplicates = set()

    for name in names:
        if name in seen:
            duplicates.add(name)
        seen.add(name)

    errors = []
    warnings = []

    if duplicates and not allow_duplicates:
        errors.append(f"Duplicate {axis_label} names detected.")
    elif duplicates:
        warnings.append(f"Duplicate {axis_label} names detected.")

    return ValidationResult(
        check_label=f"duplicate_{axis_label}_names",
        passed=len(errors) == 0,
        errors=errors,
        warnings=warnings,
        evidence={
            "n_names": len(names),
            "n_unique_names": len(seen),
            "n_duplicates": len(duplicates),
            "duplicate_examples": sorted(duplicates)[:10],
        },
    )


def check_nan_sparse_safe(x: Any) -> ValidationResult:
    errors = []
    warnings = []

    shape = _safe_shape(x)
    is_sparse = _is_sparse_matrix(x)

    if shape is None:
        errors.append("Matrix does not expose a valid 2D shape.")
        return ValidationResult(
            check_label="nan_sparse_safe",
            passed=False,
            errors=errors,
            evidence={"shape": None},
        )

    if is_sparse:
        data = x.data
        has_nan = bool(np.isnan(data).any()) if data.size else False
        checked_values = int(data.size)
        check_strategy = "sparse_data_only"
    else:
        arr = np.asarray(x)
        has_nan = bool(np.isnan(arr).any()) if arr.size else False
        checked_values = int(arr.size)
        check_strategy = "dense_array"

    if has_nan:
        errors.append("NaN values detected in matrix.")

    if not is_sparse and checked_values > 50_000_000:
        warnings.append("Dense matrix is large. NaN check may be memory-heavy.")

    return ValidationResult(
        check_label="nan_sparse_safe",
        passed=len(errors) == 0,
        errors=errors,
        warnings=warnings,
        evidence={
            "shape": shape,
            "is_sparse": is_sparse,
            "checked_values": checked_values,
            "check_strategy": check_strategy,
            "has_nan": has_nan,
        },
    )


def check_expected_sample_counts(
    observed_counts: Mapping[str, int],
    expected_counts: Mapping[str, int] | None = None,
    required_samples: list[str] | None = None,
) -> ValidationResult:
    errors = []
    warnings = []

    observed = dict(observed_counts)
    expected = dict(expected_counts or {})

    if required_samples:
        missing_required = sorted(set(required_samples) - set(observed))
        if missing_required:
            errors.append("Required samples missing from observed counts: " + ", ".join(missing_required))

    if expected:
        for sample_id, expected_count in expected.items():
            observed_count = observed.get(sample_id)
            if observed_count is None:
                errors.append(f"Expected sample {sample_id} is missing from observed counts.")
            elif observed_count != expected_count:
                errors.append(
                    f"Sample {sample_id} expected {expected_count} cells, found {observed_count}."
                )

        unexpected_samples = sorted(set(observed) - set(expected))
        if unexpected_samples:
            warnings.append("Observed samples not present in expected counts: " + ", ".join(unexpected_samples))

    zero_count_samples = sorted([sample_id for sample_id, count in observed.items() if count == 0])
    if zero_count_samples:
        errors.append("Samples with zero cells detected: " + ", ".join(zero_count_samples))

    return ValidationResult(
        check_label="expected_sample_counts",
        passed=len(errors) == 0,
        errors=errors,
        warnings=warnings,
        evidence={
            "observed_counts": observed,
            "expected_counts": expected,
            "required_samples": required_samples or [],
            "total_observed_cells": sum(observed.values()),
        },
    )


def merge_validation_results(results: list[ValidationResult]) -> ValidationResult:
    errors = []
    warnings = []
    evidence = {}

    for result in results:
        errors.extend([f"{result.check_label}: {error}" for error in result.errors])
        warnings.extend([f"{result.check_label}: {warning}" for warning in result.warnings])
        evidence[result.check_label] = result.evidence

    return ValidationResult(
        check_label="validation_suite",
        passed=len(errors) == 0,
        errors=errors,
        warnings=warnings,
        evidence=evidence,
    )


def print_validation_result(result: ValidationResult) -> None:
    print("=== Validation result ===")
    print("Check:", result.check_label)
    print("Passed:", result.passed)

    if result.errors:
        print("\nErrors:")
        for error in result.errors:
            print("-", error)

    if result.warnings:
        print("\nWarnings:")
        for warning in result.warnings:
            print("-", warning)

    print("\nEvidence:")
    for key, value in result.evidence.items():
        print(f"- {key}: {value}")
