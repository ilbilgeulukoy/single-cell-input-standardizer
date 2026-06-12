from __future__ import annotations

from dataclasses import dataclass, field
from typing import Iterable, Mapping, Any

from src.sample_id_parser import (
    FilenameSampleParse,
    parse_filename_sample_id,
    validate_filename_against_metadata,
)


@dataclass
class MetadataCheckResult:
    check_label: str
    passed: bool
    errors: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    evidence: dict[str, Any] = field(default_factory=dict)


def check_gsm_coverage(
    detected_gsms: Iterable[str],
    metadata_by_gsm: Mapping[str, Mapping[str, Any]],
) -> MetadataCheckResult:
    detected = set(detected_gsms)
    metadata = set(metadata_by_gsm)

    missing_from_metadata = sorted(detected - metadata)
    metadata_without_file = sorted(metadata - detected)

    errors = []
    warnings = []

    if missing_from_metadata:
        errors.append(
            "Detected GSMs missing from metadata: " + ", ".join(missing_from_metadata)
        )

    if metadata_without_file:
        warnings.append(
            "Metadata GSMs missing from detected files: " + ", ".join(metadata_without_file)
        )

    return MetadataCheckResult(
        check_label="gsm_coverage",
        passed=len(errors) == 0,
        errors=errors,
        warnings=warnings,
        evidence={
            "n_detected_gsms": len(detected),
            "n_metadata_gsms": len(metadata),
            "missing_from_metadata": missing_from_metadata,
            "metadata_without_file": metadata_without_file,
        },
    )


def check_required_fields(
    metadata_by_gsm: Mapping[str, Mapping[str, Any]],
    required_fields: Iterable[str],
) -> MetadataCheckResult:
    required = list(required_fields)
    errors = []
    warnings = []
    missing_by_gsm: dict[str, list[str]] = {}
    empty_by_gsm: dict[str, list[str]] = {}

    for gsm, metadata in metadata_by_gsm.items():
        missing = [field for field in required if field not in metadata]
        empty = [
            field
            for field in required
            if field in metadata and metadata[field] in {None, ""}
        ]

        if missing:
            missing_by_gsm[gsm] = missing
        if empty:
            empty_by_gsm[gsm] = empty

    if missing_by_gsm:
        errors.append("Some metadata records are missing required fields.")

    if empty_by_gsm:
        warnings.append("Some metadata records contain empty required fields.")

    return MetadataCheckResult(
        check_label="required_fields",
        passed=len(errors) == 0,
        errors=errors,
        warnings=warnings,
        evidence={
            "required_fields": required,
            "missing_by_gsm": missing_by_gsm,
            "empty_by_gsm": empty_by_gsm,
        },
    )


def check_duplicate_values(
    metadata_by_gsm: Mapping[str, Mapping[str, Any]],
    field_name: str,
    allow_duplicates: bool = True,
) -> MetadataCheckResult:
    value_to_gsms: dict[str, list[str]] = {}

    for gsm, metadata in metadata_by_gsm.items():
        value = metadata.get(field_name)
        if value is None:
            continue
        value_to_gsms.setdefault(str(value), []).append(gsm)

    duplicates = {
        value: gsms
        for value, gsms in value_to_gsms.items()
        if len(gsms) > 1
    }

    errors = []
    warnings = []

    if duplicates and not allow_duplicates:
        errors.append(f"Duplicate values detected for field '{field_name}'.")
    elif duplicates:
        warnings.append(f"Duplicate values detected for field '{field_name}'.")

    return MetadataCheckResult(
        check_label=f"duplicate_values:{field_name}",
        passed=len(errors) == 0,
        errors=errors,
        warnings=warnings,
        evidence={
            "field_name": field_name,
            "duplicates": duplicates,
            "n_unique_values": len(value_to_gsms),
        },
    )


def check_filename_metadata_consistency(
    filenames: Iterable[str],
    metadata_by_gsm: Mapping[str, Mapping[str, Any]],
) -> MetadataCheckResult:
    errors = []
    warnings = []
    parsed_files: dict[str, dict[str, Any]] = {}

    for filename in filenames:
        parsed = parse_filename_sample_id(filename)

        parsed_files[filename] = {
            "gsm_id": parsed.gsm_id,
            "parser_label": parsed.parser_label,
            "inferred_patient_id": parsed.inferred_patient_id,
            "treatment_label": parsed.treatment_label,
            "site_label": parsed.site_label,
            "warnings": parsed.warnings,
        }

        if parsed.warnings:
            warnings.extend([f"{filename}: {warning}" for warning in parsed.warnings])

        if parsed.gsm_id is None:
            warnings.append(f"{filename}: no GSM ID could be parsed from filename.")
            continue

        metadata = metadata_by_gsm.get(parsed.gsm_id)
        if metadata is None:
            errors.append(f"{filename}: parsed GSM {parsed.gsm_id} not found in metadata.")
            continue

        mismatch_warnings = validate_filename_against_metadata(parsed, metadata)
        warnings.extend([f"{filename}: {warning}" for warning in mismatch_warnings])

    return MetadataCheckResult(
        check_label="filename_metadata_consistency",
        passed=len(errors) == 0,
        errors=errors,
        warnings=warnings,
        evidence={"parsed_files": parsed_files},
    )


def check_cell_metadata_coverage(
    matrix_cell_ids: Iterable[str],
    metadata_cell_ids: Iterable[str],
    require_all_matrix_cells_in_metadata: bool = False,
) -> MetadataCheckResult:
    matrix_cells = set(matrix_cell_ids)
    metadata_cells = set(metadata_cell_ids)

    matrix_without_metadata = sorted(matrix_cells - metadata_cells)
    metadata_without_matrix = sorted(metadata_cells - matrix_cells)

    errors = []
    warnings = []

    if matrix_without_metadata:
        message = (
            f"{len(matrix_without_metadata)} matrix cells are missing from metadata."
        )
        if require_all_matrix_cells_in_metadata:
            errors.append(message)
        else:
            warnings.append(message)

    if metadata_without_matrix:
        warnings.append(
            f"{len(metadata_without_matrix)} metadata cells are missing from matrix."
        )

    return MetadataCheckResult(
        check_label="cell_metadata_coverage",
        passed=len(errors) == 0,
        errors=errors,
        warnings=warnings,
        evidence={
            "n_matrix_cells": len(matrix_cells),
            "n_metadata_cells": len(metadata_cells),
            "n_common_cells": len(matrix_cells & metadata_cells),
            "n_matrix_without_metadata": len(matrix_without_metadata),
            "n_metadata_without_matrix": len(metadata_without_matrix),
            "matrix_without_metadata_examples": matrix_without_metadata[:10],
            "metadata_without_matrix_examples": metadata_without_matrix[:10],
        },
    )


def check_project_scope_keywords(
    metadata_by_gsm: Mapping[str, Mapping[str, Any]],
    expected_keywords: Iterable[str],
    fields_to_check: Iterable[str],
) -> MetadataCheckResult:
    expected = [keyword.lower() for keyword in expected_keywords]
    fields = list(fields_to_check)

    warnings = []
    suspicious_records: dict[str, dict[str, Any]] = {}

    for gsm, metadata in metadata_by_gsm.items():
        combined_values = []

        for field in fields:
            value = metadata.get(field)
            if value is not None:
                combined_values.append(str(value).lower())

        combined_text = " ".join(combined_values)

        if expected and not any(keyword in combined_text for keyword in expected):
            suspicious_records[gsm] = {
                field: metadata.get(field)
                for field in fields
            }

    if suspicious_records:
        warnings.append(
            "Some metadata records do not contain expected project-scope keywords."
        )

    return MetadataCheckResult(
        check_label="project_scope_keywords",
        passed=True,
        errors=[],
        warnings=warnings,
        evidence={
            "expected_keywords": list(expected_keywords),
            "fields_checked": fields,
            "suspicious_records": suspicious_records,
        },
    )


def merge_check_results(results: Iterable[MetadataCheckResult]) -> MetadataCheckResult:
    result_list = list(results)
    errors = []
    warnings = []
    evidence = {}

    for result in result_list:
        errors.extend([f"{result.check_label}: {error}" for error in result.errors])
        warnings.extend([f"{result.check_label}: {warning}" for warning in result.warnings])
        evidence[result.check_label] = result.evidence

    return MetadataCheckResult(
        check_label="metadata_check_suite",
        passed=len(errors) == 0,
        errors=errors,
        warnings=warnings,
        evidence=evidence,
    )


def print_metadata_check_result(result: MetadataCheckResult) -> None:
    print("=== Metadata check result ===")
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
