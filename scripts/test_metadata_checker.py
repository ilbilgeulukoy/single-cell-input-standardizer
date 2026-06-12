from src.metadata_checker import (
    check_gsm_coverage,
    check_required_fields,
    check_duplicate_values,
    check_filename_metadata_consistency,
    check_cell_metadata_coverage,
    check_project_scope_keywords,
    merge_check_results,
)


def test_gsm_coverage_passes():
    detected = ["GSM1", "GSM2"]
    metadata = {
        "GSM1": {"patient_id": "S1"},
        "GSM2": {"patient_id": "S2"},
    }

    result = check_gsm_coverage(detected, metadata)

    print("\n=== GSM coverage passes ===")
    print(result)

    assert result.passed is True
    assert result.errors == []


def test_gsm_coverage_missing_metadata_fails():
    detected = ["GSM1", "GSM2", "GSM3"]
    metadata = {
        "GSM1": {"patient_id": "S1"},
        "GSM2": {"patient_id": "S2"},
    }

    result = check_gsm_coverage(detected, metadata)

    print("\n=== GSM coverage missing metadata ===")
    print(result)

    assert result.passed is False
    assert "GSM3" in result.errors[0]


def test_required_fields():
    metadata = {
        "GSM1": {"patient_id": "S1", "tumor_site": "Primary"},
        "GSM2": {"patient_id": "S2"},
    }

    result = check_required_fields(metadata, ["patient_id", "tumor_site"])

    print("\n=== Required fields ===")
    print(result)

    assert result.passed is False
    assert result.evidence["missing_by_gsm"]["GSM2"] == ["tumor_site"]


def test_duplicate_values_warning():
    metadata = {
        "GSM1": {"patient_id": "PT1"},
        "GSM2": {"patient_id": "PT1"},
        "GSM3": {"patient_id": "PT2"},
    }

    result = check_duplicate_values(metadata, "patient_id", allow_duplicates=True)

    print("\n=== Duplicate values warning ===")
    print(result)

    assert result.passed is True
    assert "PT1" in result.evidence["duplicates"]


def test_filename_metadata_consistency_loret():
    filenames = [
        "GSM6049610_1_N_OT_PT1_filtered_gene_bc_matrices_h5.h5",
        "GSM6049611_2_N_A_PT1_filtered_gene_bc_matrices_h5.h5",
    ]

    metadata = {
        "GSM6049610": {
            "patient_id": "N_OT_PT1",
            "tumor_treatment": "No",
            "metastasis_site": "Omentum",
        },
        "GSM6049611": {
            "patient_id": "N_A_PT1",
            "tumor_treatment": "No",
            "metastasis_site": "Ascites",
        },
    }

    result = check_filename_metadata_consistency(filenames, metadata)

    print("\n=== Filename metadata consistency Loret ===")
    print(result)

    assert result.passed is True
    assert result.errors == []
    assert result.warnings == []


def test_filename_metadata_consistency_detects_mismatch():
    filenames = [
        "GSM6049610_1_N_OT_PT1_filtered_gene_bc_matrices_h5.h5",
    ]

    metadata = {
        "GSM6049610": {
            "patient_id": "T_OT_PT1",
            "tumor_treatment": "Yes",
            "metastasis_site": "Ascites",
        },
    }

    result = check_filename_metadata_consistency(filenames, metadata)

    print("\n=== Filename metadata mismatch ===")
    print(result)

    assert result.passed is True
    assert len(result.warnings) == 3


def test_cell_metadata_coverage_ignacio_style_subset_allowed():
    matrix_cells = ["cell1", "cell2", "cell3", "cell4"]
    metadata_cells = ["cell1", "cell2"]

    result = check_cell_metadata_coverage(
        matrix_cells,
        metadata_cells,
        require_all_matrix_cells_in_metadata=False,
    )

    print("\n=== Cell metadata subset allowed ===")
    print(result)

    assert result.passed is True
    assert result.evidence["n_common_cells"] == 2
    assert result.evidence["n_matrix_without_metadata"] == 2


def test_cell_metadata_coverage_strict_fails():
    matrix_cells = ["cell1", "cell2", "cell3"]
    metadata_cells = ["cell1", "cell2"]

    result = check_cell_metadata_coverage(
        matrix_cells,
        metadata_cells,
        require_all_matrix_cells_in_metadata=True,
    )

    print("\n=== Cell metadata strict fails ===")
    print(result)

    assert result.passed is False


def test_project_scope_keywords_warning():
    metadata = {
        "GSM1": {"cancer_type": "Ovarian", "histological_subtype": "Serous"},
        "GSM2": {"cancer_type": "Cervical", "histological_subtype": "Squamous"},
    }

    result = check_project_scope_keywords(
        metadata,
        expected_keywords=["ovarian", "serous"],
        fields_to_check=["cancer_type", "histological_subtype"],
    )

    print("\n=== Project scope keywords ===")
    print(result)

    assert result.passed is True
    assert "GSM2" in result.evidence["suspicious_records"]


def test_merge_check_results():
    metadata = {
        "GSM1": {"patient_id": "S1", "tumor_site": "Primary"},
        "GSM2": {"patient_id": "S2"},
    }

    suite = merge_check_results(
        [
            check_gsm_coverage(["GSM1", "GSM2"], metadata),
            check_required_fields(metadata, ["patient_id", "tumor_site"]),
        ]
    )

    print("\n=== Merged suite ===")
    print(suite)

    assert suite.passed is False
    assert len(suite.errors) == 1


def main():
    test_gsm_coverage_passes()
    test_gsm_coverage_missing_metadata_fails()
    test_required_fields()
    test_duplicate_values_warning()
    test_filename_metadata_consistency_loret()
    test_filename_metadata_consistency_detects_mismatch()
    test_cell_metadata_coverage_ignacio_style_subset_allowed()
    test_cell_metadata_coverage_strict_fails()
    test_project_scope_keywords_warning()
    test_merge_check_results()

    print("\nAll metadata checker tests passed.")


if __name__ == "__main__":
    main()
