import numpy as np
from scipy import sparse

from src.validation import (
    check_matrix_shape,
    check_obs_var_alignment,
    check_duplicate_names,
    check_nan_sparse_safe,
    check_expected_sample_counts,
    merge_validation_results,
)


def test_matrix_shape_passes():
    x = np.zeros((3, 4))
    result = check_matrix_shape(x, expected_n_obs=3, expected_n_vars=4)

    print("\n=== Matrix shape passes ===")
    print(result)

    assert result.passed is True


def test_matrix_shape_fails():
    x = np.zeros((3, 4))
    result = check_matrix_shape(x, expected_n_obs=4, expected_n_vars=4)

    print("\n=== Matrix shape fails ===")
    print(result)

    assert result.passed is False


def test_obs_var_alignment():
    x = np.zeros((2, 3))
    result = check_obs_var_alignment(
        x,
        obs_names=["cell1", "cell2"],
        var_names=["gene1", "gene2", "gene3"],
    )

    print("\n=== Obs var alignment ===")
    print(result)

    assert result.passed is True


def test_duplicate_obs_names_fails():
    result = check_duplicate_names(
        ["cell1", "cell1", "cell2"],
        axis_label="obs",
        allow_duplicates=False,
    )

    print("\n=== Duplicate obs names ===")
    print(result)

    assert result.passed is False
    assert result.evidence["n_duplicates"] == 1


def test_duplicate_var_names_warning_allowed():
    result = check_duplicate_names(
        ["gene1", "gene1", "gene2"],
        axis_label="var",
        allow_duplicates=True,
    )

    print("\n=== Duplicate var names warning ===")
    print(result)

    assert result.passed is True
    assert len(result.warnings) == 1


def test_nan_sparse_safe_dense():
    x = np.array([[1.0, 0.0], [np.nan, 2.0]])
    result = check_nan_sparse_safe(x)

    print("\n=== NaN dense ===")
    print(result)

    assert result.passed is False


def test_nan_sparse_safe_sparse():
    x = sparse.csr_matrix(np.array([[1.0, 0.0], [0.0, 2.0]]))
    result = check_nan_sparse_safe(x)

    print("\n=== NaN sparse ===")
    print(result)

    assert result.passed is True
    assert result.evidence["check_strategy"] == "sparse_data_only"


def test_expected_sample_counts():
    observed = {"GSM1": 10, "GSM2": 20}
    expected = {"GSM1": 10, "GSM2": 20}

    result = check_expected_sample_counts(observed, expected)

    print("\n=== Expected sample counts ===")
    print(result)

    assert result.passed is True


def test_expected_sample_counts_fails():
    observed = {"GSM1": 10, "GSM2": 19}
    expected = {"GSM1": 10, "GSM2": 20}

    result = check_expected_sample_counts(observed, expected)

    print("\n=== Expected sample counts fails ===")
    print(result)

    assert result.passed is False


def test_merge_validation_results():
    x = np.zeros((2, 2))

    suite = merge_validation_results(
        [
            check_matrix_shape(x, expected_n_obs=2, expected_n_vars=2),
            check_duplicate_names(["cell1", "cell1"], "obs"),
        ]
    )

    print("\n=== Validation suite ===")
    print(suite)

    assert suite.passed is False
    assert len(suite.errors) == 1


def main():
    test_matrix_shape_passes()
    test_matrix_shape_fails()
    test_obs_var_alignment()
    test_duplicate_obs_names_fails()
    test_duplicate_var_names_warning_allowed()
    test_nan_sparse_safe_dense()
    test_nan_sparse_safe_sparse()
    test_expected_sample_counts()
    test_expected_sample_counts_fails()
    test_merge_validation_results()

    print("\nAll validation tests passed.")


if __name__ == "__main__":
    main()
