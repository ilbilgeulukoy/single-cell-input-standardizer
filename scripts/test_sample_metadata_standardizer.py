from pathlib import Path

import pandas as pd

from src.sample_metadata_standardizer import (
    apply_curated_sample_metadata_to_obs,
    validate_standard_obs_schema,
)
from src.standard_schema import STANDARD_OBS_COLUMNS


def test_apply_curated_sample_metadata_to_obs():
    obs = pd.DataFrame(
        {
            "sample_id": ["GSM1", "GSM1", "GSM2"],
            "cell_type": ["T", "B", "Tumor"],
        },
        index=["cell1", "cell2", "cell3"],
    )

    curated = {
        "GSM1": {
            "patient_id": "PT1",
            "dataset_id": "GSE_TEST",
            "cancer_type": "Ovarian",
            "tumor_site": "Primary",
            "metastasis_site": "Unknown",
            "tumor_treatment": "No",
            "cancer_site_origin": "Ovary",
            "tumour_grade": "HGSC",
            "tumour_stage": "IIIC",
            "histological_subtype": "Serous",
            "patient_ethnicity": "Unknown",
        },
        "GSM2": {
            "patient_id": "PT2",
            "dataset_id": "GSE_TEST",
            "cancer_type": "Ovarian",
            "tumor_site": "Metastasis",
            "metastasis_site": "Omentum",
            "tumor_treatment": "Yes",
            "cancer_site_origin": "Ovary",
            "tumour_grade": "HGSC",
            "tumour_stage": "IVA",
            "histological_subtype": "Serous",
            "patient_ethnicity": "Unknown",
        },
    }

    standardized_obs, result = apply_curated_sample_metadata_to_obs(
        obs=obs,
        curated_sample_metadata=curated,
        sample_id_column="sample_id",
    )

    print("\n=== Standardized obs ===")
    print(standardized_obs)
    print(result)

    assert result.passed is True
    assert list(standardized_obs.columns[: len(STANDARD_OBS_COLUMNS)]) == STANDARD_OBS_COLUMNS
    assert standardized_obs.loc["cell1", "patient_id"] == "PT1"
    assert standardized_obs.loc["cell3", "tumor_site"] == "Metastasis"
    assert "cell_type" in standardized_obs.columns

    schema_result = validate_standard_obs_schema(standardized_obs)
    assert schema_result.passed is True


def test_missing_sample_in_dictionary_fails():
    obs = pd.DataFrame(
        {
            "sample_id": ["GSM1", "GSM_MISSING"],
        },
        index=["cell1", "cell2"],
    )

    curated = {
        "GSM1": {
            "patient_id": "PT1",
            "dataset_id": "GSE_TEST",
            "cancer_type": "Ovarian",
            "tumor_site": "Primary",
            "metastasis_site": "Unknown",
            "tumor_treatment": "No",
            "cancer_site_origin": "Ovary",
            "tumour_grade": "HGSC",
            "tumour_stage": "IIIC",
            "histological_subtype": "Serous",
            "patient_ethnicity": "Unknown",
        }
    }

    standardized_obs, result = apply_curated_sample_metadata_to_obs(
        obs=obs,
        curated_sample_metadata=curated,
        sample_id_column="sample_id",
    )

    print("\n=== Missing sample result ===")
    print(result)

    assert result.passed is False
    assert "GSM_MISSING" in result.errors[0]


def main():
    test_apply_curated_sample_metadata_to_obs()
    test_missing_sample_in_dictionary_fails()
    print("\nAll sample metadata standardizer tests passed.")


if __name__ == "__main__":
    main()
