from src.sample_id_parser import (
    parse_filename_sample_id,
    summarize_barcode_suffixes,
    validate_filename_against_metadata,
)


def test_loret_filename():
    parsed = parse_filename_sample_id(
        "GSM6049610_1_N_OT_PT1_filtered_gene_bc_matrices_h5.h5"
    )

    print("\n=== Loret filename ===")
    print(parsed)

    assert parsed.parser_label == "loret_style_filename"
    assert parsed.gsm_id == "GSM6049610"
    assert parsed.sample_order == "1"
    assert parsed.treatment_label == "N"
    assert parsed.site_label == "OT"
    assert parsed.patient_label == "PT1"
    assert parsed.inferred_patient_id == "N_OT_PT1"

    metadata = {
        "patient_id": "N_OT_PT1",
        "tumor_treatment": "No",
        "metastasis_site": "Omentum",
    }

    warnings = validate_filename_against_metadata(parsed, metadata)
    assert warnings == []


def test_shen_filename():
    parsed = parse_filename_sample_id(
        "GSM5743307_Pre-NACT1A_matrix.mtx.gz"
    )

    print("\n=== Shen filename ===")
    print(parsed)

    assert parsed.parser_label == "shen_style_filename"
    assert parsed.gsm_id == "GSM5743307"
    assert parsed.inferred_patient_id == "Pre-NACT1A"
    assert parsed.treatment_label == "Pre"
    assert parsed.patient_label == "1A"


def test_simple_gsm_filename():
    parsed = parse_filename_sample_id(
        "GSM6506105_counts_Y2.txt.gz"
    )

    print("\n=== Simple GSM filename ===")
    print(parsed)

    assert parsed.parser_label == "simple_gsm_prefixed_filename"
    assert parsed.gsm_id == "GSM6506105"
    assert parsed.inferred_patient_id == "counts_Y2"


def test_barcode_suffixes():
    result = summarize_barcode_suffixes(
        [
            '"AAACCCAAGACGAGCT.1_1"',
            '"AAACCCATCTGGTGGC.1_1"',
            "AAACGAAAGAGAGAAC.1_2",
            "AAACGAAAGCCGGATA.1_10",
        ]
    )

    print("\n=== Barcode suffixes ===")
    print(result)

    assert result.parser_label == "multiple_suffixes_detected"
    assert result.suffix_counts["1"] == 2
    assert result.suffix_counts["2"] == 1
    assert result.suffix_counts["10"] == 1


def test_metadata_mismatch_warning():
    parsed = parse_filename_sample_id(
        "GSM6049610_1_N_OT_PT1_filtered_gene_bc_matrices_h5.h5"
    )
    metadata = {
        "patient_id": "T_OT_PT1",
        "tumor_treatment": "Yes",
        "metastasis_site": "Ascites",
    }

    warnings = validate_filename_against_metadata(parsed, metadata)

    print("\n=== Metadata mismatch warnings ===")
    for warning in warnings:
        print(warning)

    assert len(warnings) == 3


def main():
    test_loret_filename()
    test_shen_filename()
    test_simple_gsm_filename()
    test_barcode_suffixes()
    test_metadata_mismatch_warning()
    print("\nAll sample ID parser tests passed.")


if __name__ == "__main__":
    main()
