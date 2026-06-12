from __future__ import annotations

import subprocess
import sys


TEST_SCRIPTS = [
    "scripts/test_file_classifier.py",
    "scripts/test_data_size_policy.py",
    "scripts/test_sample_id_parser.py",
    "scripts/test_metadata_checker.py",
    "scripts/test_validation.py",
    "scripts/test_tenx_h5_standardizer.py",
    "scripts/test_tenx_mtx_standardizer.py",
    "scripts/test_count_table_standardizer.py",
    "scripts/test_global_matrix_metadata_standardizer.py",
    "scripts/test_nested_archive_extractor.py",
    "scripts/test_modality_filter.py",
    "scripts/test_gene_mapping_and_deduplication.py",
]


def main():
    failed = []

    for script in TEST_SCRIPTS:
        print("\n" + "=" * 80)
        print("Running:", script)
        print("=" * 80)

        result = subprocess.run([sys.executable, script])

        if result.returncode != 0:
            failed.append(script)

    print("\n" + "=" * 80)

    if failed:
        print("Failed tests:")
        for script in failed:
            print("-", script)
        raise SystemExit(1)

    print("All tests passed.")


if __name__ == "__main__":
    main()
