from pathlib import Path
import shutil
import tarfile

from src.nested_archive_extractor import (
    find_archives,
    safe_extract_archive,
)


TMP = Path("data/test_nested_archive_extractor")


def reset_tmp():
    if TMP.exists():
        shutil.rmtree(TMP)
    TMP.mkdir(parents=True)


def write_tar(path: Path, source_dir: Path):
    with tarfile.open(path, "w:gz") as tar:
        for file in source_dir.rglob("*"):
            if file.is_file():
                tar.add(file, arcname=file.relative_to(source_dir))


def test_extract_archive_and_detect_nested():
    reset_tmp()

    inner_dir = TMP / "inner_source"
    inner_dir.mkdir()
    (inner_dir / "matrix.mtx.gz").write_text("dummy matrix")
    (inner_dir / "features.tsv.gz").write_text("dummy features")
    (inner_dir / "barcodes.tsv.gz").write_text("dummy barcodes")

    nested_tar = TMP / "GSM1_Tumor.tar.gz"
    write_tar(nested_tar, inner_dir)

    outer_dir = TMP / "outer_source"
    outer_dir.mkdir()
    shutil.copy(nested_tar, outer_dir / nested_tar.name)

    outer_tar = TMP / "GSE_TEST_RAW.tar.gz"
    write_tar(outer_tar, outer_dir)

    archives = find_archives(TMP)
    assert outer_tar in archives

    result = safe_extract_archive(outer_tar, TMP / "extracted", overwrite=True)

    print("\n=== Nested archive extraction ===")
    print(result)

    assert result.evidence["n_extracted_files"] == 1
    assert result.evidence["n_nested_archives"] == 1
    assert result.warnings


def main():
    test_extract_archive_and_detect_nested()
    print("\nAll nested archive extractor tests passed.")


if __name__ == "__main__":
    main()
