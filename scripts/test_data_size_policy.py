from pathlib import Path
import shutil

from src.data_size_policy import evaluate_data_size


TMP = Path("data/test_data_size_policy")


def reset_tmp():
    if TMP.exists():
        shutil.rmtree(TMP)
    TMP.mkdir(parents=True)


def create_sparse_file(path: Path, size_mb: int):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as handle:
        handle.seek(size_mb * 1024 * 1024 - 1)
        handle.write(b"\0")


def test_small_input():
    reset_tmp()
    create_sparse_file(TMP / "small.csv.gz", 10)

    decision = evaluate_data_size(TMP)

    print("\n=== small input ===")
    print(decision.policy_label)
    assert decision.policy_label == "small_or_moderate_input"
    assert decision.allow_local_processing is True


def test_large_local_input():
    reset_tmp()
    create_sparse_file(TMP / "large_counts.csv.gz", 800)

    decision = evaluate_data_size(TMP)

    print("\n=== large local input ===")
    print(decision.policy_label)
    assert decision.policy_label == "large_local_with_caution"
    assert decision.allow_local_processing is True
    assert decision.allow_decompression is False


def test_server_recommended_input():
    reset_tmp()
    create_sparse_file(TMP / "large_archive.tar", 6 * 1024)

    decision = evaluate_data_size(TMP)

    print("\n=== server recommended input ===")
    print(decision.policy_label)
    assert decision.policy_label == "server_recommended"
    assert decision.allow_local_processing is False
    assert decision.allow_decompression is False


def main():
    test_small_input()
    test_large_local_input()
    test_server_recommended_input()
    print("\nAll data size policy tests passed.")


if __name__ == "__main__":
    main()
