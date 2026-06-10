import sqlite3
from pathlib import Path


DATABASE_PATH = Path("database/case_library.db")


def print_cases() -> None:
    if not DATABASE_PATH.exists():
        raise FileNotFoundError(
            "Database not found. Run: python scripts/build_case_index.py"
        )

    conn = sqlite3.connect(DATABASE_PATH)

    try:
        rows = conn.execute(
            """
            SELECT
                case_id,
                geo_accession,
                input_format,
                standardizer_module,
                recommended_processing_mode,
                final_shape
            FROM cases
            ORDER BY case_id
            """
        ).fetchall()

        print("=== Case library ===")

        for row in rows:
            print("\nCase ID:", row[0])
            print("GEO:", row[1])
            print("Input format:", row[2])
            print("Standardizer module:", row[3])
            print("Recommended mode:", row[4])
            print("Final shape:", row[5])

    finally:
        conn.close()


if __name__ == "__main__":
    print_cases()
