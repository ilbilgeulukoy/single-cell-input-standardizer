from __future__ import annotations

import sqlite3
from pathlib import Path
from typing import Any

import yaml


CASES_DIR = Path("cases")
DATABASE_DIR = Path("database")
DATABASE_PATH = DATABASE_DIR / "case_library.db"


def flatten_list(value: Any) -> str:
    if value is None:
        return ""
    if isinstance(value, list):
        return "; ".join(str(item) for item in value)
    return str(value)


def init_database(conn: sqlite3.Connection) -> None:
    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS cases (
            case_id TEXT PRIMARY KEY,
            paper_id TEXT,
            geo_accession TEXT,
            project_type TEXT,
            input_format TEXT,
            standardizer_module TEXT,
            recommended_processing_mode TEXT,
            final_shape TEXT,
            risks TEXT
        )
        """
    )

    conn.execute(
        """
        CREATE TABLE IF NOT EXISTS case_files (
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            case_id TEXT,
            file_role TEXT,
            file_pattern TEXT,
            FOREIGN KEY(case_id) REFERENCES cases(case_id)
        )
        """
    )

    conn.commit()


def extract_case_row(case_data: dict[str, Any]) -> dict[str, str]:
    detected = case_data.get("detected_input_format", {}) or {}
    final_output = case_data.get("final_output", {}) or {}
    input_source = case_data.get("input_source", {}) or {}
    tool_design = case_data.get("tool_design_implication", {}) or {}
    data_size_policy = case_data.get("data_size_policy", {}) or {}

    return {
        "case_id": str(case_data.get("case_id", "")),
        "paper_id": str(case_data.get("paper_id", "")),
        "geo_accession": str(case_data.get("geo_accession", "")),
        "project_type": str(case_data.get("project_type", "")),
        "input_format": str(detected.get("format_label", "")),
        "standardizer_module": str(tool_design.get("reusable_module_needed", "")),
        "recommended_processing_mode": str(
            input_source.get("recommended_processing_mode", "")
            or data_size_policy.get("action", "")
        ),
        "final_shape": str(final_output.get("final_shape", "")),
        "risks": flatten_list(case_data.get("potential_risks", [])),
    }


def insert_case(conn: sqlite3.Connection, case_data: dict[str, Any]) -> None:
    row = extract_case_row(case_data)

    conn.execute(
        """
        INSERT OR REPLACE INTO cases (
            case_id,
            paper_id,
            geo_accession,
            project_type,
            input_format,
            standardizer_module,
            recommended_processing_mode,
            final_shape,
            risks
        )
        VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)
        """,
        (
            row["case_id"],
            row["paper_id"],
            row["geo_accession"],
            row["project_type"],
            row["input_format"],
            row["standardizer_module"],
            row["recommended_processing_mode"],
            row["final_shape"],
            row["risks"],
        ),
    )

    conn.execute("DELETE FROM case_files WHERE case_id = ?", (row["case_id"],))

    raw_files = case_data.get("raw_files_detected", []) or []
    for filename in raw_files:
        conn.execute(
            """
            INSERT INTO case_files (case_id, file_role, file_pattern)
            VALUES (?, ?, ?)
            """,
            (row["case_id"], "raw_detected_file", str(filename)),
        )

    expected_roles = case_data.get("expected_rna_file_roles_per_sample", []) or []
    for role in expected_roles:
        conn.execute(
            """
            INSERT INTO case_files (case_id, file_role, file_pattern)
            VALUES (?, ?, ?)
            """,
            (row["case_id"], "expected_rna_role", str(role)),
        )


def build_index() -> None:
    DATABASE_DIR.mkdir(exist_ok=True)

    yaml_files = sorted(CASES_DIR.glob("*.yaml"))

    if not yaml_files:
        print("No YAML case files found.")
        return

    conn = sqlite3.connect(DATABASE_PATH)

    try:
        init_database(conn)

        for yaml_file in yaml_files:
            with yaml_file.open("r", encoding="utf-8") as f:
                case_data = yaml.safe_load(f)

            insert_case(conn, case_data)
            print(f"Indexed: {yaml_file}")

        conn.commit()

        print(f"\nCase library database written to: {DATABASE_PATH}")
        print(f"Indexed cases: {len(yaml_files)}")

    finally:
        conn.close()


if __name__ == "__main__":
    build_index()
