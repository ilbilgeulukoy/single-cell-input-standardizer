from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
import re
from collections import Counter


@dataclass
class FilenameSampleParse:
    filename: str
    gsm_id: str | None
    sample_order: str | None
    treatment_label: str | None
    site_label: str | None
    patient_label: str | None
    inferred_patient_id: str | None
    parser_label: str
    warnings: list[str] = field(default_factory=list)
    evidence: dict[str, object] = field(default_factory=dict)


@dataclass
class BarcodeSuffixParse:
    total_barcodes: int
    suffix_counts: dict[str, int]
    parser_label: str
    warnings: list[str] = field(default_factory=list)
    evidence: dict[str, object] = field(default_factory=dict)


def strip_compression_suffix(filename: str) -> str:
    suffixes = [
        ".tar.gz",
        ".csv.gz",
        ".txt.gz",
        ".tsv.gz",
        ".mtx.gz",
        ".h5ad.gz",
        ".gz",
        ".h5ad",
        ".h5",
        ".csv",
        ".txt",
        ".tsv",
        ".mtx",
        ".tar",
    ]

    name = filename
    for suffix in suffixes:
        if name.endswith(suffix):
            return name[: -len(suffix)]

    return name


def extract_gsm_id(text: str) -> str | None:
    match = re.search(r"(GSM\d+)", text)
    if match:
        return match.group(1)
    return None


def parse_loret_style_filename(filename: str) -> FilenameSampleParse:
    base = strip_compression_suffix(Path(filename).name)
    tokens = base.split("_")

    warnings: list[str] = []
    gsm_id = extract_gsm_id(base)

    sample_order = None
    treatment_label = None
    site_label = None
    patient_label = None
    inferred_patient_id = None

    if len(tokens) >= 5:
        sample_order = tokens[1]
        treatment_label = tokens[2]
        site_label = tokens[3]
        patient_label = tokens[4]
        inferred_patient_id = "_".join([treatment_label, site_label, patient_label])
    else:
        warnings.append("Filename does not contain enough tokens for Loret-style parsing.")

    if treatment_label not in {"N", "T", None}:
        warnings.append("Treatment token is not N or T.")

    if patient_label is not None and not re.match(r"PT\d+", patient_label):
        warnings.append("Patient token does not match PT<number> pattern.")

    return FilenameSampleParse(
        filename=filename,
        gsm_id=gsm_id,
        sample_order=sample_order,
        treatment_label=treatment_label,
        site_label=site_label,
        patient_label=patient_label,
        inferred_patient_id=inferred_patient_id,
        parser_label="loret_style_filename",
        warnings=warnings,
        evidence={"tokens": tokens},
    )


def parse_shen_style_filename(filename: str) -> FilenameSampleParse:
    base = strip_compression_suffix(Path(filename).name)
    tokens = base.split("_")

    warnings: list[str] = []
    gsm_id = extract_gsm_id(base)

    inferred_patient_id = None
    sample_order = None
    treatment_label = None
    site_label = None
    patient_label = None

    if len(tokens) >= 2:
        inferred_patient_id = tokens[1]
        match = re.match(r"(Pre|Post)-NACT(.+)", inferred_patient_id)
        if match:
            treatment_label = match.group(1)
            patient_label = match.group(2)
        else:
            warnings.append("Could not parse Pre/Post-NACT label from filename.")
    else:
        warnings.append("Filename does not contain enough tokens for Shen-style parsing.")

    return FilenameSampleParse(
        filename=filename,
        gsm_id=gsm_id,
        sample_order=sample_order,
        treatment_label=treatment_label,
        site_label=site_label,
        patient_label=patient_label,
        inferred_patient_id=inferred_patient_id,
        parser_label="shen_style_filename",
        warnings=warnings,
        evidence={"tokens": tokens},
    )


def parse_simple_gsm_prefixed_filename(filename: str) -> FilenameSampleParse:
    base = strip_compression_suffix(Path(filename).name)
    tokens = base.split("_")
    gsm_id = extract_gsm_id(base)

    sample_label = None
    if len(tokens) >= 2:
        sample_label = "_".join(tokens[1:])

    return FilenameSampleParse(
        filename=filename,
        gsm_id=gsm_id,
        sample_order=None,
        treatment_label=None,
        site_label=None,
        patient_label=None,
        inferred_patient_id=sample_label,
        parser_label="simple_gsm_prefixed_filename",
        warnings=[],
        evidence={"tokens": tokens},
    )


def parse_filename_sample_id(filename: str) -> FilenameSampleParse:
    base = strip_compression_suffix(Path(filename).name)
    tokens = base.split("_")

    if len(tokens) >= 5 and tokens[2] in {"N", "T"} and re.match(r"PT\d+", tokens[4]):
        return parse_loret_style_filename(filename)

    if len(tokens) >= 2 and re.match(r"(Pre|Post)-NACT", tokens[1]):
        return parse_shen_style_filename(filename)

    if extract_gsm_id(base) is not None:
        return parse_simple_gsm_prefixed_filename(filename)

    return FilenameSampleParse(
        filename=filename,
        gsm_id=None,
        sample_order=None,
        treatment_label=None,
        site_label=None,
        patient_label=None,
        inferred_patient_id=None,
        parser_label="unrecognized_filename_pattern",
        warnings=["Could not infer sample ID from filename."],
        evidence={"tokens": tokens},
    )


def parse_barcode_suffix(barcode: str) -> str | None:
    clean = barcode.strip().strip('"').strip("'")

    if "_" in clean:
        suffix = clean.rsplit("_", 1)[1]
        if suffix:
            return suffix.strip('"').strip("'")

    if "." in clean:
        suffix = clean.rsplit(".", 1)[1]
        if suffix.isdigit():
            return suffix

    return None


def summarize_barcode_suffixes(barcodes: list[str]) -> BarcodeSuffixParse:
    suffixes = []
    missing = 0

    for barcode in barcodes:
        suffix = parse_barcode_suffix(barcode)
        if suffix is None:
            missing += 1
        else:
            suffixes.append(suffix)

    counts = dict(sorted(Counter(suffixes).items(), key=lambda item: item[0]))

    warnings: list[str] = []
    if missing:
        warnings.append(f"{missing} barcodes did not contain a recognizable suffix.")

    if len(counts) == 0:
        parser_label = "no_suffix_detected"
    elif len(counts) == 1:
        parser_label = "single_suffix_detected"
    else:
        parser_label = "multiple_suffixes_detected"

    return BarcodeSuffixParse(
        total_barcodes=len(barcodes),
        suffix_counts=counts,
        parser_label=parser_label,
        warnings=warnings,
        evidence={"missing_suffix_count": missing},
    )


def validate_filename_against_metadata(
    parsed: FilenameSampleParse,
    metadata: dict[str, object],
) -> list[str]:
    warnings: list[str] = []

    curated_patient_id = metadata.get("patient_id")
    curated_treatment = metadata.get("tumor_treatment")
    curated_site = metadata.get("metastasis_site")

    if parsed.inferred_patient_id and curated_patient_id:
        if parsed.inferred_patient_id != curated_patient_id:
            warnings.append(
                f"Filename-derived patient_id '{parsed.inferred_patient_id}' does not match metadata patient_id '{curated_patient_id}'."
            )

    if parsed.treatment_label in {"N", "T"} and curated_treatment:
        expected = "No" if parsed.treatment_label == "N" else "Yes"
        if expected != curated_treatment:
            warnings.append(
                f"Filename treatment token '{parsed.treatment_label}' implies '{expected}' but metadata has '{curated_treatment}'."
            )

    site_map = {
        "A": "Ascites",
        "PER": "Peritoneum",
        "OM": "Omentum",
        "OT": "Omentum",
        "BL": "Bladder",
    }

    if parsed.site_label in site_map and curated_site:
        expected_site = site_map[parsed.site_label]
        if expected_site != curated_site:
            warnings.append(
                f"Filename site token '{parsed.site_label}' implies '{expected_site}' but metadata has '{curated_site}'."
            )

    return warnings


def print_filename_parse(parsed: FilenameSampleParse) -> None:
    print("=== Filename sample parse ===")
    print("Filename:", parsed.filename)
    print("Parser:", parsed.parser_label)
    print("GSM ID:", parsed.gsm_id)
    print("Sample order:", parsed.sample_order)
    print("Treatment label:", parsed.treatment_label)
    print("Site label:", parsed.site_label)
    print("Patient label:", parsed.patient_label)
    print("Inferred patient ID:", parsed.inferred_patient_id)

    if parsed.warnings:
        print("\nWarnings:")
        for warning in parsed.warnings:
            print("-", warning)

    print("\nEvidence:")
    for key, value in parsed.evidence.items():
        print(f"- {key}: {value}")


if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Parse sample identity from filenames or barcodes.")
    parser.add_argument("--filename", help="Filename to parse.")
    parser.add_argument("--barcodes", nargs="*", help="Barcodes to parse for suffixes.")

    args = parser.parse_args()

    if args.filename:
        print_filename_parse(parse_filename_sample_id(args.filename))

    if args.barcodes:
        result = summarize_barcode_suffixes(args.barcodes)
        print("\n=== Barcode suffix parse ===")
        print("Total barcodes:", result.total_barcodes)
        print("Parser:", result.parser_label)
        print("Suffix counts:", result.suffix_counts)
        if result.warnings:
            print("Warnings:")
            for warning in result.warnings:
                print("-", warning)
