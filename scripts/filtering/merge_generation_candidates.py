#!/usr/bin/env python
"""Merge filtered generated molecules into a de-duplicated candidate pool."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


MERGED_COLUMNS = [
    "candidate_id",
    "source_models",
    "source_row_count",
    "generated_smiles",
    "canonical_smiles",
    "parent_smiles",
    "parent_inchi_key",
    "mw",
    "tpsa",
    "clogp",
    "hbd",
    "hba",
    "rotatable_bonds",
    "formal_charge",
    "sa_score",
    "lipinski_rule_of_five_violations",
    "passes_lipinski_rule_of_five",
    "max_clinical_similarity",
    "max_patent_proxy_similarity",
]


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def row_sort_key(row: dict[str, str]) -> tuple[float, float, float]:
    clinical = float(row.get("max_clinical_similarity") or 1.0)
    patent = float(row.get("max_patent_proxy_similarity") or 1.0)
    sa_score = float(row.get("sa_score") or 99.0)
    return clinical, patent, sa_score


def merge_candidate_files(inputs: list[tuple[str, Path]], output_csv: str | Path) -> dict[str, object]:
    output_csv = Path(output_csv)
    grouped: dict[str, list[tuple[str, dict[str, str]]]] = {}
    input_rows = 0
    passing_rows = 0

    for source_model, path in inputs:
        for row in read_csv(Path(path)):
            input_rows += 1
            if row.get("passes_all_filters") != "true":
                continue
            key = row.get("parent_inchi_key") or row.get("parent_smiles") or row.get("canonical_smiles")
            if not key:
                continue
            passing_rows += 1
            grouped.setdefault(key, []).append((source_model, row))

    merged_rows: list[dict[str, object]] = []
    for idx, (_key, grouped_rows) in enumerate(sorted(grouped.items()), start=1):
        best_source, best_row = sorted(grouped_rows, key=lambda item: row_sort_key(item[1]))[0]
        source_models = sorted({source for source, _row in grouped_rows})
        merged_rows.append(
            {
                "candidate_id": f"FXIA-GEN-{idx:06d}",
                "source_models": ";".join(source_models),
                "source_row_count": len(grouped_rows),
                **best_row,
            }
        )

    write_csv(output_csv, merged_rows, MERGED_COLUMNS)
    summary = {
        "input_rows": input_rows,
        "passing_rows": passing_rows,
        "merged_rows": len(merged_rows),
        "output_csv": str(output_csv.resolve()),
        "inputs": [{"source_model": model, "path": str(Path(path).resolve())} for model, path in inputs],
    }
    output_csv.with_suffix(".summary.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    return summary


def parse_input(value: str) -> tuple[str, Path]:
    if "=" not in value:
        raise argparse.ArgumentTypeError("inputs must use source_model=path")
    model, path = value.split("=", 1)
    return model, Path(path)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", action="append", required=True, type=parse_input, dest="inputs")
    parser.add_argument("--output", required=True, dest="output_csv")
    args = parser.parse_args()
    summary = merge_candidate_files(args.inputs, args.output_csv)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
