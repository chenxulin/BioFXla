#!/usr/bin/env python3
"""Rank docked FXIa hits after ADMET-AI prediction and export MD candidates."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from collections import Counter
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from scripts.docking.select_docking_hits import choose_balanced_subset, float_or_none, write_csv, write_selected_sdf


OUTPUT_PREFIX = "top_hits"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--docking-csv", required=True, type=Path)
    parser.add_argument("--admet-csv", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--top-n", default=50, type=int)
    parser.add_argument("--max-per-scaffold", default=2, type=int)
    parser.add_argument("--max-per-cluster", default=1, type=int)
    parser.add_argument("--write-sdf", action="store_true")
    parser.add_argument("--summary-json", type=Path)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def read_csv_by_id(path: Path) -> dict[str, dict[str, str]]:
    rows = read_csv(path)
    return {row["candidate_id"]: row for row in rows if row.get("candidate_id")}


def clamp01(value: float) -> float:
    return max(0.0, min(1.0, value))


def high_good(value: float | None, low: float, high: float, default: float = 0.5) -> float:
    if value is None:
        return default
    if high == low:
        return default
    return clamp01((value - low) / (high - low))


def low_good(value: float | None, good: float, bad: float, default: float = 0.5) -> float:
    if value is None:
        return default
    if bad == good:
        return default
    return clamp01((bad - value) / (bad - good))


def field_float(row: dict[str, str], *names: str) -> float | None:
    for name in names:
        value = float_or_none(row.get(name))
        if value is not None:
            return value
    return None


def count_risk_flags(row: dict[str, str]) -> tuple[int, list[str]]:
    flags: list[str] = []
    for field in ["hERG", "AMES", "DILI", "ClinTox", "Carcinogens_Lagunin"]:
        value = field_float(row, field)
        if value is not None and value >= 0.50:
            flags.append(f"{field}_ge_0.50")

    for field in ["CYP1A2_Veith", "CYP2C19_Veith", "CYP2C9_Veith", "CYP2D6_Veith", "CYP3A4_Veith"]:
        value = field_float(row, field)
        if value is not None and value >= 0.70:
            flags.append(f"{field}_ge_0.70")

    for field in ["PAINS_alert", "BRENK_alert", "NIH_alert"]:
        value = field_float(row, field)
        if value is not None and value > 0:
            flags.append(f"{field}_present")

    pgp = field_float(row, "Pgp_Broccatelli")
    if pgp is not None and pgp >= 0.70:
        flags.append("Pgp_Broccatelli_ge_0.70")

    return len(flags), flags


def docking_component(row: dict[str, str]) -> float:
    cnn = high_good(field_float(row, "cnn_score"), 0.50, 1.00)
    affinity = low_good(field_float(row, "minimized_affinity_kcal_mol"), -10.0, -5.0)
    return 0.70 * cnn + 0.30 * affinity


def physchem_component(row: dict[str, str]) -> float:
    mw = field_float(row, "molecular_weight", "mw")
    tpsa = field_float(row, "tpsa")
    logp = field_float(row, "logP", "clogp")
    hbd = field_float(row, "hydrogen_bond_donors", "hbd")
    hba = field_float(row, "hydrogen_bond_acceptors", "hba")
    lipinski = field_float(row, "Lipinski", "lipinski_rule_of_five_violations")
    qed = field_float(row, "QED")

    scores = [
        1.0 if mw is not None and 350.0 <= mw <= 650.0 else 0.5,
        1.0 if tpsa is not None and 60.0 <= tpsa <= 140.0 else 0.5,
        1.0 if logp is not None and 1.0 <= logp <= 5.0 else 0.5,
        1.0 if hbd is not None and hbd <= 3 else 0.5,
        1.0 if hba is not None and hba <= 12 else 0.5,
        low_good(lipinski, 0.0, 2.0),
        high_good(qed, 0.25, 0.85),
    ]
    return sum(scores) / len(scores)


def admet_component(row: dict[str, str]) -> float:
    absorption_scores = [
        high_good(field_float(row, "HIA_Hou"), 0.0, 1.0),
        high_good(field_float(row, "Bioavailability_Ma"), 0.0, 1.0),
        high_good(field_float(row, "PAMPA_NCATS"), 0.0, 1.0),
        high_good(field_float(row, "Solubility_AqSolDB"), -6.0, -2.0),
        high_good(field_float(row, "Caco2_Wang"), -6.0, -4.0),
    ]
    absorption = sum(absorption_scores) / len(absorption_scores)

    tox_fields = ["hERG", "AMES", "DILI", "ClinTox", "Carcinogens_Lagunin"]
    tox = sum(low_good(field_float(row, field), 0.0, 1.0) for field in tox_fields) / len(tox_fields)

    cyp_fields = ["CYP1A2_Veith", "CYP2C19_Veith", "CYP2C9_Veith", "CYP2D6_Veith", "CYP3A4_Veith"]
    cyp = sum(low_good(field_float(row, field), 0.0, 1.0) for field in cyp_fields) / len(cyp_fields)

    alert_penalty = 0.0
    for field in ["PAINS_alert", "BRENK_alert", "NIH_alert"]:
        value = field_float(row, field)
        if value is not None and value > 0:
            alert_penalty += 0.08

    return clamp01(0.45 * absorption + 0.40 * tox + 0.15 * cyp - alert_penalty)


def similarity_component(row: dict[str, str]) -> float:
    clinical = high_good(field_float(row, "max_clinical_similarity"), 0.0, 0.45)
    patent = high_good(field_float(row, "max_patent_proxy_similarity"), 0.0, 0.50)
    return clamp01(1.0 - (0.50 * clinical + 0.50 * patent))


def score_row(row: dict[str, str]) -> None:
    dock = docking_component(row)
    physchem = physchem_component(row)
    admet = admet_component(row)
    similarity = similarity_component(row)
    risk_count, risk_flags = count_risk_flags(row)
    risk_penalty = min(0.35, 0.08 * risk_count)
    composite = clamp01(0.35 * dock + 0.40 * admet + 0.15 * physchem + 0.10 * similarity - risk_penalty)

    row["docking_component_score"] = f"{dock:.6f}"
    row["admet_component_score"] = f"{admet:.6f}"
    row["physchem_component_score"] = f"{physchem:.6f}"
    row["similarity_component_score"] = f"{similarity:.6f}"
    row["admet_risk_count"] = str(risk_count)
    row["admet_risk_flags"] = ";".join(risk_flags) if risk_flags else "none"
    row["composite_score"] = f"{composite:.6f}"


def merge_rows(docking_rows: list[dict[str, str]], admet_rows_by_id: dict[str, dict[str, str]]) -> list[dict[str, str]]:
    merged_rows: list[dict[str, str]] = []
    for docking_row in docking_rows:
        candidate_id = docking_row.get("candidate_id", "")
        admet_row = admet_rows_by_id.get(candidate_id)
        if not admet_row:
            continue
        merged = dict(docking_row)
        for key, value in admet_row.items():
            if key == "candidate_id":
                continue
            if key == "smiles":
                merged["admet_ai_smiles"] = value
            elif key not in merged or merged[key] == "":
                merged[key] = value
            else:
                merged[f"admet_ai_{key}"] = value
        score_row(merged)
        merged_rows.append(merged)
    return merged_rows


def ranking_key(row: dict[str, str]) -> tuple[float, int, int, str]:
    composite = field_float(row, "composite_score") or 0.0
    risk_count = int(field_float(row, "admet_risk_count") or 0)
    docking_rank = int(field_float(row, "docking_rank") or 999_999_999)
    return (-composite, risk_count, docking_rank, row.get("candidate_id", ""))


def add_md_selection_flags(ranked_rows: list[dict[str, str]], selected_rows: list[dict[str, str]]) -> None:
    selected_ids = {row["candidate_id"] for row in selected_rows}
    selected_rank = {row["candidate_id"]: str(idx) for idx, row in enumerate(selected_rows, start=1)}
    for rank, row in enumerate(ranked_rows, start=1):
        row["admet_rank"] = str(rank)
        row["selected_for_md"] = "true" if row["candidate_id"] in selected_ids else "false"
        row["md_selection_rank"] = selected_rank.get(row["candidate_id"], "")


def summarize(
    docking_rows: list[dict[str, str]],
    admet_rows_by_id: dict[str, dict[str, str]],
    ranked_rows: list[dict[str, str]],
    selected_rows: list[dict[str, str]],
) -> dict[str, object]:
    risk_counts = Counter(row.get("admet_risk_count", "") for row in ranked_rows)
    source_counts = Counter(row.get("source_models", "") for row in selected_rows)
    missing_admet = sorted(
        row.get("candidate_id", "")
        for row in docking_rows
        if row.get("candidate_id", "") not in admet_rows_by_id
    )
    return {
        "docking_rows": len(docking_rows),
        "admet_rows": len(admet_rows_by_id),
        "ranked_rows": len(ranked_rows),
        "missing_admet_count": len(missing_admet),
        "missing_admet_candidate_ids": missing_admet[:100],
        "selected_count": len(selected_rows),
        "selected_unique_murcko_scaffolds": len({row.get("murcko_scaffold", "") for row in selected_rows}),
        "selected_unique_ecfp_clusters": len({row.get("ecfp_cluster_id", "") for row in selected_rows}),
        "selected_source_model_counts": dict(sorted(source_counts.items())),
        "admet_risk_count_distribution": dict(sorted(risk_counts.items())),
        "top10_selected": [
            {
                "candidate_id": row.get("candidate_id"),
                "admet_rank": row.get("admet_rank"),
                "docking_rank": row.get("docking_rank"),
                "composite_score": row.get("composite_score"),
                "cnn_score": row.get("cnn_score"),
                "admet_component_score": row.get("admet_component_score"),
                "admet_risk_count": row.get("admet_risk_count"),
                "admet_risk_flags": row.get("admet_risk_flags"),
                "murcko_scaffold": row.get("murcko_scaffold"),
                "source_models": row.get("source_models"),
            }
            for row in selected_rows[:10]
        ],
    }


def rank_admet_docking_hits(
    docking_csv: Path,
    admet_csv: Path,
    output_dir: Path,
    top_n: int = 50,
    max_per_scaffold: int = 2,
    max_per_cluster: int = 1,
    write_sdf: bool = False,
    summary_json: Path | None = None,
) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    docking_rows = read_csv(docking_csv)
    admet_rows_by_id = read_csv_by_id(admet_csv)
    ranked_rows = merge_rows(docking_rows, admet_rows_by_id)
    ranked_rows.sort(key=ranking_key)
    selected_rows = choose_balanced_subset(
        ranked_rows,
        target_n=top_n,
        max_per_scaffold=max_per_scaffold,
        max_per_cluster=max_per_cluster,
    )
    add_md_selection_flags(ranked_rows, selected_rows)

    extra_fields = [
        "admet_ai_smiles",
        "docking_component_score",
        "admet_component_score",
        "physchem_component_score",
        "similarity_component_score",
        "admet_risk_count",
        "admet_risk_flags",
        "composite_score",
        "admet_rank",
        "selected_for_md",
        "md_selection_rank",
    ]
    all_fields = list(ranked_rows[0].keys()) if ranked_rows else extra_fields
    fieldnames = list(dict.fromkeys([*all_fields, *extra_fields]))

    ranked_csv = output_dir / f"{OUTPUT_PREFIX}_admet_ranked.csv"
    md_csv = output_dir / f"{OUTPUT_PREFIX}_md_candidates.csv"
    write_csv(ranked_csv, ranked_rows, fieldnames)
    write_csv(md_csv, selected_rows, fieldnames)

    summary = summarize(docking_rows, admet_rows_by_id, ranked_rows, selected_rows)
    summary.update(
        {
            "docking_csv": str(docking_csv),
            "admet_csv": str(admet_csv),
            "ranked_csv": str(ranked_csv),
            "md_candidates_csv": str(md_csv),
            "selection_rules": {
                "top_n": top_n,
                "max_per_scaffold": max_per_scaffold,
                "max_per_cluster": max_per_cluster,
                "composite_score": "0.35*docking + 0.40*ADMET + 0.15*physchem + 0.10*similarity - ADMET risk penalty",
                "diversity": "greedy selection after composite ranking by Murcko scaffold and ECFP cluster caps",
            },
        }
    )

    if write_sdf:
        md_sdf = output_dir / f"{OUTPUT_PREFIX}_md_candidates.sdf"
        summary["md_candidates_sdf"] = str(md_sdf)
        summary["md_candidates_sdf_molecules"] = write_selected_sdf(md_sdf, selected_rows)

    summary_path = summary_json or output_dir / f"{OUTPUT_PREFIX}_admet_rank_summary.json"
    summary["summary_json"] = str(summary_path)
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    return summary


def main() -> int:
    args = parse_args()
    summary = rank_admet_docking_hits(
        docking_csv=args.docking_csv,
        admet_csv=args.admet_csv,
        output_dir=args.output_dir,
        top_n=args.top_n,
        max_per_scaffold=args.max_per_scaffold,
        max_per_cluster=args.max_per_cluster,
        write_sdf=args.write_sdf,
        summary_json=args.summary_json,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
