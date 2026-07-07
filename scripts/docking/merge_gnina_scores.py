#!/usr/bin/env python3
"""Merge GNINA pose SDF files into ranked docking score tables."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

from rdkit import Chem, RDLogger


RDLogger.DisableLog("rdApp.*")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--poses-dir", required=True, type=Path)
    parser.add_argument("--candidate-csv", required=True, type=Path)
    parser.add_argument("--ligand-manifest-csv", required=True, type=Path)
    parser.add_argument("--output-csv", required=True, type=Path)
    parser.add_argument("--top-hits-csv", required=True, type=Path)
    parser.add_argument("--summary-json", required=True, type=Path)
    parser.add_argument("--top-n", default=500, type=int)
    parser.add_argument("--clinical-sim-max", default=0.45, type=float)
    parser.add_argument("--patent-sim-max", default=0.50, type=float)
    return parser.parse_args()


def read_csv_by_id(path: Path) -> dict[str, dict[str, str]]:
    with path.open() as handle:
        return {row["candidate_id"]: row for row in csv.DictReader(handle)}


def prop(mol: Chem.Mol, key: str, default: str = "") -> str:
    return mol.GetProp(key) if mol.HasProp(key) else default


def prop_float(mol: Chem.Mol, key: str) -> float | None:
    if not mol.HasProp(key):
        return None
    try:
        return float(mol.GetProp(key))
    except ValueError:
        return None


def row_float(row: dict[str, str], key: str) -> float | None:
    try:
        return float(row.get(key, ""))
    except ValueError:
        return None


def main() -> int:
    args = parse_args()
    args.output_csv.parent.mkdir(parents=True, exist_ok=True)
    args.top_hits_csv.parent.mkdir(parents=True, exist_ok=True)
    args.summary_json.parent.mkdir(parents=True, exist_ok=True)

    candidate_rows = read_csv_by_id(args.candidate_csv)
    manifest_rows = read_csv_by_id(args.ligand_manifest_csv)
    pose_files = sorted(args.poses_dir.glob("*_gnina_poses.sdf"))

    records: list[dict[str, object]] = []
    unreadable_mols = 0
    duplicate_ids: list[str] = []
    seen_ids: set[str] = set()

    for pose_file in pose_files:
        chunk_id = pose_file.name.replace("fxia_filtered_all_chunk_", "").split("_", 1)[0]
        for mol in Chem.SDMolSupplier(str(pose_file), removeHs=False, sanitize=False):
            if mol is None:
                unreadable_mols += 1
                continue
            candidate_id = prop(mol, "candidate_id") or prop(mol, "_Name")
            if candidate_id in seen_ids:
                duplicate_ids.append(candidate_id)
                continue
            seen_ids.add(candidate_id)

            candidate = candidate_rows.get(candidate_id, {})
            manifest = manifest_rows.get(candidate_id, {})
            rec: dict[str, object] = {
                "candidate_id": candidate_id,
                "chunk_id": chunk_id,
                "source_models": candidate.get("source_models") or prop(mol, "source_models"),
                "minimized_affinity_kcal_mol": prop_float(mol, "minimizedAffinity"),
                "cnn_score": prop_float(mol, "CNNscore"),
                "cnn_affinity": prop_float(mol, "CNNaffinity"),
                "cnn_vs": prop_float(mol, "CNN_VS"),
                "mw": row_float(candidate, "mw") if candidate else prop_float(mol, "mw"),
                "tpsa": row_float(candidate, "tpsa") if candidate else prop_float(mol, "tpsa"),
                "clogp": row_float(candidate, "clogp") if candidate else prop_float(mol, "clogp"),
                "hbd": row_float(candidate, "hbd") if candidate else prop_float(mol, "hbd"),
                "hba": row_float(candidate, "hba") if candidate else prop_float(mol, "hba"),
                "rotatable_bonds": row_float(candidate, "rotatable_bonds")
                if candidate
                else prop_float(mol, "rotatable_bonds"),
                "formal_charge": row_float(candidate, "formal_charge"),
                "sa_score": row_float(candidate, "sa_score"),
                "lipinski_rule_of_five_violations": row_float(candidate, "lipinski_rule_of_five_violations"),
                "max_clinical_similarity": row_float(candidate, "max_clinical_similarity")
                if candidate
                else prop_float(mol, "max_clinical_similarity"),
                "max_patent_proxy_similarity": row_float(candidate, "max_patent_proxy_similarity")
                if candidate
                else prop_float(mol, "max_patent_proxy_similarity"),
                "parent_inchi_key": candidate.get("parent_inchi_key") or manifest.get("parent_inchi_key") or prop(mol, "parent_inchi_key"),
                "parent_smiles": candidate.get("parent_smiles") or manifest.get("parent_smiles", ""),
                "canonical_smiles": candidate.get("canonical_smiles", ""),
                "generated_smiles": candidate.get("generated_smiles", ""),
                "pose_sdf": str(pose_file),
            }
            records.append(rec)

    def sort_key(row: dict[str, object]) -> tuple[float, float, str]:
        cnn_score = row["cnn_score"] if isinstance(row["cnn_score"], float) else -999.0
        affinity = row["minimized_affinity_kcal_mol"] if isinstance(row["minimized_affinity_kcal_mol"], float) else 999.0
        return (-cnn_score, affinity, str(row["candidate_id"]))

    records.sort(key=sort_key)
    for rank, row in enumerate(records, start=1):
        row["docking_rank"] = rank

    top_hits = [
        row
        for row in records
        if isinstance(row["max_clinical_similarity"], float)
        and isinstance(row["max_patent_proxy_similarity"], float)
        and row["max_clinical_similarity"] < args.clinical_sim_max
        and row["max_patent_proxy_similarity"] < args.patent_sim_max
    ][: args.top_n]

    fieldnames = [
        "docking_rank",
        "candidate_id",
        "chunk_id",
        "source_models",
        "minimized_affinity_kcal_mol",
        "cnn_score",
        "cnn_affinity",
        "cnn_vs",
        "mw",
        "tpsa",
        "clogp",
        "hbd",
        "hba",
        "rotatable_bonds",
        "formal_charge",
        "sa_score",
        "lipinski_rule_of_five_violations",
        "max_clinical_similarity",
        "max_patent_proxy_similarity",
        "parent_inchi_key",
        "parent_smiles",
        "canonical_smiles",
        "generated_smiles",
        "pose_sdf",
    ]

    for path, rows in ((args.output_csv, records), (args.top_hits_csv, top_hits)):
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames, lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)

    manifest_ids = set(manifest_rows)
    affinity_values = [
        row["minimized_affinity_kcal_mol"]
        for row in records
        if isinstance(row["minimized_affinity_kcal_mol"], float)
    ]
    summary = {
        "pose_files": len(pose_files),
        "scored_records": len(records),
        "manifest_records": len(manifest_rows),
        "missing_manifest_candidate_ids": sorted(manifest_ids - seen_ids)[:100],
        "missing_manifest_candidate_count": len(manifest_ids - seen_ids),
        "unreadable_molecules": unreadable_mols,
        "duplicate_candidate_count": len(duplicate_ids),
        "duplicate_candidate_examples": duplicate_ids[:20],
        "sort_order": "cnn_score desc, minimized_affinity_kcal_mol asc",
        "top_hits_count": len(top_hits),
        "top_hits_similarity_filter": {
            "max_clinical_similarity_lt": args.clinical_sim_max,
            "max_patent_proxy_similarity_lt": args.patent_sim_max,
        },
        "best_cnn_score": records[0]["cnn_score"] if records else None,
        "best_minimized_affinity_kcal_mol": min(affinity_values) if affinity_values else None,
        "top10_all_scored": records[:10],
        "top10_top_hits": top_hits[:10],
        "output_csv": str(args.output_csv),
        "top_hits_csv": str(args.top_hits_csv),
    }
    args.summary_json.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
