#!/usr/bin/env python3
"""Cluster GNINA top hits, flag simple pose issues, and export ADMET-AI input."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from functools import lru_cache
from pathlib import Path

from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold


RDLogger.DisableLog("rdApp.*")


OUTPUT_PREFIX = "top_hits"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-csv", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--admet-n", default=0, type=int, help="Rows for ADMET-AI input; <=0 keeps all valid hits.")
    parser.add_argument("--pre-admet-md-n", default=0, type=int, help="Optional exploratory pre-ADMET MD subset.")
    parser.add_argument("--fingerprint-sim-threshold", default=0.65, type=float)
    parser.add_argument("--max-per-scaffold-admet", default=5, type=int)
    parser.add_argument("--max-per-scaffold-md", default=2, type=int)
    parser.add_argument("--max-per-cluster-admet", default=3, type=int)
    parser.add_argument("--max-per-cluster-md", default=1, type=int)
    parser.add_argument("--write-sdf", action="store_true")
    parser.add_argument("--summary-json", type=Path)
    return parser.parse_args()


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def float_or_none(value: object) -> float | None:
    try:
        text = str(value).strip()
        if not text:
            return None
        return float(text)
    except (TypeError, ValueError):
        return None


def int_or_none(value: object) -> int | None:
    try:
        text = str(value).strip()
        if not text:
            return None
        return int(float(text))
    except (TypeError, ValueError):
        return None


def analysis_smiles(row: dict[str, str]) -> str:
    for field in ("parent_smiles", "canonical_smiles", "generated_smiles"):
        value = row.get(field, "").strip()
        if value:
            return value
    return ""


def canonicalize_smiles(smiles: str) -> tuple[Chem.Mol | None, str]:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, ""
    return mol, Chem.MolToSmiles(mol, isomericSmiles=True)


def scaffold_smiles(mol: Chem.Mol) -> tuple[str, str]:
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold.GetNumAtoms() == 0:
        canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
        return f"acyclic:{canonical}", "acyclic"

    murcko = Chem.MolToSmiles(scaffold, isomericSmiles=True)
    generic_mol = MurckoScaffold.MakeScaffoldGeneric(scaffold)
    generic = Chem.MolToSmiles(generic_mol, isomericSmiles=False) if generic_mol.GetNumAtoms() else ""
    return murcko, generic


def sort_key(row: dict[str, str]) -> tuple[float, float, int, str]:
    cnn_score = float_or_none(row.get("cnn_score"))
    affinity = float_or_none(row.get("minimized_affinity_kcal_mol"))
    rank = int_or_none(row.get("docking_rank"))
    return (
        -(cnn_score if cnn_score is not None else -999.0),
        affinity if affinity is not None else 999.0,
        rank if rank is not None else 999_999_999,
        row.get("candidate_id", ""),
    )


@lru_cache(maxsize=256)
def pose_mols_by_id(pose_sdf: str) -> dict[str, Chem.Mol]:
    path = Path(pose_sdf)
    if not path.exists():
        return {}

    mols: dict[str, Chem.Mol] = {}
    for mol in Chem.SDMolSupplier(str(path), removeHs=False, sanitize=False):
        if mol is None:
            continue
        candidate_id = mol.GetProp("candidate_id") if mol.HasProp("candidate_id") else mol.GetProp("_Name")
        if candidate_id:
            mols[candidate_id] = mol
    return mols


def min_nonbonded_heavy_distance(mol: Chem.Mol) -> float | None:
    if not mol.GetNumConformers():
        return None
    conf = mol.GetConformer()
    bonded = {tuple(sorted((bond.GetBeginAtomIdx(), bond.GetEndAtomIdx()))) for bond in mol.GetBonds()}
    heavy_atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1]
    min_distance: float | None = None

    for left_idx, atom_i in enumerate(heavy_atoms):
        pos_i = conf.GetAtomPosition(atom_i)
        for atom_j in heavy_atoms[left_idx + 1 :]:
            if tuple(sorted((atom_i, atom_j))) in bonded:
                continue
            pos_j = conf.GetAtomPosition(atom_j)
            distance = pos_i.Distance(pos_j)
            min_distance = distance if min_distance is None else min(min_distance, distance)
    return min_distance


def pose_qc(row: dict[str, str]) -> tuple[bool, list[str]]:
    flags: list[str] = []
    pose_sdf = row.get("pose_sdf", "").strip()
    candidate_id = row.get("candidate_id", "").strip()
    if not pose_sdf:
        return True, ["pose_not_checked"]

    pose_path = Path(pose_sdf)
    if not pose_path.exists():
        return False, ["pose_sdf_missing"]

    mol = pose_mols_by_id(pose_sdf).get(candidate_id)
    if mol is None:
        return False, ["candidate_pose_missing"]
    if not mol.GetNumConformers():
        return False, ["pose_has_no_conformer"]

    conf = mol.GetConformer()
    xs: list[float] = []
    ys: list[float] = []
    zs: list[float] = []
    for idx in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(idx)
        coords = (pos.x, pos.y, pos.z)
        if not all(math.isfinite(value) for value in coords):
            flags.append("pose_has_nonfinite_coordinate")
            continue
        xs.append(pos.x)
        ys.append(pos.y)
        zs.append(pos.z)

    if not xs:
        flags.append("pose_has_no_finite_coordinates")
    else:
        max_span = max(max(xs) - min(xs), max(ys) - min(ys), max(zs) - min(zs))
        if max_span < 0.1:
            flags.append("pose_collapsed_coordinates")
        if max_span > 35.0:
            flags.append("pose_span_gt_35a")

    min_distance = min_nonbonded_heavy_distance(mol)
    if min_distance is not None and min_distance < 0.70:
        flags.append("intramolecular_heavy_atom_clash")

    fail_flags = {
        "pose_has_nonfinite_coordinate",
        "pose_has_no_finite_coordinates",
        "pose_collapsed_coordinates",
        "pose_span_gt_35a",
        "intramolecular_heavy_atom_clash",
    }
    return not any(flag in fail_flags for flag in flags), flags or ["ok"]


def assign_ecfp_clusters(
    rows: list[dict[str, str]],
    fingerprints: dict[str, DataStructs.cDataStructs.ExplicitBitVect],
    similarity_threshold: float,
) -> dict[str, int]:
    cluster_centers: list[tuple[int, DataStructs.cDataStructs.ExplicitBitVect]] = []
    assignments: dict[str, int] = {}

    for row in rows:
        candidate_id = row["candidate_id"]
        fp = fingerprints.get(candidate_id)
        if fp is None:
            assignments[candidate_id] = -1
            continue

        assigned_cluster: int | None = None
        for cluster_id, center_fp in cluster_centers:
            if DataStructs.TanimotoSimilarity(fp, center_fp) >= similarity_threshold:
                assigned_cluster = cluster_id
                break

        if assigned_cluster is None:
            assigned_cluster = len(cluster_centers) + 1
            cluster_centers.append((assigned_cluster, fp))
        assignments[candidate_id] = assigned_cluster

    return assignments


def choose_balanced_subset(
    rows: list[dict[str, str]],
    target_n: int,
    max_per_scaffold: int,
    max_per_cluster: int,
) -> list[dict[str, str]]:
    selected: list[dict[str, str]] = []
    scaffold_counts: Counter[str] = Counter()
    cluster_counts: Counter[str] = Counter()

    for row in rows:
        if row.get("pose_qc_pass") != "true":
            continue
        scaffold = row.get("murcko_scaffold", "")
        cluster = row.get("ecfp_cluster_id", "-1")
        if scaffold_counts[scaffold] >= max_per_scaffold:
            continue
        if cluster_counts[cluster] >= max_per_cluster:
            continue
        selected.append(row)
        scaffold_counts[scaffold] += 1
        cluster_counts[cluster] += 1
        if len(selected) >= target_n:
            break

    if len(selected) >= target_n:
        return selected

    selected_ids = {row["candidate_id"] for row in selected}
    for row in rows:
        if row["candidate_id"] in selected_ids or row.get("pose_qc_pass") != "true":
            continue
        scaffold = row.get("murcko_scaffold", "")
        if scaffold_counts[scaffold] >= max_per_scaffold:
            continue
        selected.append(row)
        scaffold_counts[scaffold] += 1
        selected_ids.add(row["candidate_id"])
        if len(selected) >= target_n:
            break

    if len(selected) >= target_n:
        return selected

    for row in rows:
        if row["candidate_id"] in selected_ids or row.get("pose_qc_pass") != "true":
            continue
        selected.append(row)
        selected_ids.add(row["candidate_id"])
        if len(selected) >= target_n:
            break
    return selected


def choose_admet_input_rows(
    rows: list[dict[str, str]],
    target_n: int,
    max_per_scaffold: int,
    max_per_cluster: int,
) -> list[dict[str, str]]:
    valid_rows = [
        row
        for row in rows
        if row.get("structure_qc_flags") == "ok" and row.get("analysis_smiles", "").strip()
    ]
    if target_n <= 0 or target_n >= len(valid_rows):
        return valid_rows
    return choose_balanced_subset(
        valid_rows,
        target_n=target_n,
        max_per_scaffold=max_per_scaffold,
        max_per_cluster=max_per_cluster,
    )


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def write_admet_ai_input(path: Path, rows: list[dict[str, str]]) -> None:
    fieldnames = [
        "candidate_id",
        "smiles",
        "docking_rank",
        "source_models",
        "cnn_score",
        "cnn_affinity",
        "minimized_affinity_kcal_mol",
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
        "murcko_scaffold",
        "generic_murcko_scaffold",
        "ecfp_cluster_id",
        "pose_qc_pass",
        "pose_qc_flags",
        "parent_inchi_key",
        "parent_smiles",
        "pose_sdf",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            output_row = dict(row)
            output_row["smiles"] = row.get("analysis_smiles", "")
            writer.writerow(output_row)


def write_selected_sdf(path: Path, rows: list[dict[str, str]]) -> int:
    path.parent.mkdir(parents=True, exist_ok=True)
    written = 0
    with Chem.SDWriter(str(path)) as writer:
        for row in rows:
            pose_sdf = row.get("pose_sdf", "").strip()
            candidate_id = row.get("candidate_id", "").strip()
            if not pose_sdf:
                continue
            mol = pose_mols_by_id(pose_sdf).get(candidate_id)
            if mol is None:
                continue
            for key, value in row.items():
                mol.SetProp(key, str(value))
            writer.write(mol)
            written += 1
    return written


def annotate_rows(rows: list[dict[str, str]], fingerprint_sim_threshold: float) -> list[dict[str, str]]:
    sorted_rows = sorted(rows, key=sort_key)
    generator = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fingerprints: dict[str, DataStructs.cDataStructs.ExplicitBitVect] = {}

    for row in sorted_rows:
        smiles = analysis_smiles(row)
        mol, canonical = canonicalize_smiles(smiles)
        row["analysis_smiles"] = canonical or smiles
        if mol is None:
            row["murcko_scaffold"] = ""
            row["generic_murcko_scaffold"] = ""
            row["structure_qc_flags"] = "invalid_smiles"
            row["pose_qc_pass"] = "false"
            row["pose_qc_flags"] = "invalid_smiles"
            continue

        murcko, generic = scaffold_smiles(mol)
        row["murcko_scaffold"] = murcko
        row["generic_murcko_scaffold"] = generic
        row["structure_qc_flags"] = "ok"
        fingerprints[row["candidate_id"]] = generator.GetFingerprint(mol)

        pose_pass, pose_flags = pose_qc(row)
        row["pose_qc_pass"] = "true" if pose_pass else "false"
        row["pose_qc_flags"] = ";".join(pose_flags)

    cluster_assignments = assign_ecfp_clusters(sorted_rows, fingerprints, fingerprint_sim_threshold)
    for row in sorted_rows:
        row["ecfp_cluster_id"] = str(cluster_assignments.get(row["candidate_id"], -1))

    return sorted_rows


def add_selection_flags(
    annotated_rows: list[dict[str, str]],
    admet_rows: list[dict[str, str]],
    pre_admet_md_rows: list[dict[str, str]],
) -> None:
    admet_ids = {row["candidate_id"] for row in admet_rows}
    md_ids = {row["candidate_id"] for row in pre_admet_md_rows}
    admet_rank = {row["candidate_id"]: str(idx) for idx, row in enumerate(admet_rows, start=1)}
    md_rank = {row["candidate_id"]: str(idx) for idx, row in enumerate(pre_admet_md_rows, start=1)}

    for row in annotated_rows:
        candidate_id = row["candidate_id"]
        row["selected_for_admet"] = "true" if candidate_id in admet_ids else "false"
        row["selected_for_pre_admet_md"] = "true" if candidate_id in md_ids else "false"
        row["admet_selection_rank"] = admet_rank.get(candidate_id, "")
        row["pre_admet_md_selection_rank"] = md_rank.get(candidate_id, "")


def summarize(
    rows: list[dict[str, str]],
    admet_rows: list[dict[str, str]],
    pre_admet_md_rows: list[dict[str, str]],
) -> dict[str, object]:
    qc_counts = Counter(row.get("pose_qc_pass", "false") for row in rows)
    flag_counts: Counter[str] = Counter()
    for row in rows:
        for flag in row.get("pose_qc_flags", "").split(";"):
            if flag:
                flag_counts[flag] += 1

    source_counts: dict[str, int] = defaultdict(int)
    for row in pre_admet_md_rows:
        source_counts[row.get("source_models", "")] += 1

    return {
        "input_rows": len(rows),
        "pose_qc_pass_count": qc_counts.get("true", 0),
        "pose_qc_fail_count": qc_counts.get("false", 0),
        "pose_qc_flag_counts": dict(sorted(flag_counts.items())),
        "unique_murcko_scaffolds": len({row.get("murcko_scaffold", "") for row in rows if row.get("murcko_scaffold")}),
        "unique_ecfp_clusters": len({row.get("ecfp_cluster_id", "") for row in rows if row.get("ecfp_cluster_id") not in {"", "-1"}}),
        "admet_input_count": len(admet_rows),
        "admet_input_unique_murcko_scaffolds": len({row.get("murcko_scaffold", "") for row in admet_rows}),
        "admet_subset_count": len(admet_rows),
        "admet_unique_murcko_scaffolds": len({row.get("murcko_scaffold", "") for row in admet_rows}),
        "pre_admet_md_subset_count": len(pre_admet_md_rows),
        "pre_admet_md_unique_murcko_scaffolds": len({row.get("murcko_scaffold", "") for row in pre_admet_md_rows}),
        "pre_admet_md_source_model_counts": dict(sorted(source_counts.items())),
        "pre_admet_md_top_candidates": [
            {
                "candidate_id": row.get("candidate_id"),
                "docking_rank": row.get("docking_rank"),
                "cnn_score": row.get("cnn_score"),
                "minimized_affinity_kcal_mol": row.get("minimized_affinity_kcal_mol"),
                "murcko_scaffold": row.get("murcko_scaffold"),
                "source_models": row.get("source_models"),
            }
            for row in pre_admet_md_rows[:10]
        ],
    }


def select_docking_hits(
    input_csv: Path,
    output_dir: Path,
    admet_n: int = 0,
    pre_admet_md_n: int = 0,
    fingerprint_sim_threshold: float = 0.65,
    max_per_scaffold_admet: int = 5,
    max_per_scaffold_md: int = 2,
    max_per_cluster_admet: int = 3,
    max_per_cluster_md: int = 1,
    write_sdf: bool = False,
    summary_json: Path | None = None,
) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    rows = read_csv(input_csv)
    annotated_rows = annotate_rows(rows, fingerprint_sim_threshold)
    admet_rows = choose_admet_input_rows(
        annotated_rows,
        target_n=admet_n,
        max_per_scaffold=max_per_scaffold_admet,
        max_per_cluster=max_per_cluster_admet,
    )
    pre_admet_md_rows = (
        choose_balanced_subset(
            annotated_rows,
            target_n=pre_admet_md_n,
            max_per_scaffold=max_per_scaffold_md,
            max_per_cluster=max_per_cluster_md,
        )
        if pre_admet_md_n > 0
        else []
    )
    add_selection_flags(annotated_rows, admet_rows, pre_admet_md_rows)

    extra_fields = [
        "analysis_smiles",
        "murcko_scaffold",
        "generic_murcko_scaffold",
        "ecfp_cluster_id",
        "structure_qc_flags",
        "pose_qc_pass",
        "pose_qc_flags",
        "selected_for_admet",
        "selected_for_pre_admet_md",
        "admet_selection_rank",
        "pre_admet_md_selection_rank",
    ]
    fieldnames = list(dict.fromkeys([*rows[0].keys(), *extra_fields])) if rows else extra_fields
    clustered_csv = output_dir / f"{OUTPUT_PREFIX}_clustered.csv"
    admet_ai_csv = output_dir / f"{OUTPUT_PREFIX}_admet_ai_input.csv"
    write_csv(clustered_csv, annotated_rows, fieldnames)
    write_admet_ai_input(admet_ai_csv, admet_rows)
    pre_admet_md_csv = output_dir / f"{OUTPUT_PREFIX}_pre_admet_md_candidates.csv"
    if pre_admet_md_rows:
        write_csv(pre_admet_md_csv, pre_admet_md_rows, fieldnames)

    summary = summarize(annotated_rows, admet_rows, pre_admet_md_rows)
    summary.update(
        {
            "input_csv": str(input_csv),
            "clustered_csv": str(clustered_csv),
            "admet_ai_input_csv": str(admet_ai_csv),
            "selection_rules": {
                "fingerprint": "ECFP4/Morgan radius 2, 2048 bits",
                "fingerprint_sim_threshold": fingerprint_sim_threshold,
                "admet_n": admet_n,
                "pre_admet_md_n": pre_admet_md_n,
                "max_per_scaffold_admet": max_per_scaffold_admet,
                "max_per_scaffold_md": max_per_scaffold_md,
                "max_per_cluster_admet": max_per_cluster_admet,
                "max_per_cluster_md": max_per_cluster_md,
            },
        }
    )
    if pre_admet_md_rows:
        summary["pre_admet_md_candidates_csv"] = str(pre_admet_md_csv)

    if write_sdf and pre_admet_md_rows:
        pre_admet_md_sdf = output_dir / f"{OUTPUT_PREFIX}_pre_admet_md_candidates.sdf"
        summary["pre_admet_md_candidates_sdf"] = str(pre_admet_md_sdf)
        summary["pre_admet_md_candidates_sdf_molecules"] = write_selected_sdf(pre_admet_md_sdf, pre_admet_md_rows)

    summary_path = summary_json or output_dir / f"{OUTPUT_PREFIX}_selection_summary.json"
    summary["summary_json"] = str(summary_path)
    summary_path.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    return summary


def main() -> int:
    args = parse_args()
    summary = select_docking_hits(
        input_csv=args.input_csv,
        output_dir=args.output_dir,
        admet_n=args.admet_n,
        pre_admet_md_n=args.pre_admet_md_n,
        fingerprint_sim_threshold=args.fingerprint_sim_threshold,
        max_per_scaffold_admet=args.max_per_scaffold_admet,
        max_per_scaffold_md=args.max_per_scaffold_md,
        max_per_cluster_admet=args.max_per_cluster_admet,
        max_per_cluster_md=args.max_per_cluster_md,
        write_sdf=args.write_sdf,
        summary_json=args.summary_json,
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
