#!/usr/bin/env python
"""Filter and calibrate FXIa generated molecules with RDKit rules."""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from pathlib import Path
from statistics import median

from rdkit import Chem, DataStructs
from rdkit import RDLogger
from rdkit.Chem import AllChem, Crippen, Descriptors, Lipinski, rdMolDescriptors
from rdkit.Chem.MolStandardize import rdMolStandardize


ROOT = Path(__file__).resolve().parents[2]
RDLogger.DisableLog("rdApp.warning")
SASCORE_DIR = ROOT / "tools" / "REINVENT4" / "reinvent_plugins" / "components" / "SAScore"
sys.path.insert(0, str(SASCORE_DIR))
import sascorer  # noqa: E402


FILTER_COLUMNS = [
    "generated_smiles",
    "canonical_smiles",
    "parent_smiles",
    "parent_inchi_key",
    "input_smiles",
    "reinvent_tanimoto",
    "nll",
    "valid",
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
    "passes_property_filters",
    "passes_clinical_similarity_filter",
    "passes_patent_proxy_soft_filter",
    "passes_hard_filters",
    "passes_all_filters",
    "failure_reasons",
]

DEFAULT_THRESHOLDS = {
    "mw_min": 350.0,
    "mw_max": 650.0,
    "tpsa_min": 60.0,
    "tpsa_max": 140.0,
    "clogp_min": 1.0,
    "clogp_max": 5.0,
    "hbd_max": 3,
    "hba_max": 12,
    "rotatable_bonds_max": 10,
    "formal_charge_min": -1,
    "formal_charge_max": 1,
    "sa_score_max": 4.0,
    "lipinski_mw_max": 500.0,
    "lipinski_clogp_max": 5.0,
    "lipinski_hbd_max": 5,
    "lipinski_hba_max": 10,
    "clinical_similarity_max": 0.45,
    "patent_proxy_similarity_max": 0.50,
    "patent_proxy_soft_similarity_max": 0.85,
}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, rows: list[dict[str, object]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parent_mol(smiles: str) -> Chem.Mol | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    parent = rdMolStandardize.FragmentParent(mol)
    parent.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(parent)
    return parent


def smiles_fingerprint(smiles: str):
    mol = parent_mol(smiles)
    if mol is None:
        return None
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)


def load_exclusion_fps(path: Path) -> list:
    fps = []
    if not path.exists():
        return fps
    for line in path.read_text().splitlines():
        smiles = line.strip().split()[0] if line.strip() else ""
        if not smiles:
            continue
        fp = smiles_fingerprint(smiles)
        if fp is not None:
            fps.append(fp)
    return fps


def max_similarity(fp, exclusion_fps: list) -> float:
    if fp is None or not exclusion_fps:
        return 0.0
    return max(DataStructs.BulkTanimotoSimilarity(fp, exclusion_fps), default=0.0)


def mol_props(mol: Chem.Mol) -> dict[str, float | int | str]:
    parent = rdMolStandardize.FragmentParent(mol)
    parent.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(parent)
    canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
    parent_smiles = Chem.MolToSmiles(parent, isomericSmiles=True)
    mw = round(Descriptors.MolWt(parent), 3)
    clogp = round(Crippen.MolLogP(parent), 3)
    hbd = Lipinski.NumHDonors(parent)
    hba = Lipinski.NumHAcceptors(parent)
    lipinski_violations = sum(
        [
            mw > DEFAULT_THRESHOLDS["lipinski_mw_max"],
            clogp > DEFAULT_THRESHOLDS["lipinski_clogp_max"],
            hbd > DEFAULT_THRESHOLDS["lipinski_hbd_max"],
            hba > DEFAULT_THRESHOLDS["lipinski_hba_max"],
        ]
    )
    return {
        "canonical_smiles": canonical,
        "parent_smiles": parent_smiles,
        "parent_inchi_key": Chem.MolToInchiKey(parent),
        "mw": mw,
        "tpsa": round(rdMolDescriptors.CalcTPSA(parent), 3),
        "clogp": clogp,
        "hbd": hbd,
        "hba": hba,
        "rotatable_bonds": Lipinski.NumRotatableBonds(parent),
        "formal_charge": sum(atom.GetFormalCharge() for atom in parent.GetAtoms()),
        "sa_score": round(float(sascorer.calculateScore(parent)), 3),
        "lipinski_rule_of_five_violations": int(lipinski_violations),
        "passes_lipinski_rule_of_five": str(lipinski_violations == 0).lower(),
    }


def parse_float(value: str) -> float | str:
    try:
        return float(value)
    except (TypeError, ValueError):
        return ""


def assess_filters(props: dict[str, object], clinical_sim: float, patent_sim: float, thresholds: dict[str, float]) -> tuple[bool, bool, bool, list[str]]:
    failures: list[str] = []
    checks = [
        ("mw_low", float(props["mw"]) >= thresholds["mw_min"]),
        ("mw_high", float(props["mw"]) <= thresholds["mw_max"]),
        ("tpsa_low", float(props["tpsa"]) >= thresholds["tpsa_min"]),
        ("tpsa_high", float(props["tpsa"]) <= thresholds["tpsa_max"]),
        ("clogp_low", float(props["clogp"]) >= thresholds["clogp_min"]),
        ("clogp_high", float(props["clogp"]) <= thresholds["clogp_max"]),
        ("hbd_high", int(props["hbd"]) <= thresholds["hbd_max"]),
        ("hba_high", int(props["hba"]) <= thresholds["hba_max"]),
        ("rotatable_bonds_high", int(props["rotatable_bonds"]) <= thresholds["rotatable_bonds_max"]),
        ("formal_charge_low", int(props["formal_charge"]) >= thresholds["formal_charge_min"]),
        ("formal_charge_high", int(props["formal_charge"]) <= thresholds["formal_charge_max"]),
        ("sa_score_high", float(props["sa_score"]) <= thresholds["sa_score_max"]),
    ]
    for reason, passed in checks:
        if not passed:
            failures.append(reason)

    clinical_similarity_pass = True
    if clinical_sim >= thresholds["clinical_similarity_max"]:
        clinical_similarity_pass = False
        failures.append("clinical_similarity_high")
    patent_proxy_soft_pass = patent_sim < thresholds["patent_proxy_soft_similarity_max"]
    if not patent_proxy_soft_pass:
        failures.append("patent_proxy_similarity_high_soft")

    property_pass = not any(
        reason
        not in {"clinical_similarity_high", "patent_proxy_similarity_high_soft"}
        for reason in failures
    )
    return property_pass, clinical_similarity_pass, patent_proxy_soft_pass, failures


def percentile(values: list[float], pct: float) -> float | None:
    if not values:
        return None
    ordered = sorted(values)
    rank = (len(ordered) - 1) * pct
    lower = math.floor(rank)
    upper = math.ceil(rank)
    if lower == upper:
        return ordered[int(rank)]
    return ordered[lower] * (upper - rank) + ordered[upper] * (rank - lower)


def numeric_summary(rows: list[dict[str, object]], field: str) -> dict[str, float | None]:
    values = [float(row[field]) for row in rows if row.get(field) not in {"", None}]
    return {
        "min": round(min(values), 3) if values else None,
        "p10": round(percentile(values, 0.10), 3) if values else None,
        "median": round(median(values), 3) if values else None,
        "p90": round(percentile(values, 0.90), 3) if values else None,
        "max": round(max(values), 3) if values else None,
    }


def calibration_recommendation(rows: list[dict[str, object]]) -> dict[str, object]:
    passed = [row for row in rows if row.get("passes_property_filters") == "true"]
    basis = passed or rows
    return {
        "basis": "property_pass_rows" if passed else "all_valid_rows",
        "recommended_gate": {
            "mw": [DEFAULT_THRESHOLDS["mw_min"], DEFAULT_THRESHOLDS["mw_max"]],
            "tpsa": [DEFAULT_THRESHOLDS["tpsa_min"], DEFAULT_THRESHOLDS["tpsa_max"]],
            "clogp": [DEFAULT_THRESHOLDS["clogp_min"], DEFAULT_THRESHOLDS["clogp_max"]],
            "formal_charge": [DEFAULT_THRESHOLDS["formal_charge_min"], DEFAULT_THRESHOLDS["formal_charge_max"]],
            "sa_score_max": DEFAULT_THRESHOLDS["sa_score_max"],
            "clinical_similarity_max": DEFAULT_THRESHOLDS["clinical_similarity_max"],
            "patent_example_similarity_max_when_curated": DEFAULT_THRESHOLDS["patent_proxy_similarity_max"],
            "patent_proxy_soft_similarity_max_for_current_seed_proxy": DEFAULT_THRESHOLDS["patent_proxy_soft_similarity_max"],
        },
        "observed_distribution": {field: numeric_summary(basis, field) for field in ["mw", "tpsa", "clogp", "formal_charge", "sa_score", "lipinski_rule_of_five_violations", "max_clinical_similarity", "max_patent_proxy_similarity"]},
    }


def filter_reinvent4_outputs(
    input_csv: str | Path,
    manifest_path: str | Path,
    output_dir: str | Path,
    thresholds: dict[str, float] | None = None,
    output_prefix: str | None = None,
) -> dict[str, object]:
    input_csv = Path(input_csv)
    manifest_path = Path(manifest_path)
    output_dir = Path(output_dir)
    thresholds = dict(DEFAULT_THRESHOLDS if thresholds is None else thresholds)
    output_prefix = output_prefix or input_csv.stem

    manifest = json.loads(manifest_path.read_text())
    clinical_fps = load_exclusion_fps(Path(manifest["outputs"]["clinical_exclusion"]))
    patent_fps = load_exclusion_fps(Path(manifest["outputs"]["patent_proxy_exclusion"]))

    input_rows = read_csv(input_csv)
    output_rows: list[dict[str, object]] = []
    seen_parent_keys: set[str] = set()
    duplicate_count = 0

    for row in input_rows:
        smiles = row.get("SMILES") or row.get("smiles") or row.get("sampled_smiles") or ""
        mol = Chem.MolFromSmiles(smiles)
        base = {
            "generated_smiles": smiles,
            "input_smiles": row.get("Input_SMILES", ""),
            "reinvent_tanimoto": parse_float(row.get("Tanimoto", "")),
            "nll": parse_float(row.get("NLL", "")),
        }
        if mol is None:
            output_rows.append({**base, "valid": "false", "failure_reasons": "invalid_smiles"})
            continue
        props = mol_props(mol)
        fp = AllChem.GetMorganFingerprintAsBitVect(parent_mol(smiles), 2, nBits=2048)
        clinical_sim = round(max_similarity(fp, clinical_fps), 3)
        patent_sim = round(max_similarity(fp, patent_fps), 3)
        property_pass, clinical_similarity_pass, patent_proxy_soft_pass, failures = assess_filters(props, clinical_sim, patent_sim, thresholds)
        lipinski_pass = props["passes_lipinski_rule_of_five"] == "true"
        if not lipinski_pass:
            failures.append("lipinski_rule_of_five_failed")
        if props["parent_inchi_key"] in seen_parent_keys:
            duplicate_count += 1
            failures.append("duplicate_parent_inchikey")
        seen_parent_keys.add(str(props["parent_inchi_key"]))
        output_rows.append(
            {
                **base,
                **props,
                "valid": "true",
                "max_clinical_similarity": clinical_sim,
                "max_patent_proxy_similarity": patent_sim,
                "passes_property_filters": str(property_pass).lower(),
                "passes_clinical_similarity_filter": str(clinical_similarity_pass).lower(),
                "passes_patent_proxy_soft_filter": str(patent_proxy_soft_pass).lower(),
                "passes_hard_filters": str(property_pass and lipinski_pass and clinical_similarity_pass and "duplicate_parent_inchikey" not in failures).lower(),
                "passes_all_filters": str(property_pass and lipinski_pass and clinical_similarity_pass and patent_proxy_soft_pass and "duplicate_parent_inchikey" not in failures).lower(),
                "failure_reasons": ";".join(failures),
            }
        )

    valid_rows = [row for row in output_rows if row.get("valid") == "true"]
    pass_rows = [row for row in valid_rows if row.get("passes_all_filters") == "true"]
    hard_pass_rows = [row for row in valid_rows if row.get("passes_hard_filters") == "true"]
    summary = {
        "input_csv": str(input_csv.resolve()),
        "total_rows": len(input_rows),
        "valid_rows": len(valid_rows),
        "invalid_rows": len(input_rows) - len(valid_rows),
        "unique_parent_inchikeys": len(seen_parent_keys),
        "duplicate_parent_inchikey_rows": duplicate_count,
        "property_pass_rows": sum(1 for row in valid_rows if row.get("passes_property_filters") == "true"),
        "clinical_similarity_pass_rows": sum(1 for row in valid_rows if row.get("passes_clinical_similarity_filter") == "true"),
        "patent_proxy_soft_pass_rows": sum(1 for row in valid_rows if row.get("passes_patent_proxy_soft_filter") == "true"),
        "lipinski_rule_of_five_pass_rows": sum(1 for row in valid_rows if row.get("passes_lipinski_rule_of_five") == "true"),
        "hard_filter_pass_rows": len(hard_pass_rows),
        "all_filter_pass_rows": len(pass_rows),
        "output_prefix": output_prefix,
        "thresholds": thresholds,
        "calibration": calibration_recommendation(valid_rows),
    }

    output_dir.mkdir(parents=True, exist_ok=True)
    filtered_csv = output_dir / f"{output_prefix}_filtered.csv"
    summary_json = output_dir / f"{output_prefix}_calibration.json"
    summary_md = output_dir / f"{output_prefix}_calibration.md"
    write_csv(filtered_csv, output_rows, FILTER_COLUMNS)
    summary_json.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    summary_md.write_text(markdown_summary(summary))
    return summary


def markdown_summary(summary: dict[str, object]) -> str:
    calibration = summary["calibration"]
    gate = calibration["recommended_gate"]
    dist = calibration["observed_distribution"]
    lines = [
        "# FXIa Mol2Mol Pilot Calibration",
        "",
        f"- Total REINVENT rows: {summary['total_rows']}",
        f"- Valid generated rows: {summary['valid_rows']}",
        f"- Unique parent InChIKeys: {summary['unique_parent_inchikeys']}",
        f"- Rows passing property filters: {summary['property_pass_rows']}",
        f"- Rows passing clinical similarity hard filter: {summary['clinical_similarity_pass_rows']}",
        f"- Rows passing current patent-proxy soft filter: {summary['patent_proxy_soft_pass_rows']}",
        f"- Rows passing property + clinical hard filters: {summary['hard_filter_pass_rows']}",
        f"- Rows passing all current filters including patent-proxy soft filter: {summary['all_filter_pass_rows']}",
        "",
        "## Recommended Gate For First Scale-Up",
        "",
        f"- MW: {gate['mw'][0]}-{gate['mw'][1]}",
        f"- TPSA: {gate['tpsa'][0]}-{gate['tpsa'][1]}",
        f"- cLogP: {gate['clogp'][0]}-{gate['clogp'][1]}",
        f"- Formal charge: {gate['formal_charge'][0]} to {gate['formal_charge'][1]}",
        f"- SA score: <= {gate['sa_score_max']}",
        f"- ECFP4 Tanimoto to clinical anchors: < {gate['clinical_similarity_max']}",
        f"- ECFP4 Tanimoto to curated patent examples when available: < {gate['patent_example_similarity_max_when_curated']}",
        f"- Current seed-derived patent-proxy soft triage: < {gate['patent_proxy_soft_similarity_max_for_current_seed_proxy']}",
        "",
        "## Observed Pilot Distributions",
        "",
        "| metric | min | p10 | median | p90 | max |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for field, values in dist.items():
        lines.append(f"| {field} | {values['min']} | {values['p10']} | {values['median']} | {values['p90']} | {values['max']} |")
    lines.extend(
        [
            "",
            "## Caveats",
            "",
            "- Patent-proxy similarity currently uses a broad public seed proxy and is a soft triage field, not a legal/IP rejection rule.",
            "- These filters are a first scale-up gate, not evidence of FXIa activity, selectivity, synthesizability, or freedom to operate.",
            "- Clinical candidates remain calibration/exclusion anchors only.",
        ]
    )
    return "\n".join(lines) + "\n"


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, dest="input_csv")
    parser.add_argument("--manifest", required=True, dest="manifest_path")
    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--output-prefix")
    args = parser.parse_args()
    summary = filter_reinvent4_outputs(args.input_csv, args.manifest_path, args.output_dir, output_prefix=args.output_prefix)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
