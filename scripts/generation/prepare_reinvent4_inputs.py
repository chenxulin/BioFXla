#!/usr/bin/env python
"""Prepare source-backed FXIa inputs for REINVENT4 pilot and scale-up runs."""

from __future__ import annotations

import argparse
import csv
import json
import re
from pathlib import Path
from typing import Iterable

from rdkit import Chem
from rdkit.Chem import BRICS
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.MolStandardize import rdMolStandardize


ROOT = Path(__file__).resolve().parents[2]
PROCESSED = ROOT / "data" / "processed"
REINVENT4 = ROOT / "tools" / "REINVENT4"
PRIORS = REINVENT4 / "priors"
MOL2MOL_PRIOR = PRIORS / "pubchem_ecfp4_with_count_with_rank_reinvent4_dict_voc.prior"
DENOVO_PRIOR = PRIORS / "reinvent_pubchem.prior"
LIBINVENT_PRIOR = PRIORS / "libinvent.prior"
LINKINVENT_PRIOR = PRIORS / "linkinvent.prior"

CLINICAL_IDS = {"CHEMBL5427131", "CHEMBL4112929"}
CLINICAL_NAME_TOKENS = {
    "asundexian",
    "bay 2433334",
    "bay-2433334",
    "milvexian",
    "bms-986177",
    "jnj-70033093",
}
PRIMARY_ENDPOINTS = {"Ki", "IC50", "Kd"}
PRIMARY_RELATIONS = {"=", "<", "<="}
LIBINVENT_SUPPORTED_BRACKET_TOKENS = {
    "[*]",
    "[*:0]",
    "[*:1]",
    "[N+]",
    "[N-]",
    "[N]",
    "[O-]",
    "[O]",
    "[S+]",
    "[n+]",
    "[nH]",
    "[s+]",
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


def write_lines(path: Path, lines: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def parent_mol(smiles: str) -> Chem.Mol | None:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None
    parent = rdMolStandardize.FragmentParent(mol)
    parent.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(parent)
    return parent


def parent_smiles(smiles: str) -> str:
    mol = parent_mol(smiles)
    if mol is None:
        return ""
    return Chem.MolToSmiles(mol, isomericSmiles=True)


def parent_key(smiles: str) -> str:
    mol = parent_mol(smiles)
    if mol is None:
        return ""
    return Chem.MolToInchiKey(mol)


def murcko_smiles(smiles: str) -> str:
    mol = parent_mol(smiles)
    if mol is None:
        return ""
    scaffold = MurckoScaffold.GetScaffoldForMol(mol)
    if scaffold is None or scaffold.GetNumAtoms() == 0:
        return ""
    return Chem.MolToSmiles(scaffold, isomericSmiles=True)


def is_clinical_anchor(row: dict[str, str], clinical_parent_keys: set[str]) -> bool:
    molecule_id = row.get("molecule_id", "")
    if molecule_id in CLINICAL_IDS:
        return True
    text = f"{row.get('name', '')} {row.get('synonyms', '')} {row.get('notes', '')}".lower()
    if any(token in text for token in CLINICAL_NAME_TOKENS):
        return True
    return bool(row.get("parent_inchi_key") and row["parent_inchi_key"] in clinical_parent_keys)


def get_clinical_rows(rows: list[dict[str, str]]) -> list[dict[str, str]]:
    clinical: dict[str, dict[str, str]] = {}
    for row in rows:
        molecule_id = row.get("molecule_id", "")
        text = f"{row.get('name', '')} {row.get('synonyms', '')} {row.get('notes', '')}".lower()
        if molecule_id in CLINICAL_IDS or any(token in text for token in CLINICAL_NAME_TOKENS):
            key = row.get("parent_inchi_key") or row.get("canonical_smiles") or molecule_id
            if key and key not in clinical:
                clinical[key] = row
    return list(clinical.values())


def potency_value(row: dict[str, str]) -> float | None:
    try:
        value = float(row.get("activity_value") or "")
    except ValueError:
        return None
    if value <= 0:
        return None
    return value


def candidate_rows(rows: list[dict[str, str]], clinical_parent_keys: set[str]) -> list[dict[str, object]]:
    candidates: list[dict[str, object]] = []
    seen_activity_records: set[str] = set()
    for row in rows:
        if row.get("source") != "ChEMBL":
            continue
        if is_clinical_anchor(row, clinical_parent_keys):
            continue
        if row.get("activity_type") not in PRIMARY_ENDPOINTS:
            continue
        if row.get("activity_unit", "").lower() != "nm":
            continue
        if row.get("activity_relation") not in PRIMARY_RELATIONS:
            continue
        value = potency_value(row)
        if value is None:
            continue
        smiles = parent_smiles(row.get("canonical_smiles", ""))
        key = row.get("parent_inchi_key") or parent_key(smiles)
        scaffold = murcko_smiles(smiles)
        if not smiles or not key or not scaffold:
            continue
        record_key = row.get("source_record_id", "")
        if record_key in seen_activity_records:
            continue
        seen_activity_records.add(record_key)
        enriched = dict(row)
        enriched.update({"parent_smiles": smiles, "parent_key": key, "murcko_smiles": scaffold, "activity_nm": value})
        candidates.append(enriched)
    return candidates


def best_rows_by_molecule(candidates: list[dict[str, object]]) -> list[dict[str, object]]:
    best: dict[str, dict[str, object]] = {}
    for row in candidates:
        key = str(row["parent_key"])
        if key not in best or float(row["activity_nm"]) < float(best[key]["activity_nm"]):
            best[key] = row
    return list(best.values())


def select_mol2mol_seeds(candidates: list[dict[str, object]], target_count: int = 40) -> list[dict[str, object]]:
    rows = sorted(best_rows_by_molecule(candidates), key=lambda row: (float(row["activity_nm"]), str(row["molecule_id"])))
    chosen: list[dict[str, object]] = []
    used_scaffolds: set[str] = set()
    for row in rows:
        scaffold = str(row["murcko_smiles"])
        if scaffold in used_scaffolds:
            continue
        chosen.append(row)
        used_scaffolds.add(scaffold)
        if len(chosen) >= target_count:
            return chosen
    for row in rows:
        if row in chosen:
            continue
        chosen.append(row)
        if len(chosen) >= target_count:
            break
    return chosen


def make_libinvent_template(scaffold_smiles: str) -> str:
    mol = Chem.MolFromSmiles(scaffold_smiles)
    if mol is None or mol.GetNumAtoms() < 3:
        return ""
    atoms = [
        atom.GetIdx()
        for atom in mol.GetAtoms()
        if atom.GetAtomicNum() > 1 and atom.GetTotalNumHs() > 0 and atom.GetDegree() <= 3
    ]
    if len(atoms) < 2:
        atoms = [atom.GetIdx() for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1 and atom.GetDegree() <= 3]
    if len(atoms) < 2:
        return ""
    first, second = atoms[0], atoms[-1]
    rw_mol = Chem.RWMol(mol)
    for map_num, anchor_idx in enumerate([first, second]):
        dummy = Chem.Atom(0)
        dummy.SetAtomMapNum(map_num)
        dummy_idx = rw_mol.AddAtom(dummy)
        rw_mol.AddBond(anchor_idx, dummy_idx, Chem.BondType.SINGLE)
    template = rw_mol.GetMol()
    try:
        Chem.SanitizeMol(template)
    except Exception:
        return ""
    smiles = Chem.MolToSmiles(template, isomericSmiles=False)
    bracket_tokens = set(re.findall(r"\[[^\]]+\]", smiles))
    if bracket_tokens - LIBINVENT_SUPPORTED_BRACKET_TOKENS:
        return ""
    return smiles if smiles.count("*") == 2 else ""


def select_libinvent_templates(seeds: list[dict[str, object]], min_count: int = 10) -> list[dict[str, str]]:
    templates: list[dict[str, str]] = []
    seen: set[str] = set()
    for row in seeds:
        scaffold = str(row["murcko_smiles"])
        if scaffold in seen:
            continue
        template = make_libinvent_template(scaffold)
        if not template:
            continue
        templates.append(
            {
                "scaffold_smiles": scaffold,
                "template_smiles": template,
                "source_molecule_id": str(row["molecule_id"]),
                "source_activity_nm": str(row["activity_nm"]),
                "confidence": "low",
                "notes": "heuristic two-point LibInvent template; attachment vectors require medchem review",
            }
        )
        seen.add(scaffold)
        if len(templates) >= min_count:
            break
    return templates


def normalize_brics_fragment(fragment: str) -> str:
    normalized = re.sub(r"\[\d+\*\]", "*", fragment)
    return normalized if normalized.count("*") == 1 else ""


def select_linkinvent_pairs(seeds: list[dict[str, object]], max_count: int = 10) -> list[dict[str, str]]:
    pairs: list[dict[str, str]] = []
    seen: set[str] = set()
    for row in seeds:
        mol = Chem.MolFromSmiles(str(row["parent_smiles"]))
        if mol is None:
            continue
        fragments = [normalize_brics_fragment(fragment) for fragment in BRICS.BRICSDecompose(mol)]
        valid = []
        for fragment in fragments:
            if not fragment:
                continue
            parsed = Chem.MolFromSmiles(fragment)
            if parsed is None or parsed.GetNumHeavyAtoms() < 5:
                continue
            valid.append(fragment)
        valid = sorted(set(valid), key=lambda frag: Chem.MolFromSmiles(frag).GetNumHeavyAtoms(), reverse=True)
        if len(valid) < 2:
            continue
        pair = f"{valid[0]}|{valid[1]}"
        if pair in seen:
            continue
        pairs.append(
            {
                "warhead_pair": pair,
                "source_molecule_id": str(row["molecule_id"]),
                "source_activity_nm": str(row["activity_nm"]),
                "confidence": "low",
                "notes": "BRICS-derived LinkInvent pair; requires binding-mode and synthetic review before scale-up",
            }
        )
        seen.add(pair)
        if len(pairs) >= max_count:
            break
    return pairs


def toml_sampling(
    *,
    device: str,
    json_out_config: Path,
    model_file: Path,
    output_file: Path,
    num_smiles: int,
    smiles_file: Path | None = None,
    sample_strategy: str | None = None,
    temperature: float | None = None,
    randomize_smiles: bool = True,
) -> str:
    lines = [
        'run_type = "sampling"',
        f'device = "{device}"',
        f'json_out_config = "{json_out_config}"',
        "",
        "[parameters]",
        f'model_file = "{model_file}"',
    ]
    if smiles_file is not None:
        lines.append(f'smiles_file = "{smiles_file}"')
    if sample_strategy is not None:
        lines.append(f'sample_strategy = "{sample_strategy}"')
    if temperature is not None:
        lines.append(f"temperature = {temperature}")
    lines.extend(
        [
            f'output_file = "{output_file}"',
            f"num_smiles = {num_smiles}",
            "unique_molecules = true",
            f"randomize_smiles = {str(randomize_smiles).lower()}",
            "",
        ]
    )
    return "\n".join(lines)


def prepare_generation_inputs(output_dir: str | Path = ROOT / "results" / "reinvent4_fxia_pilot") -> dict[str, object]:
    output_dir = Path(output_dir).resolve()
    molecule_rows = read_csv(PROCESSED / "fxia_molecule_seed.csv")

    clinical_rows = get_clinical_rows(molecule_rows)
    clinical_parent_keys = {row.get("parent_inchi_key", "") for row in clinical_rows if row.get("parent_inchi_key")}
    candidates = candidate_rows(molecule_rows, clinical_parent_keys)
    seeds = select_mol2mol_seeds(candidates, target_count=40)
    templates = select_libinvent_templates(seeds, min_count=10)
    link_pairs = select_linkinvent_pairs(seeds, max_count=10)

    potent_proxy = sorted(best_rows_by_molecule(candidates), key=lambda row: float(row["activity_nm"]))[:100]
    patent_proxy_rows = clinical_rows + [dict(row) for row in potent_proxy]

    mol2mol_dir = output_dir / "mol2mol"
    libinvent_dir = output_dir / "libinvent"
    linkinvent_dir = output_dir / "linkinvent"
    exclusion_dir = output_dir / "exclusions"
    config_dir = output_dir / "configs"
    raw_dir = output_dir / "raw_scale"

    write_lines(mol2mol_dir / "fxia_mol2mol_pilot.smi", [str(row["parent_smiles"]) for row in seeds])
    write_lines(libinvent_dir / "fxia_libinvent_scaffolds.smi", [row["template_smiles"] for row in templates])
    write_lines(linkinvent_dir / "fxia_linkinvent_warheads.smi", [row["warhead_pair"] for row in link_pairs])
    write_lines(exclusion_dir / "fxia_clinical_similarity_exclusion.smi", [parent_smiles(row["canonical_smiles"]) for row in clinical_rows])
    write_lines(
        exclusion_dir / "fxia_patent_similarity_exclusion.smi",
        [parent_smiles(row["canonical_smiles"]) if row.get("canonical_smiles") else str(row.get("parent_smiles", "")) for row in patent_proxy_rows],
    )

    seed_columns = [
        "molecule_id",
        "parent_smiles",
        "parent_key",
        "murcko_smiles",
        "activity_type",
        "activity_nm",
        "assay_type",
        "source_record_id",
        "source_url",
        "doi",
        "confidence",
        "notes",
    ]
    write_csv(mol2mol_dir / "fxia_mol2mol_pilot_metadata.csv", seeds, seed_columns)
    write_csv(libinvent_dir / "fxia_libinvent_scaffolds_metadata.csv", templates, list(templates[0].keys()) if templates else [])
    write_csv(linkinvent_dir / "fxia_linkinvent_warheads_metadata.csv", link_pairs, list(link_pairs[0].keys()) if link_pairs else [])

    clinical_meta = [
        {
            "molecule_id": row.get("molecule_id", ""),
            "name": row.get("name", ""),
            "synonyms": row.get("synonyms", ""),
            "parent_smiles": parent_smiles(row.get("canonical_smiles", "")),
            "parent_inchi_key": row.get("parent_inchi_key", ""),
            "source_url": row.get("source_url", ""),
            "notes": row.get("notes", ""),
        }
        for row in clinical_rows
    ]
    write_csv(
        exclusion_dir / "fxia_clinical_similarity_exclusion_metadata.csv",
        clinical_meta,
        ["molecule_id", "name", "synonyms", "parent_smiles", "parent_inchi_key", "source_url", "notes"],
    )

    patent_meta = []
    for row in patent_proxy_rows:
        smiles = parent_smiles(row.get("canonical_smiles", "")) if row.get("canonical_smiles") else str(row.get("parent_smiles", ""))
        patent_meta.append(
            {
                "molecule_id": row.get("molecule_id", ""),
                "name": row.get("name", ""),
                "parent_smiles": smiles,
                "activity_type": row.get("activity_type", ""),
                "activity_nm": row.get("activity_nm", row.get("activity_value", "")),
                "source_url": row.get("source_url", ""),
                "confidence": "low" if row not in clinical_rows else row.get("confidence", ""),
                "notes": "proxy exclusion; replace or augment with curated patent examples before legal/IP conclusions",
            }
        )
    write_csv(
        exclusion_dir / "fxia_patent_similarity_exclusion_metadata.csv",
        patent_meta,
        ["molecule_id", "name", "parent_smiles", "activity_type", "activity_nm", "source_url", "confidence", "notes"],
    )

    mol2mol_pilot_config = config_dir / "mol2mol_pilot.toml"
    mol2mol_scale_config = config_dir / "mol2mol_scale.toml"
    libinvent_scale_config = config_dir / "libinvent_scale.toml"
    linkinvent_scale_config = config_dir / "linkinvent_scale.toml"
    denovo_scale_config = config_dir / "reinvent_denovo_scale.toml"

    config_dir.mkdir(parents=True, exist_ok=True)
    mol2mol_pilot_config.write_text(
        toml_sampling(
            device="cpu",
            json_out_config=config_dir / "mol2mol_pilot.json",
            model_file=MOL2MOL_PRIOR,
            smiles_file=mol2mol_dir / "fxia_mol2mol_pilot.smi",
            output_file=mol2mol_dir / "mol2mol_pilot.csv",
            num_smiles=25,
            sample_strategy="beamsearch",
        )
    )
    mol2mol_scale_config.write_text(
        toml_sampling(
            device="cuda:0",
            json_out_config=config_dir / "mol2mol_scale.json",
            model_file=MOL2MOL_PRIOR,
            smiles_file=mol2mol_dir / "fxia_mol2mol_pilot.smi",
            output_file=raw_dir / "mol2mol_scale.csv",
            num_smiles=2500,
            sample_strategy="multinomial",
            temperature=1.0,
        )
    )
    libinvent_scale_config.write_text(
        toml_sampling(
            device="cuda:0",
            json_out_config=config_dir / "libinvent_scale.json",
            model_file=LIBINVENT_PRIOR,
            smiles_file=libinvent_dir / "fxia_libinvent_scaffolds.smi",
            output_file=raw_dir / "libinvent_scale.csv",
            num_smiles=10000,
            sample_strategy="multinomial",
            temperature=1.0,
        )
    )
    linkinvent_scale_config.write_text(
        toml_sampling(
            device="cuda:0",
            json_out_config=config_dir / "linkinvent_scale.json",
            model_file=LINKINVENT_PRIOR,
            smiles_file=linkinvent_dir / "fxia_linkinvent_warheads.smi",
            output_file=raw_dir / "linkinvent_scale.csv",
            num_smiles=10000,
            sample_strategy="multinomial",
            temperature=1.0,
        )
    )
    denovo_scale_config.write_text(
        toml_sampling(
            device="cuda:0",
            json_out_config=config_dir / "reinvent_denovo_scale.json",
            model_file=DENOVO_PRIOR,
            output_file=raw_dir / "reinvent_denovo_scale.csv",
            num_smiles=100000,
            randomize_smiles=True,
        )
    )

    summary = {
        "mol2mol_seed_count": len(seeds),
        "libinvent_template_count": len(templates),
        "linkinvent_pair_count": len(link_pairs),
        "clinical_exclusion_count": len(clinical_rows),
        "patent_proxy_exclusion_count": len(patent_meta),
        "planned_scale_raw_counts": {
            "mol2mol": len(seeds) * 2500,
            "libinvent": len(templates) * 10000,
            "linkinvent": len(link_pairs) * 10000,
            "reinvent_denovo": 100000,
        },
        "local_priors": {
            "mol2mol": str(MOL2MOL_PRIOR.resolve()),
            "libinvent": str(LIBINVENT_PRIOR.resolve()),
            "linkinvent": str(LINKINVENT_PRIOR.resolve()),
            "reinvent_denovo": str(DENOVO_PRIOR.resolve()),
        },
        "input_source_tables": [
            str((PROCESSED / "fxia_molecule_seed.csv").resolve()),
            str((PROCESSED / "fxia_scaffold_seed.csv").resolve()),
        ],
        "outputs": {
            "mol2mol_pilot_smiles": str((mol2mol_dir / "fxia_mol2mol_pilot.smi").resolve()),
            "libinvent_scaffolds": str((libinvent_dir / "fxia_libinvent_scaffolds.smi").resolve()),
            "linkinvent_warheads": str((linkinvent_dir / "fxia_linkinvent_warheads.smi").resolve()),
            "clinical_exclusion": str((exclusion_dir / "fxia_clinical_similarity_exclusion.smi").resolve()),
            "patent_proxy_exclusion": str((exclusion_dir / "fxia_patent_similarity_exclusion.smi").resolve()),
            "mol2mol_pilot_config": str(mol2mol_pilot_config.resolve()),
            "mol2mol_scale_config": str(mol2mol_scale_config.resolve()),
            "libinvent_scale_config": str(libinvent_scale_config.resolve()),
            "linkinvent_scale_config": str(linkinvent_scale_config.resolve()),
            "reinvent_denovo_scale_config": str(denovo_scale_config.resolve()),
        },
        "caveats": [
            "Clinical candidates are exclusion/calibration anchors only, not direct design templates.",
            "Patent exclusion set is a public-data proxy until curated patent examples are collected.",
            "LibInvent and LinkInvent attachment points are heuristic and require medchem/PDB review before scale-up.",
            "Scale configs are starting points; run only after pilot filtering calibrates thresholds.",
        ],
    }
    (output_dir / "generation_input_manifest.json").write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n")
    (output_dir / "README.md").write_text(readme_text(summary))
    return summary


def readme_text(summary: dict[str, object]) -> str:
    counts = summary["planned_scale_raw_counts"]
    return f"""# FXIa REINVENT4 Pilot Package

This package is generated from `data/processed/fxia_molecule_seed.csv` and keeps clinical candidates as exclusion/calibration anchors only.

## Current Counts

- Mol2Mol non-clinical pilot seeds: {summary["mol2mol_seed_count"]}
- LibInvent heuristic scaffold templates: {summary["libinvent_template_count"]}
- LinkInvent BRICS-derived warhead pairs: {summary["linkinvent_pair_count"]}
- Clinical similarity exclusion anchors: {summary["clinical_exclusion_count"]}
- Patent-proxy similarity exclusions: {summary["patent_proxy_exclusion_count"]}

## Pilot Command

```bash
mamba run -n aidd reinvent -l results/reinvent4_fxia_pilot/mol2mol/mol2mol_pilot.log results/reinvent4_fxia_pilot/configs/mol2mol_pilot.toml
```

## Scale-Up Configs

- Mol2Mol: `configs/mol2mol_scale.toml` ({counts["mol2mol"]} planned raw molecules)
- LibInvent: `configs/libinvent_scale.toml` ({counts["libinvent"]} planned raw molecules)
- LinkInvent: `configs/linkinvent_scale.toml` ({counts["linkinvent"]} planned raw molecules)
- REINVENT de novo: `configs/reinvent_denovo_scale.toml` ({counts["reinvent_denovo"]} planned raw molecules)

Scale-up command:

```bash
mamba run -n aidd bash scripts/generation/run_reinvent4_scaleup.sh
```

This repository-level script starts or resumes all four paths, skips paths whose logs already contain `Finished REINVENT`, and writes PID/stdout/stderr files under `raw_scale/`. The generated configs in this package are inputs to that script.

Do not treat the scale configs as medicinal chemistry recommendations. Run them only after the pilot output has been filtered for MW, TPSA, cLogP, charge, SA, clinical-candidate similarity, patent-proxy similarity, and reactive/PAINS liabilities.
"""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", default=ROOT / "results" / "reinvent4_fxia_pilot")
    args = parser.parse_args()
    summary = prepare_generation_inputs(args.output_dir)
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
