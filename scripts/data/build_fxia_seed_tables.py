#!/usr/bin/env python
"""Build source-backed seed tables for the FXIa scaffold atlas.

This script intentionally keeps assay rows separate. It standardizes structures
with RDKit for grouping, but it does not collapse Ki/IC50/Kd values across assay
contexts.
"""

from __future__ import annotations

import csv
import glob
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path

from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.MolStandardize import rdMolStandardize


ROOT = Path(__file__).resolve().parents[2]
RAW = ROOT / "data" / "raw_public"
OUT = ROOT / "data" / "processed"

MOLECULE_COLUMNS = [
    "molecule_id",
    "source",
    "source_record_id",
    "name",
    "synonyms",
    "canonical_smiles",
    "inchi_key",
    "parent_inchi_key",
    "stereochemistry_defined",
    "salt_stripped",
    "target",
    "activity_type",
    "activity_value",
    "activity_unit",
    "activity_relation",
    "p_activity",
    "assay_type",
    "assay_description",
    "species",
    "selectivity_fx_a",
    "selectivity_thrombin",
    "selectivity_kallikrein",
    "selectivity_trypsin",
    "functional_aptt",
    "functional_pt",
    "adme_notes",
    "pk_notes",
    "tox_notes",
    "source_url",
    "doi",
    "confidence",
    "notes",
]

PDB_INTERACTION_COLUMNS = [
    "pdb_id",
    "ligand_id",
    "ligand_name",
    "ligand_smiles",
    "resolution",
    "protein_species",
    "binding_site_region",
    "s1_occupancy",
    "s2_s4_occupancy",
    "solvent_front_exit_vector",
    "key_residue_contacts",
    "water_mediated_contacts",
    "pose_quality_notes",
    "source_url",
    "doi",
]

CLINICAL_TRIAL_COLUMNS = [
    "nct_id",
    "query_term",
    "representative_molecule",
    "brief_title",
    "overall_status",
    "phase",
    "start_date",
    "completion_date",
    "interventions",
    "source_url",
    "confidence",
    "notes",
]

SCAFFOLD_COLUMNS = [
    "scaffold_id",
    "scaffold_name",
    "murcko_smiles",
    "generic_murcko_smiles",
    "source_classes",
    "representative_molecules",
    "best_fxia_activity",
    "median_fxia_activity",
    "activity_tier",
    "selectivity_evidence",
    "pdb_evidence",
    "sar_depth",
    "patent_crowding",
    "synthetic_accessibility",
    "generation_suitability",
    "recommended_generation_methods",
    "key_exit_vectors",
    "known_liabilities",
    "priority_score",
    "decision",
    "notes",
]


def load_json(path: Path) -> dict:
    with path.open() as fh:
        return json.load(fh)


def canonicalize(smiles: str | None) -> dict[str, str | bool | None]:
    if not smiles:
        return {
            "canonical_smiles": None,
            "inchi_key": None,
            "parent_inchi_key": None,
            "salt_stripped": None,
            "stereochemistry_defined": None,
            "murcko_smiles": None,
            "generic_murcko_smiles": None,
            "notes": "missing_smiles",
        }

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {
            "canonical_smiles": smiles,
            "inchi_key": None,
            "parent_inchi_key": None,
            "salt_stripped": None,
            "stereochemistry_defined": None,
            "murcko_smiles": None,
            "generic_murcko_smiles": None,
            "notes": "rdkit_parse_failed",
        }

    canonical = Chem.MolToSmiles(mol, isomericSmiles=True)
    inchi_key = Chem.MolToInchiKey(mol)
    parent = rdMolStandardize.FragmentParent(mol)
    parent.UpdatePropertyCache(strict=False)
    Chem.GetSymmSSSR(parent)
    parent_smiles = Chem.MolToSmiles(parent, isomericSmiles=True)
    parent_inchi_key = Chem.MolToInchiKey(parent)
    salt_stripped = parent_smiles != canonical
    stereo_defined = any(atom.HasProp("_CIPCode") for atom in parent.GetAtoms())

    murcko = MurckoScaffold.GetScaffoldForMol(parent)
    murcko_smiles = Chem.MolToSmiles(murcko, isomericSmiles=True) if murcko else ""
    generic = MurckoScaffold.MakeScaffoldGeneric(murcko) if murcko else None
    generic_smiles = Chem.MolToSmiles(generic, isomericSmiles=False) if generic else ""

    return {
        "canonical_smiles": canonical,
        "inchi_key": inchi_key,
        "parent_inchi_key": parent_inchi_key,
        "salt_stripped": str(salt_stripped).lower(),
        "stereochemistry_defined": str(stereo_defined).lower(),
        "murcko_smiles": murcko_smiles,
        "generic_murcko_smiles": generic_smiles,
        "notes": "",
    }


def p_activity_from_row(row: dict) -> str:
    if row.get("pchembl_value"):
        return str(row["pchembl_value"])
    try:
        value = float(row.get("standard_value") or "")
    except ValueError:
        return ""
    units = (row.get("standard_units") or "").lower()
    if value <= 0:
        return ""
    if units == "nm":
        return f"{-math.log10(value * 1e-9):.3f}"
    if units in {"um", "\u00b5m"}:
        return f"{-math.log10(value * 1e-6):.3f}"
    if units == "m":
        return f"{-math.log10(value):.3f}"
    return ""


def activity_tier(best_nm: float | None, n_compounds: int) -> str:
    if best_nm is None:
        return "Tier D"
    if best_nm <= 10 and n_compounds >= 2:
        return "Tier A"
    if best_nm <= 100 and n_compounds >= 2:
        return "Tier B"
    if best_nm <= 1000:
        return "Tier C"
    return "Tier D"


def write_csv(path: Path, columns: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def build_chembl_activity_tables() -> tuple[list[dict], list[dict]]:
    files = sorted(RAW.glob("chembl_fxia_activities_page*.json"))
    molecule_rows: list[dict] = []
    scaffold_groups: dict[str, list[dict]] = defaultdict(list)

    for path in files:
        data = load_json(path)
        for act in data.get("activities", []):
            endpoint = act.get("standard_type") or act.get("type") or ""
            smiles = act.get("canonical_smiles")
            chem = canonicalize(smiles)
            notes = [chem["notes"]] if chem["notes"] else []
            if endpoint not in {"Ki", "IC50", "Kd"}:
                notes.append("non_primary_potency_endpoint_retained_but_not_scaffold_scored")

            source_record_id = "|".join(
                filter(
                    None,
                    [
                        act.get("activity_id") and f"activity:{act.get('activity_id')}",
                        act.get("assay_chembl_id"),
                        act.get("document_chembl_id"),
                    ],
                )
            )
            row = {
                "molecule_id": act.get("molecule_chembl_id") or "",
                "source": "ChEMBL",
                "source_record_id": source_record_id,
                "name": act.get("molecule_pref_name") or "",
                "synonyms": "",
                "canonical_smiles": chem["canonical_smiles"] or smiles or "",
                "inchi_key": chem["inchi_key"] or "",
                "parent_inchi_key": chem["parent_inchi_key"] or "",
                "stereochemistry_defined": chem["stereochemistry_defined"] or "",
                "salt_stripped": chem["salt_stripped"] or "",
                "target": "CHEMBL2820; human coagulation factor XI/F11; UniProt P03951",
                "activity_type": endpoint,
                "activity_value": act.get("standard_value") or act.get("value") or "",
                "activity_unit": act.get("standard_units") or act.get("units") or "",
                "activity_relation": act.get("standard_relation") or act.get("relation") or "",
                "p_activity": p_activity_from_row(act),
                "assay_type": act.get("assay_type") or "",
                "assay_description": act.get("assay_description") or "",
                "species": act.get("target_organism") or "",
                "selectivity_fx_a": "",
                "selectivity_thrombin": "",
                "selectivity_kallikrein": "",
                "selectivity_trypsin": "",
                "functional_aptt": "",
                "functional_pt": "",
                "adme_notes": "",
                "pk_notes": "",
                "tox_notes": "",
                "source_url": (
                    f"https://www.ebi.ac.uk/chembl/activity_report_card/"
                    f"{act.get('activity_id')}"
                    if act.get("activity_id")
                    else "https://www.ebi.ac.uk/chembl/target_report_card/CHEMBL2820/"
                ),
                "doi": "",
                "confidence": "medium",
                "notes": "; ".join(notes),
                "_murcko_smiles": chem["murcko_smiles"] or "",
                "_generic_murcko_smiles": chem["generic_murcko_smiles"] or "",
            }
            molecule_rows.append(row)

            if row["_murcko_smiles"]:
                scaffold_groups[row["_murcko_smiles"]].append(row)

    scaffold_rows: list[dict] = []
    for idx, (murcko, rows) in enumerate(
        sorted(scaffold_groups.items(), key=lambda item: len(item[1]), reverse=True),
        start=1,
    ):
        scored_nm = []
        for row in rows:
            if row["activity_type"] not in {"Ki", "IC50", "Kd"}:
                continue
            if row["activity_relation"] not in {"=", "<", "<="}:
                continue
            if row["activity_unit"] != "nM":
                continue
            try:
                scored_nm.append(float(row["activity_value"]))
            except ValueError:
                continue
        best = min(scored_nm) if scored_nm else None
        median = statistics.median(scored_nm) if scored_nm else None
        reps = []
        seen = set()
        for row in rows:
            mol_id = row["molecule_id"]
            if mol_id and mol_id not in seen:
                reps.append(mol_id)
                seen.add(mol_id)
            if len(reps) >= 5:
                break
        n_compounds = len({row["parent_inchi_key"] for row in rows if row["parent_inchi_key"]})
        sar_depth = "low_sar_depth" if n_compounds < 3 else "chembl_activity_cluster"
        scaffold_rows.append(
            {
                "scaffold_id": f"FXIA-SCAF-{idx:04d}",
                "scaffold_name": "ChEMBL Murcko scaffold cluster",
                "murcko_smiles": murcko,
                "generic_murcko_smiles": rows[0]["_generic_murcko_smiles"],
                "source_classes": "ChEMBL activity records",
                "representative_molecules": ";".join(reps),
                "best_fxia_activity": f"{best:g} nM" if best is not None else "needs_verification",
                "median_fxia_activity": f"{median:g} nM" if median is not None else "needs_verification",
                "activity_tier": activity_tier(best, n_compounds),
                "selectivity_evidence": "needs_verification",
                "pdb_evidence": "needs_verification",
                "sar_depth": sar_depth,
                "patent_crowding": "needs_patent_review",
                "synthetic_accessibility": "needs_verification",
                "generation_suitability": "needs_verification",
                "recommended_generation_methods": "needs_verification",
                "key_exit_vectors": "needs_verification",
                "known_liabilities": "selectivity and assay context not yet reviewed",
                "priority_score": "",
                "decision": "hold_for_more_data",
                "notes": f"{len(rows)} ChEMBL assay rows; {n_compounds} unique parent InChIKeys. Primary potency endpoints only used for activity tier.",
            }
        )

    return molecule_rows, scaffold_rows


def build_clinical_representatives() -> list[dict]:
    reps = [
        {
            "name": "Asundexian",
            "synonyms": "BAY 2433334; BAY-2433334",
            "source_file": "pubchem_asundexian_properties.json",
            "chembl_file": "chembl_molecule_asundexian_search.json",
            "chembl_id": "CHEMBL5427131",
            "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/135206011",
            "notes": "clinical-candidate scaffold anchor; patent-crowding assumed high until patent review",
        },
        {
            "name": "Milvexian",
            "synonyms": "BMS-986177; JNJ-70033093",
            "source_file": "pubchem_milvexian_properties.json",
            "chembl_file": "chembl_molecule_milvexian_search.json",
            "chembl_id": "CHEMBL4112929",
            "source_url": "https://pubchem.ncbi.nlm.nih.gov/compound/118277544",
            "notes": "clinical-candidate scaffold anchor; patent-crowding assumed high until patent review",
        },
    ]
    rows = []
    for rep in reps:
        props = load_json(RAW / rep["source_file"])["PropertyTable"]["Properties"][0]
        chem = canonicalize(props.get("SMILES"))
        rows.append(
            {
                "molecule_id": rep["chembl_id"],
                "source": "PubChem; ChEMBL",
                "source_record_id": f"PubChem CID {props.get('CID')}; {rep['chembl_id']}",
                "name": rep["name"],
                "synonyms": rep["synonyms"],
                "canonical_smiles": chem["canonical_smiles"] or props.get("SMILES", ""),
                "inchi_key": chem["inchi_key"] or props.get("InChIKey", ""),
                "parent_inchi_key": chem["parent_inchi_key"] or "",
                "stereochemistry_defined": chem["stereochemistry_defined"] or "",
                "salt_stripped": chem["salt_stripped"] or "",
                "target": "FXIa / F11",
                "activity_type": "",
                "activity_value": "",
                "activity_unit": "",
                "activity_relation": "",
                "p_activity": "",
                "assay_type": "",
                "assay_description": "",
                "species": "",
                "selectivity_fx_a": "",
                "selectivity_thrombin": "",
                "selectivity_kallikrein": "",
                "selectivity_trypsin": "",
                "functional_aptt": "",
                "functional_pt": "",
                "adme_notes": "",
                "pk_notes": "",
                "tox_notes": "",
                "source_url": rep["source_url"],
                "doi": "",
                "confidence": "high",
                "notes": rep["notes"],
            }
        )
    return rows


def build_pdb_id_table() -> list[dict]:
    target = load_json(RAW / "chembl_target_factor_xia_search.json")
    rows = []
    seen = set()
    for tgt in target.get("targets", []):
        if tgt.get("target_chembl_id") != "CHEMBL2820":
            continue
        for component in tgt.get("target_components", []):
            for xref in component.get("target_component_xrefs", []):
                if xref.get("xref_src_db") != "PDB":
                    continue
                pdb_id = xref.get("xref_id")
                if pdb_id and pdb_id not in seen:
                    seen.add(pdb_id)
                    rows.append(
                        {
                            "pdb_id": pdb_id,
                            "source": "ChEMBL target cross-reference",
                            "target": "CHEMBL2820; UniProt P03951",
                            "source_url": f"https://www.rcsb.org/structure/{pdb_id}",
                            "confidence": "medium",
                            "notes": "binding ligand/resolution/interaction metadata still needs RCSB/PDB extraction",
                        }
                    )
    return sorted(rows, key=lambda row: row["pdb_id"])


def build_pdb_interaction_table() -> list[dict]:
    rows: list[dict] = []
    chemcomp_cache: dict[str, dict] = {}

    for entry_path in sorted((RAW / "rcsb_entries").glob("*.json")):
        entry = load_json(entry_path)
        pdb_id = entry.get("rcsb_id") or entry_path.stem
        resolution = entry.get("rcsb_entry_info", {}).get("resolution_combined") or []
        resolution_text = ";".join(str(x) for x in resolution)
        doi = entry.get("rcsb_primary_citation", {}).get("pdbx_database_id_DOI") or ""
        entity_ids = (
            entry.get("rcsb_entry_container_identifiers", {}).get("non_polymer_entity_ids")
            or []
        )

        for entity_id in entity_ids:
            nonpoly_path = RAW / "rcsb_nonpolymer" / f"{pdb_id}_{entity_id}.json"
            if not nonpoly_path.exists():
                rows.append(
                    {
                        "pdb_id": pdb_id,
                        "ligand_id": "",
                        "ligand_name": "",
                        "ligand_smiles": "",
                        "resolution": resolution_text,
                        "protein_species": "Homo sapiens; inferred from CHEMBL2820/P03951 cross-reference",
                        "binding_site_region": "needs_verification",
                        "s1_occupancy": "needs_verification",
                        "s2_s4_occupancy": "needs_verification",
                        "solvent_front_exit_vector": "needs_verification",
                        "key_residue_contacts": "needs_verification",
                        "water_mediated_contacts": "needs_verification",
                        "pose_quality_notes": f"missing RCSB nonpolymer entity file for entity {entity_id}",
                        "source_url": f"https://www.rcsb.org/structure/{pdb_id}",
                        "doi": doi,
                    }
                )
                continue

            nonpoly = load_json(nonpoly_path)
            ligand_id = nonpoly.get("pdbx_entity_nonpoly", {}).get("comp_id") or ""
            ligand_name = nonpoly.get("pdbx_entity_nonpoly", {}).get("name") or ""
            ligand_smiles = ""
            ligand_type = ""
            if ligand_id:
                if ligand_id not in chemcomp_cache:
                    chemcomp_path = RAW / "rcsb_chemcomp" / f"{ligand_id}.json"
                    chemcomp_cache[ligand_id] = (
                        load_json(chemcomp_path) if chemcomp_path.exists() else {}
                    )
                chemcomp = chemcomp_cache[ligand_id]
                descriptor = chemcomp.get("rcsb_chem_comp_descriptor", {})
                ligand_smiles = (
                    descriptor.get("SMILES_stereo")
                    or descriptor.get("SMILES")
                    or descriptor.get("InChI")
                    or ""
                )
                ligand_type = chemcomp.get("chem_comp", {}).get("type") or ""

            pose_notes = "RCSB metadata only; binding-site contacts need structure interaction extraction"
            if ligand_type:
                pose_notes += f"; chem_comp_type={ligand_type}"

            rows.append(
                {
                    "pdb_id": pdb_id,
                    "ligand_id": ligand_id,
                    "ligand_name": ligand_name,
                    "ligand_smiles": ligand_smiles,
                    "resolution": resolution_text,
                    "protein_species": "Homo sapiens; inferred from CHEMBL2820/P03951 cross-reference",
                    "binding_site_region": "needs_verification",
                    "s1_occupancy": "needs_verification",
                    "s2_s4_occupancy": "needs_verification",
                    "solvent_front_exit_vector": "needs_verification",
                    "key_residue_contacts": "needs_verification",
                    "water_mediated_contacts": "needs_verification",
                    "pose_quality_notes": pose_notes,
                    "source_url": f"https://www.rcsb.org/structure/{pdb_id}",
                    "doi": doi,
                }
            )

    return rows


def build_clinical_trial_table() -> list[dict]:
    files = {
        "asundexian": "clinicaltrials_asundexian.json",
        "BAY 2433334": "clinicaltrials_bay_2433334.json",
        "milvexian": "clinicaltrials_milvexian.json",
        "BMS-986177": "clinicaltrials_bms_986177.json",
    }
    rows_by_nct: dict[str, dict] = {}

    for query_term, filename in files.items():
        path = RAW / filename
        if not path.exists():
            continue
        data = load_json(path)
        for study in data.get("studies", []):
            protocol = study.get("protocolSection", {})
            ident = protocol.get("identificationModule", {})
            status = protocol.get("statusModule", {})
            design = protocol.get("designModule", {})
            arms = protocol.get("armsInterventionsModule", {})
            nct_id = ident.get("nctId") or ""
            if not nct_id:
                continue
            interventions = []
            for intervention in arms.get("interventions", []) or []:
                name = intervention.get("name")
                if name:
                    interventions.append(name)
            molecule = (
                "Asundexian"
                if "asundexian" in query_term.lower() or "bay" in query_term.lower()
                else "Milvexian"
            )
            row = rows_by_nct.setdefault(
                nct_id,
                {
                    "nct_id": nct_id,
                    "query_term": query_term,
                    "representative_molecule": molecule,
                    "brief_title": ident.get("briefTitle") or "",
                    "overall_status": status.get("overallStatus") or "",
                    "phase": ";".join(design.get("phases") or []),
                    "start_date": (status.get("startDateStruct") or {}).get("date") or "",
                    "completion_date": (status.get("completionDateStruct") or {}).get("date")
                    or "",
                    "interventions": ";".join(sorted(set(interventions))),
                    "source_url": f"https://clinicaltrials.gov/study/{nct_id}",
                    "confidence": "high",
                    "notes": "ClinicalTrials.gov API v2 record; status should be rechecked before final clinical-stage claims",
                },
            )
            if query_term not in row["query_term"].split(";"):
                row["query_term"] += f";{query_term}"
            if molecule not in row["representative_molecule"].split(";"):
                row["representative_molecule"] += f";{molecule}"
            if interventions:
                merged = set(filter(None, row["interventions"].split(";")))
                merged.update(interventions)
                row["interventions"] = ";".join(sorted(merged))

    return sorted(rows_by_nct.values(), key=lambda row: row["nct_id"])


def main() -> None:
    molecule_rows, scaffold_rows = build_chembl_activity_tables()
    clinical_rows = build_clinical_representatives()
    pdb_rows = build_pdb_id_table()
    pdb_interaction_rows = build_pdb_interaction_table()
    clinical_trial_rows = build_clinical_trial_table()

    write_csv(OUT / "fxia_molecule_seed.csv", MOLECULE_COLUMNS, molecule_rows + clinical_rows)
    write_csv(OUT / "fxia_scaffold_seed.csv", SCAFFOLD_COLUMNS, scaffold_rows)
    write_csv(
        OUT / "fxia_pdb_ids_seed.csv",
        ["pdb_id", "source", "target", "source_url", "confidence", "notes"],
        pdb_rows,
    )
    write_csv(
        OUT / "fxia_pdb_interactions_seed.csv",
        PDB_INTERACTION_COLUMNS,
        pdb_interaction_rows,
    )
    write_csv(
        OUT / "fxia_clinical_trials_seed.csv",
        CLINICAL_TRIAL_COLUMNS,
        clinical_trial_rows,
    )

    summary = {
        "chembl_activity_rows": len(molecule_rows),
        "clinical_representative_rows": len(clinical_rows),
        "murcko_scaffold_clusters": len(scaffold_rows),
        "pdb_ids_from_chembl_target": len(pdb_rows),
        "pdb_nonpolymer_ligand_rows": len(pdb_interaction_rows),
        "clinical_trial_rows": len(clinical_trial_rows),
    }
    with (OUT / "fxia_seed_summary.json").open("w") as fh:
        json.dump(summary, fh, indent=2, sort_keys=True)


if __name__ == "__main__":
    main()
