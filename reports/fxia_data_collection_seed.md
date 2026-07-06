# FXIa Public Data Collection Seed

Generated: 2026-07-06

This is a first-pass, source-backed seed dataset for the FXIa scaffold atlas. It is
not yet a ranked patent-breakthrough scaffold recommendation. Patent, BindingDB,
full SAR extraction, and protein-ligand interaction annotation remain open.

## Collected Sources

- ChEMBL target and activity data: human coagulation factor XI / F11,
  `CHEMBL2820`, UniProt `P03951`.
  Source: https://www.ebi.ac.uk/chembl/target_report_card/CHEMBL2820/
- PubChem representative molecule identifiers:
  Asundexian PubChem CID `135206011`;
  Milvexian PubChem CID `118277544`.
- ClinicalTrials.gov API v2 records for `asundexian`, `BAY 2433334`,
  `milvexian`, and `BMS-986177`.
- RCSB PDB entry, non-polymer entity, and chemcomp metadata for PDB IDs linked
  from ChEMBL target `CHEMBL2820`.

## Output Tables

- `data/processed/fxia_molecule_seed.csv`
  - ChEMBL activity rows plus PubChem/ChEMBL clinical-candidate anchors.
  - Activity rows keep endpoint type, value, unit, relation, assay type,
    assay description, species, source URL, and RDKit-derived InChIKey fields.
- `data/processed/fxia_scaffold_seed.csv`
  - RDKit Bemis-Murcko scaffold clusters from ChEMBL activity structures.
  - Scaffold tiering only uses `Ki`, `IC50`, or `Kd` rows with nM units and
    `=`, `<`, or `<=` relations.
  - All scaffolds are currently `hold_for_more_data`; selectivity, PDB evidence,
    patent crowding, and synthesis have not been adjudicated.
- `data/processed/fxia_pdb_ids_seed.csv`
  - PDB IDs cross-referenced from ChEMBL target `CHEMBL2820`.
- `data/processed/fxia_pdb_interactions_seed.csv`
  - RCSB metadata-derived PDB ligand rows with ligand ID, ligand name, ligand
    SMILES when available, resolution, DOI, and source URL.
  - Binding-site contacts and S1/S2/S4 occupancy are marked
    `needs_verification`.
- `data/processed/fxia_clinical_trials_seed.csv`
  - ClinicalTrials.gov study status, phase, dates, and interventions.
- `data/processed/fxia_seed_summary.json`
  - Machine-readable counts for the current collection.

## Known Limits

- BindingDB has not yet been collected.
- Patent families and Markush/example proximity have not yet been collected.
- SAR papers have not yet been curated into structured series.
- RCSB interaction rows are metadata only; residue contacts require structure
  interaction extraction, e.g. PLIP, ProLIF, or curated literature review.
- ChEMBL document DOI fields are not yet expanded from document records.
- Clinical status should be rechecked before any final stage/status claim.

## Data Gaps

- `missing_activity_data`: BindingDB and literature tables not yet integrated.
- `missing_selectivity_data`: FXa, thrombin/FIIa, kallikrein, trypsin, plasmin,
  FIXa, and FVIIa counter-screen fields are empty pending extraction.
- `missing_pdb_evidence`: ligand contact annotations and pocket occupancy are
  `needs_verification`.
- `missing_patent_review`: no patent families have been scored; clinical
  chemotypes should be treated as crowded anchors only.
- `missing_synthesis_route`: synthetic accessibility is not yet assessed.
