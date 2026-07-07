# FXIa REINVENT4 Generation Package

This package is generated from `data/processed/fxia_molecule_seed.csv` and keeps clinical candidates as exclusion/calibration anchors only.

## Current Counts

- Mol2Mol non-clinical pilot seeds: 40
- LibInvent heuristic scaffold templates: 10
- LinkInvent BRICS-derived warhead pairs: 9
- Clinical similarity exclusion anchors: 2
- Patent-proxy similarity exclusions: 102

## Pilot Command

```bash
mamba run -n aidd reinvent -l results/reinvent4_fxia_pilot/mol2mol/mol2mol_pilot.log results/reinvent4_fxia_pilot/configs/mol2mol_pilot.toml
```

## Scale-Up Configs

- Mol2Mol: `configs/mol2mol_scale.toml` (100000 planned raw molecules)
- LibInvent: `configs/libinvent_scale.toml` (100000 planned raw molecules)
- LinkInvent: `configs/linkinvent_scale.toml` (90000 planned raw molecules)
- REINVENT de novo: `configs/reinvent_denovo_scale.toml` (100000 planned raw molecules)

Scale-up command:

```bash
mamba run -n aidd bash scripts/generation/run_reinvent4_scaleup.sh
```

This repository-level script starts or resumes all four paths, skips paths whose logs already contain `Finished REINVENT`, and writes PID/stdout/stderr files under `raw_scale/`. The generated `run_scaleup_commands.sh` is retained only as provenance for the initial package generation.

## Downstream Status

The scale-up outputs have been filtered and merged into:

```text
results/candidates/fxia_generation_pool_filtered.csv
```

That candidate pool has been docked with GNINA, reduced to 500 docking top hits, predicted with ADMET-AI, and ranked into 50 final MD candidates:

```text
results/admet/admet_ai/top_hits_ranked/top_hits_md_candidates.csv
results/admet/admet_ai/top_hits_ranked/top_hits_md_candidates.sdf
```

Do not treat the generation configs or generated molecules as medicinal chemistry recommendations. The retained molecules are computational hypotheses after property, similarity, docking, and ADMET-AI triage.
