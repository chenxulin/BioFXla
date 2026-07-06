# FXIa REINVENT4 Pilot Package

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

Do not treat the scale configs as medicinal chemistry recommendations. Run them only after the pilot output has been filtered for MW, TPSA, cLogP, charge, SA, clinical-candidate similarity, patent-proxy similarity, and reactive/PAINS liabilities.
