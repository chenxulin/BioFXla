# BioFXIa

FXIa small-molecule scaffold discovery workflow.

## Pipeline

1. Collect and standardize public FXIa data.
2. Generate molecules with Mol2Mol, LibInvent, LinkInvent, and REINVENT de novo.
3. Filter generated molecules with RDKit structure standardization, Lipinski Rule-of-Five, clinical-candidate similarity, and patent-proxy similarity.
4. Merge filtered molecules into a de-duplicated candidate pool.
5. Dock candidate molecules with GNINA and rank by score/pose quality.
6. Advance selected compounds to MD, ADMET prediction, and FEP.

## Current Stage

Implemented:

- FXIa public seed table construction.
- REINVENT4 input preparation.
- REINVENT4 scale-up sampling helpers.
- RDKit property and Lipinski Rule-of-Five filtering.
- Multi-model candidate-pool merging.

Not yet implemented:

- Receptor preparation and GNINA docking workflow.
- MD setup/analysis.
- ADMET prediction workflow.
- FEP setup.

## Main Commands

Build processed seed tables:

```bash
mamba run -n aidd python scripts/data/build_fxia_seed_tables.py
```

Prepare REINVENT4 inputs:

```bash
mamba run -n aidd python scripts/generation/prepare_reinvent4_inputs.py --output-dir results/reinvent4_fxia_pilot
```

Run or resume REINVENT4 scale-up:

```bash
mamba run -n aidd bash scripts/generation/run_reinvent4_scaleup.sh
```

Filter one generated CSV:

```bash
mamba run -n aidd python scripts/filtering/filter_generated_molecules.py \
  --input results/reinvent4_fxia_pilot/raw_scale/mol2mol_scale.csv \
  --manifest results/reinvent4_fxia_pilot/generation_input_manifest.json \
  --output-dir results/reinvent4_fxia_pilot/filtered_scale/mol2mol_scale \
  --output-prefix mol2mol_scale
```

Merge filtered candidate pools:

```bash
mamba run -n aidd python scripts/filtering/merge_generation_candidates.py \
  --input mol2mol=results/reinvent4_fxia_pilot/filtered_scale/mol2mol_scale/mol2mol_scale_filtered.csv \
  --output results/candidates/fxia_generation_pool_filtered.csv
```

## Notes

- Clinical candidates are calibration/exclusion anchors only, not direct design templates.
- Patent-proxy filtering is technical triage, not a legal freedom-to-operate conclusion.
- Docking, MD, ADMET, and FEP outputs must remain clearly separated from source-backed activity data.
